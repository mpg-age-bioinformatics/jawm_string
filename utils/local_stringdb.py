"""Versioned, shared STRING downloads and an indexed offline API adapter.

Only exact aliases are used. Ambiguous aliases choose the lexically first STRING
protein ID (after preferring an exact preferred name); this is deliberately
deterministic and is recorded in provenance, not claimed to reproduce API ranking.
"""
import csv
from contextlib import contextmanager, closing
import gzip
import hashlib
import json
import os
import re
from pathlib import Path
import shutil
import socket
import sqlite3
import tempfile
import time
from urllib.request import urlopen

VERSION = '12.5'
SCHEMA = 1
CHANNELS = dict(nscore='neighborhood', fscore='fusion', pscore='cooccurence',
                ascore='coexpression', escore='experimental', dscore='database',
                tscore='textmining')
MAPPING_POLICY = 'exact alias; prefer exact preferred name; lexical STRING ID tie-break'


def checksum(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


@contextmanager
def download_lock(folder, timeout=86400, poll=2):
    """Atomic cross-host lock. Never remove another process's/stale lock.

    A killed downloader may leave download.lock behind: inspect its PID/host and
    remove it manually only after confirming the owner is no longer running.
    """
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    lock = folder / 'download.lock'
    deadline = time.monotonic() + timeout
    announced = False
    while True:
        try:
            fd = os.open(lock, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
            break
        except FileExistsError:
            if not announced:
                print('Waiting for STRING download lock: ' + str(lock), flush=True)
                announced = True
            if time.monotonic() >= deadline:
                raise TimeoutError('STRING download lock still exists: ' + str(lock)
                                   + '; inspect owner before removing a stale lock.')
            time.sleep(poll)
    try:
        with os.fdopen(fd, 'w') as stream:
            json.dump(dict(pid=os.getpid(), host=socket.gethostname(), started=time.time()), stream)
        yield
    finally:
        lock.unlink()


def _ready(folder, species, network_type):
    stem = f'{species}.{network_type}'
    try:
        manifest = json.loads((folder / (stem + '.manifest.json')).read_text())
        if (manifest['schema'], manifest['version'], manifest['species'], manifest['network_type']) != (SCHEMA, VERSION, species, network_type):
            return None
        for name, entry in manifest['files'].items():
            if (folder / name).stat().st_size != entry['size']:
                return None
        database = folder / manifest['database']
        with closing(sqlite3.connect(database.as_uri() + '?mode=ro', uri=True)) as db:
            if db.execute('SELECT value FROM metadata WHERE key="schema"').fetchone() != (str(SCHEMA),):
                return None
        return manifest
    except (OSError, ValueError, KeyError, sqlite3.Error):
        return None


def _download(url, target):
    print('Downloading ' + url, flush=True)
    for attempt in range(3):
        try:
            with urlopen(url, timeout=180) as response, target.open('wb') as output:
                shutil.copyfileobj(response, output, 1024 * 1024)
                output.flush()
                os.fsync(output.fileno())
            # Read to EOF to verify the gzip trailer/CRC before publication.
            with gzip.open(target, 'rb') as stream:
                for _ in iter(lambda: stream.read(1024 * 1024), b''):
                    pass
            return
        except (OSError, EOFError):
            if attempt == 2:
                raise
            time.sleep(2 ** attempt)


def _build(database, files, species, network_type):
    prefix = str(species) + '.'
    with closing(sqlite3.connect(database)) as db, db:
        db.executescript('''
            PRAGMA journal_mode=OFF;
            CREATE TABLE metadata(key TEXT PRIMARY KEY, value TEXT);
            CREATE TABLE proteins(id TEXT PRIMARY KEY, name TEXT, annotation TEXT);
            CREATE TABLE aliases(alias TEXT, protein TEXT, PRIMARY KEY(alias, protein)) WITHOUT ROWID;
            CREATE TABLE edges(a TEXT, b TEXT, score REAL, nscore REAL, fscore REAL,
                pscore REAL, ascore REAL, escore REAL, dscore REAL, tscore REAL,
                PRIMARY KEY(a,b)) WITHOUT ROWID;
        ''')
        db.execute('INSERT INTO metadata VALUES (?,?)', ('schema', str(SCHEMA)))
        with gzip.open(files['protein.info'], 'rt') as stream:
            rows = csv.DictReader(stream, delimiter='\t')
            for row in rows:
                protein = row['#string_protein_id']
                if not protein.startswith(prefix):
                    raise ValueError('Incorrect species in STRING protein file.')
                db.execute('INSERT INTO proteins VALUES (?,?,?)', (protein, row['preferred_name'], row['annotation']))
                db.executemany('INSERT OR IGNORE INTO aliases VALUES (?,?)',
                               [(protein, protein), (row['preferred_name'], protein)])
        with gzip.open(files['protein.aliases'], 'rt') as stream:
            for row in csv.DictReader(stream, delimiter='\t'):
                protein = row['#string_protein_id']
                if not protein.startswith(prefix):
                    raise ValueError('Incorrect species in STRING aliases file.')
                db.execute('INSERT OR IGNORE INTO aliases VALUES (?,?)', (row['alias'], protein))
        with gzip.open(files['links'], 'rt') as stream:
            header = stream.readline().split()
            channels = CHANNELS.values() if network_type == 'functional' else ['experimental', 'database', 'textmining']
            required = {'protein1', 'protein2', 'combined_score', *channels}
            if not required.issubset(header):
                raise ValueError('Missing columns in STRING detailed links file.')
            for line in stream:
                values = line.split()
                if len(values) != len(header):
                    raise ValueError('Malformed STRING links row.')
                row = dict(zip(header, values))
                a, b = sorted((row['protein1'], row['protein2']))
                if not a.startswith(prefix) or not b.startswith(prefix):
                    raise ValueError('Incorrect species in STRING links file.')
                scores = [int(row.get(key, 0)) for key in ['combined_score', *CHANNELS.values()]]
                if any(value < 0 or value > 1000 for value in scores):
                    raise ValueError('Invalid STRING score.')
                db.execute('INSERT OR IGNORE INTO edges VALUES (?,?,?,?,?,?,?,?,?,?)',
                           (a, b, *(value / 1000 for value in scores)))
        db.execute('CREATE INDEX edges_b ON edges(b)')
        if not db.execute('SELECT count(*) FROM proteins').fetchone()[0]:
            raise ValueError('Empty STRING protein database.')
        if db.execute('PRAGMA quick_check').fetchone() != ('ok',):
            raise ValueError('Invalid STRING SQLite database.')


def ensure_database(root, species, network_type='functional', timeout=86400, poll=2):
    if type(species) is not int or species <= 0:
        raise ValueError('A positive integer STRING species ID is required.')
    if network_type not in ('functional', 'physical'):
        raise ValueError('Unsupported STRING network type.')
    folder = Path(root).expanduser().resolve() / VERSION
    # Always take the lock, even when files appear ready: a different workflow
    # may currently be downloading or publishing this version.
    with download_lock(folder, timeout=timeout, poll=poll):
        ready = _ready(folder, species, network_type)
        if ready:
            return folder, ready
        kinds = {'protein.info': 'protein.info', 'protein.aliases': 'protein.aliases',
                 'links': 'protein.links.detailed' if network_type == 'functional'
                 else 'protein.physical.links.detailed'}
        with tempfile.TemporaryDirectory(prefix='.download-', dir=folder) as scratch:
            scratch = Path(scratch)
            files, entries = {}, {}
            for role, kind in kinds.items():
                name = f'{species}.{kind}.v{VERSION}.txt.gz'
                url = f'https://stringdb-downloads.org/download/{kind}.v{VERSION}/{name}'
                path = scratch / name
                _download(url, path)
                files[role] = path
                entries[name] = dict(url=url, sha256=checksum(path), size=path.stat().st_size)
            stem = f'{species}.{network_type}'
            name = stem + '.sqlite'
            database = scratch / name
            print('Indexing STRING database: ' + str(folder / name), flush=True)
            _build(database, files, species, network_type)
            entries[name] = dict(sha256=checksum(database), size=database.stat().st_size)
            manifest = dict(schema=SCHEMA, version=VERSION, species=species,
                            network_type=network_type, database=name, files=entries,
                            mapping_policy=MAPPING_POLICY)
            marker = scratch / (stem + '.manifest.json')
            marker.write_text(json.dumps(manifest, indent=2) + '\n')
            # Publish the readiness manifest last, while holding the lock.
            for path in [*files.values(), database, marker]:
                os.replace(path, folder / path.name)
            if not _ready(folder, species, network_type):
                raise ValueError('Published STRING database failed validation.')
            return folder, manifest


class LocalStringClient:
    def __init__(self, root, species, progress=print, network_type='functional'):
        self.species, self.progress, self.keys = species, progress, []
        self.folder, self.manifest = ensure_database(root, species, network_type)
        self.database = self.folder / self.manifest['database']
        self.network_type = network_type

    def request(self, method, identifiers, **parameters):
        records = []
        # Published indexes are never edited in place; atomic replacement lets
        # readers use immutable snapshots without per-query network-FS locks.
        with closing(sqlite3.connect(self.database.as_uri() + '?mode=ro&immutable=1', uri=True)) as db:
            db.execute('PRAGMA cache_size=-65536')
            if method == 'get_string_ids':
                for index, identifier in enumerate(identifiers):
                    row = db.execute('''SELECT p.id,p.name,p.annotation FROM aliases a
                        JOIN proteins p ON p.id=a.protein WHERE a.alias=?
                        ORDER BY (p.name=?) DESC,p.id LIMIT 1''', (identifier, identifier)).fetchone()
                    if row:
                        records.append(dict(queryIndex=index, queryItem=identifier,
                                            stringId=row[0], preferredName=row[1],
                                            annotation=row[2], ncbiTaxonId=self.species))
            elif method == 'interaction_partners':
                if parameters.get('network_type', 'functional') != self.network_type:
                    raise ValueError('Local STRING database network type mismatch.')
                seen = set()
                threshold = float(parameters.get('required_score', 0)) / 1000
                for identifier in identifiers:
                    for row in db.execute('''SELECT e.*,p.name,q.name FROM edges e
                        JOIN proteins p ON p.id=e.a JOIN proteins q ON q.id=e.b
                        WHERE (e.a=? OR e.b=?) AND e.score>=?''', (identifier, identifier, threshold)):
                        if row[:2] in seen:
                            continue
                        seen.add(row[:2])
                        records.append(dict(zip(['stringId_A','stringId_B','score',*CHANNELS,
                                                 'preferredName_A','preferredName_B'], row)))
                        if len(records) > 200000:
                            raise ValueError('Too many STRING edges; raise confidence or reduce seeds.')
            else:
                raise ValueError('Unknown STRING method.')
        return records


def resolve_species(taxon, organism, reference):
    if str(taxon).strip():
        if not str(taxon).strip().isdigit() or int(taxon) <= 0:
            raise ValueError('Invalid STRING taxonomy ID.')
        return int(taxon)
    normalize = lambda value: re.sub(r'[_\s]+', ' ', value.strip().lower())
    matches = set()
    with gzip.open(reference, 'rt') as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            if normalize(organism) in {normalize(row['STRING_name_compact']), normalize(row['official_name_NCBI'])}:
                matches.add(int(row['#taxon_id']))
    if len(matches) != 1:
        raise ValueError('Set an unambiguous organism or NCBI taxonomy ID for STRING.')
    return matches.pop()
