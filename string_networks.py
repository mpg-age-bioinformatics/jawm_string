import jawm
import logging
import subprocess
import uuid
import os
from pathlib import Path

logger = logging.getLogger("jawm_string")

ROOT = Path(__file__).resolve().parent
STRING_IMAGE = "mpgagebioinformatics/networkx:3.5"


string_network = jawm.Process(
    name="string_network",
    script=r'''#!/usr/bin/env python3
"""STRING 12.5 PPI workflow adapted from ppis.ipynb; never executes notebook cells."""

import csv
import gzip
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import re
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen, build_opener, HTTPRedirectHandler
import networkx as nx
import numpy as np
import pandas as pd

URL = 'https://version-12-5.string-db.org/api'
NOTEBOOK_SHA256 = '9594735cf295863aa4da0fc5aa4d7276acc83dc34cad6347215123dbe2bcc632'
CHANNELS = ['nscore', 'fscore', 'pscore', 'ascore', 'escore', 'dscore', 'tscore']
EDGE_COLUMNS = [
    'stringId_A',
    'stringId_B',
    'preferredName_A',
    'preferredName_B',
    'score',
    *CHANNELS,
]
DEFAULTS = {
    'padj': 0.05,
    'min_abs_log2fc': 0.0,
    'top_n': 1000,
    'rank_mode': 'absolute',
    'expansion_fraction': 0.25,
    'network_type': 'functional',
    'min_combined_score': 0.7,
    'evidence_min': {},
    'evidence_mode': 'any',
}
MAX_EDGES = 200000


@lru_cache(maxsize=1)
def species_names():
    names = {}
    with gzip.open(REFERENCE_DIR / 'species.v12.5.tsv.gz', 'rt') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            for key in ['STRING_name_compact', 'official_name_NCBI']:
                name = re.sub('[_\\s]+', ' ', row[key].strip().lower())
                names.setdefault(name, set()).add(int(row['#taxon_id']))
    return {name: next(iter(ids)) for name, ids in names.items() if len(ids) == 1}


def species_hint(metadata):
    value = metadata.get('taxon_id') or metadata.get('ncbi_taxon_id')
    if (
        value is not None
        and (not isinstance(value, bool))
        and re.fullmatch('[0-9]+', str(value).strip())
    ):
        try:
            identifier = int(value)
            if identifier > 0:
                return identifier
        except (ValueError, TypeError):
            pass
    return species_names().get(
        re.sub('[_\\s]+', ' ', str(metadata.get('organism') or '').strip().lower())
    )


def validate(settings, metadata):
    if not isinstance(settings, dict) or set(settings) - set(DEFAULTS) - {'species'}:
        raise ValueError('Invalid STRING settings.')
    result = dict(DEFAULTS, **settings)
    if type(result['top_n']) is not int or not 1 <= result['top_n'] <= 2000:
        raise ValueError('Select between 1 and 2000 seed genes.')
    for key, minimum, maximum in [
        ('padj', 1e-300, 1),
        ('min_abs_log2fc', 0, 1000),
        ('expansion_fraction', 0, 1),
        ('min_combined_score', 0, 1),
    ]:
        value = result[key]
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or (not math.isfinite(value))
            or (not minimum <= value <= maximum)
        ):
            raise ValueError('Invalid ' + key + '.')
    if (
        result['rank_mode'] not in ['absolute', 'signed', 'log2fc', 'padj']
        or result['network_type'] not in ['functional', 'physical']
        or result['evidence_mode'] not in ['any', 'all']
    ):
        raise ValueError('Invalid STRING network or ranking mode.')
    if not isinstance(result['evidence_min'], dict) or set(
        result['evidence_min']
    ) - set(CHANNELS):
        raise ValueError('Unknown evidence channel.')
    for value in result['evidence_min'].values():
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or (not math.isfinite(value))
            or (not 0 < value <= 1)
        ):
            raise ValueError('Evidence thresholds must be between zero and one.')
    hint = species_hint(metadata)
    if result.get('species') and result['species'] != hint:
        raise ValueError('The chosen taxonomy ID does not match the dataset organism.')
    result['species'] = hint
    if type(result['species']) is not int or not 1 <= result['species'] <= 10**9:
        raise ValueError(
            'Set the dataset organism’s NCBI taxonomy ID before requesting STRING.'
        )
    if result['top_n'] * (1 + result['expansion_fraction']) > 3000:
        raise ValueError('Choose settings with at most 3000 seed plus added proteins.')
    return result


def ranking_score(fc, padj, mode):
    significance = -np.log10(np.clip(padj, 1e-300, 1))
    if mode == 'log2fc':
        return fc
    if mode == 'padj':
        return significance
    score = fc * significance
    return np.abs(score) if mode == 'absolute' else score


def rank_genes(frame, spec):
    data = frame.copy()
    if not {'ensembl_gene_id', 'log2FoldChange', 'padj'}.issubset(data):
        raise ValueError(
            'The contrast lacks gene IDs, fold changes or adjusted p-values.'
        )
    data['ensembl_gene_id'] = data['ensembl_gene_id'].astype('string')
    for column in ['log2FoldChange', 'padj']:
        data[column] = pd.to_numeric(data[column], errors='coerce')
    valid = (
        data['ensembl_gene_id'].notna()
        & data['ensembl_gene_id'].str.len().gt(0)
        & np.isfinite(data.log2FoldChange)
        & data.padj.between(0, 1)
    )
    selected = data.loc[
        valid
        & (data.padj < spec['padj'])
        & (data.log2FoldChange.abs() >= spec['min_abs_log2fc'])
    ].copy()
    selected['signed_score'] = selected.log2FoldChange * -np.log10(
        selected.padj.clip(lower=1e-300)
    )
    selected['ranking_score'] = ranking_score(
        selected.log2FoldChange, selected.padj, spec['rank_mode']
    )
    selected = selected.sort_values(
        ['ranking_score', 'padj', 'ensembl_gene_id'],
        ascending=[False, True, True],
        kind='stable',
    ).drop_duplicates('ensembl_gene_id')
    eligible = len(selected)
    selected = selected.head(spec['top_n']).copy()
    selected['rank'] = np.arange(1, len(selected) + 1)
    selected['direction'] = np.where(
        selected.log2FoldChange > 0,
        'up',
        np.where(selected.log2FoldChange < 0, 'down', 'unchanged'),
    )
    return (
        data,
        selected,
        {
            'input_rows': len(data),
            'invalid_rows': int((~valid).sum()),
            'de_genes': eligible,
            'selected_genes': len(selected),
        },
    )


class NoRedirect(HTTPRedirectHandler):

    def redirect_request(self, *args, **kwargs):
        return None


class StringClient:

    def __init__(self, cache, species, progress=lambda message: None):
        self.cache = Path(cache)
        self.cache.mkdir(parents=True, exist_ok=True)
        self.species = species
        self.progress = progress
        self.keys = []
        self.opener = build_opener(NoRedirect())

    def request(self, method, identifiers, **parameters):
        if not identifiers:
            return []
        if method not in ['get_string_ids', 'interaction_partners']:
            raise ValueError('Unknown STRING method.')
        payload = {
            'identifiers': '\r'.join(identifiers),
            'species': self.species,
            'caller_identity': 'bioinformatics.studio.rnaseq',
            **parameters,
        }
        key = hashlib.sha256(
            json.dumps([URL, method, payload], sort_keys=True).encode()
        ).hexdigest()
        self.keys.append(key)
        path = self.cache / (key + '.json')
        if path.is_file():
            return json.loads(path.read_text())
        for attempt in range(4):
            time.sleep(1 if attempt == 0 else 2**attempt)
            self.progress('Querying STRING ' + method)
            try:
                request = Request(
                    URL + '/json/' + method,
                    data=urlencode(payload).encode(),
                    headers={'Content-Type': 'application/x-www-form-urlencoded'},
                )
                with self.opener.open(request, timeout=180) as response:
                    content = response.read(64 * 1024 * 1024 + 1)
                    if len(content) > 64 * 1024 * 1024:
                        raise ValueError(
                            'STRING response exceeds the safe network size. Increase the confidence threshold or reduce seeds.'
                        )
                    value = json.loads(content)
                if not isinstance(value, list):
                    raise ValueError('Unexpected STRING response.')
                temporary = path.with_suffix('.tmp')
                temporary.write_text(json.dumps(value))
                temporary.replace(path)
                return value
            except HTTPError as error:
                if error.code == 404 and method == 'get_string_ids':
                    try:
                        body = json.loads(error.read(65536))
                        errors = body if isinstance(body, list) else [body]
                        if errors and all(isinstance(item, dict) and item.get('Error') == 'nothing found' for item in errors):
                            return []
                    except (ValueError, OSError):
                        pass
                if error.code not in [429, 500, 502, 503, 504] or attempt == 3:
                    raise ValueError(
                        'STRING request failed (HTTP ' + str(error.code) + ').'
                    ) from None
            except (URLError, TimeoutError, OSError):
                if attempt == 3:
                    raise ValueError(
                        'STRING is unavailable. Retry this network later.'
                    ) from None
        raise ValueError('STRING request failed.')


def map_ids(frame, client):
    ids = sorted(set((str(g) for g in frame.ensembl_gene_id.dropna() if str(g))))
    if any(
        (
            len(gene) > 500 or any((character in gene for character in '\r\n\t'))
            for gene in ids
        )
    ):
        raise ValueError('Invalid stable gene identifiers in the contrast.')
    if len(ids) > 200000:
        raise ValueError('Too many gene identifiers for an interactive STRING run.')
    records = []
    seen = set()
    for start in range(0, len(ids), 500):
        batch = ids[start : start + 500]
        client.progress(
            f'Mapping gene IDs {start + 1}–{min(start + 500, len(ids))} of {len(ids)}'
        )
        for record in client.request('get_string_ids', batch, echo_query=1):
            try:
                if not re.fullmatch(
                    '[0-9]+', str(record['queryIndex'])
                ) or not re.fullmatch('[0-9]+', str(record['ncbiTaxonId'])):
                    raise ValueError()
                index = int(record['queryIndex'])
                taxon = int(record['ncbiTaxonId'])
                protein = record['stringId']
                if (
                    not 0 <= index < len(batch)
                    or taxon != client.species
                    or (not str(protein).startswith(str(client.species) + '.'))
                ):
                    raise ValueError()
                gene = batch[index]
                if gene in seen:
                    raise ValueError('Ambiguous STRING mapping for ' + gene + '.')
                seen.add(gene)
                records.append(
                    {
                        'ensembl_gene_id': gene,
                        'stringId': str(protein),
                        'preferredName': str(record.get('preferredName') or protein),
                        'annotation': str(record.get('annotation') or ''),
                    }
                )
            except (KeyError, TypeError, ValueError) as error:
                raise ValueError(
                    'Invalid or ambiguous STRING species/mapping response.'
                ) from error
    return pd.DataFrame(
        records, columns=['ensembl_gene_id', 'stringId', 'preferredName', 'annotation']
    )


def filter_edges(records, spec):
    if not records:
        return pd.DataFrame(columns=EDGE_COLUMNS)
    frame = pd.DataFrame(records)
    if not set(EDGE_COLUMNS).issubset(frame):
        raise ValueError('STRING evidence fields are missing.')
    frame = frame[EDGE_COLUMNS].copy()
    for suffix in ['A', 'B']:
        frame['preferredName_' + suffix] = (
            frame['preferredName_' + suffix]
            .fillna(frame['stringId_' + suffix])
            .astype(str)
        )
    for column in ['score', *CHANNELS]:
        frame[column] = pd.to_numeric(frame[column], errors='coerce')
        if not frame[column].between(0, 1).all():
            raise ValueError('Invalid STRING evidence score scale.')
    for column in ['stringId_A', 'stringId_B']:
        if (
            not frame[column]
            .map(
                lambda value: isinstance(value, str)
                and value.startswith(str(spec['species']) + '.')
            )
            .all()
        ):
            raise ValueError('STRING returned an interaction from another species.')
    keep = (frame.score >= spec['min_combined_score']) & (
        frame.stringId_A != frame.stringId_B
    )
    if spec['evidence_min']:
        gates = pd.DataFrame(
            {name: frame[name] >= value for name, value in spec['evidence_min'].items()}
        )
        keep &= (
            gates.any(axis=1) if spec['evidence_mode'] == 'any' else gates.all(axis=1)
        )
    frame = frame.loc[keep].copy()
    frame['pair'] = [
        tuple(sorted(pair)) for pair in zip(frame.stringId_A, frame.stringId_B)
    ]
    return (
        frame.sort_values('score', ascending=False, kind='stable')
        .drop_duplicates('pair')
        .drop(columns='pair')
        .reset_index(drop=True)
    )


def partner_edges(ids, client, spec):
    records = []
    for start in range(0, len(ids), 500):
        records.extend(
            client.request(
                'interaction_partners',
                ids[start : start + 500],
                required_score=math.ceil(spec['min_combined_score'] * 1000),
                network_type=spec['network_type'],
            )
        )
        if len(records) > MAX_EDGES:
            raise ValueError(
                'Too many STRING edges; raise the confidence threshold or reduce seed genes.'
            )
    return filter_edges(records, spec)


def save_table(frame, path):
    frame.to_csv(path, sep='\t', index=False)


def build_network(frame, spec, folder, client, sample_columns=(), provenance=None):
    folder = Path(folder)
    data, top, summary = rank_genes(frame, spec)
    save_table(top, folder / 'selected_genes.tsv')
    mapping = (
        map_ids(data, client)
        if len(top)
        else pd.DataFrame(
            columns=['ensembl_gene_id', 'stringId', 'preferredName', 'annotation']
        )
    )
    save_table(mapping, folder / 'string_mapping.tsv')
    mapped = top.merge(mapping, on='ensembl_gene_id', how='left', validate='one_to_one')
    mapped['mapping_status'] = np.where(mapped.stringId.notna(), 'mapped', 'unmapped')
    save_table(mapped, folder / 'mapping.tsv')
    seeds = (
        mapped.dropna(subset=['stringId'])
        .sort_values('rank')
        .drop_duplicates('stringId')
    )
    seed_ids = set(seeds.stringId)
    budget = math.floor(len(seed_ids) * spec['expansion_fraction'])
    edges = partner_edges(sorted(seed_ids), client, spec)
    names = dict(zip(mapping.stringId, mapping.preferredName))
    links = []
    for row in edges.to_dict('records'):
        names[row['stringId_A']] = row['preferredName_A']
        names[row['stringId_B']] = row['preferredName_B']
        for seed, candidate in [
            (row['stringId_A'], row['stringId_B']),
            (row['stringId_B'], row['stringId_A']),
        ]:
            if seed in seed_ids and candidate not in seed_ids:
                links.append((candidate, seed, row['score']))
    candidates = (
        pd.DataFrame(links, columns=['stringId', 'seed', 'score'])
        .drop_duplicates(['stringId', 'seed'])
        .groupby('stringId', as_index=False)
        .agg(seed_connections=('seed', 'nunique'), confidence_sum=('score', 'sum'))
    )
    candidates = candidates.sort_values(
        ['seed_connections', 'confidence_sum', 'stringId'],
        ascending=[False, False, True],
        kind='stable',
    )
    candidates['expansion_rank'] = np.arange(1, len(candidates) + 1)
    candidates['selected'] = candidates.expansion_rank <= budget
    candidates['preferredName'] = candidates.stringId.map(names)
    save_table(candidates, folder / 'expansion_candidates.tsv')
    added_ids = set(candidates.loc[candidates.selected, 'stringId'])
    final_ids = seed_ids | added_ids
    additional = partner_edges(sorted(added_ids), client, spec)
    combined = filter_edges(
        edges.to_dict('records') + additional.to_dict('records'), spec
    )
    final_edges = combined.loc[
        combined.stringId_A.isin(final_ids) & combined.stringId_B.isin(final_ids)
    ].copy()
    final_edges['interaction'] = spec['network_type']
    final_edges['edge_role'] = [
        (
            'seed-seed'
            if a in seed_ids and b in seed_ids
            else 'added-added' if a in added_ids and b in added_ids else 'seed-added'
        )
        for a, b in zip(final_edges.stringId_A, final_edges.stringId_B)
    ]
    graph = nx.Graph(
        species=spec['species'],
        network_type=spec['network_type'],
        settings_json=json.dumps(spec, sort_keys=True),
        **provenance or {},
    )
    for row in seeds.to_dict('records'):
        graph.add_node(
            row['stringId'],
            name=str(row['preferredName']),
            role='seed',
            rank=int(row['rank']),
        )
    for row in candidates.loc[candidates.selected].to_dict('records'):
        graph.add_node(
            row['stringId'],
            name=str(row['preferredName']),
            role='added',
            confidence_sum=float(row['confidence_sum']),
        )
    expression = data.merge(
        mapping, on='ensembl_gene_id', how='inner', validate='many_to_one'
    )
    expression = expression.loc[expression.stringId.isin(final_ids)].copy()
    expression = expression.merge(
        top[['ensembl_gene_id', 'rank']], on='ensembl_gene_id', how='left'
    )
    expression['has_log2FC'] = np.isfinite(expression.log2FoldChange)
    expression = expression.sort_values(
        ['rank', 'has_log2FC', 'padj', 'ensembl_gene_id'],
        ascending=[True, False, True, True],
        na_position='last',
        kind='stable',
    )
    expression['representative'] = ~expression.stringId.duplicated()
    gene_groups = expression.groupby('stringId').ensembl_gene_id.agg(
        lambda values: ';'.join(dict.fromkeys(values))
    )
    numeric_columns = [
        name
        for name in [
            'baseMean',
            'log2FoldChange',
            'lfcSE',
            'pvalue',
            'padj',
            *sample_columns,
        ]
        if name in data
    ]
    for node in graph:
        graph.nodes[node].update(
            expression_status='not_in_results', direction='unknown'
        )
    for row in expression.loc[expression.representative].to_dict('records'):
        attrs = graph.nodes[row['stringId']]
        attrs.update(
            expression_gene_id=str(row['ensembl_gene_id']),
            gene_ids=str(gene_groups[row['stringId']]),
            expression_status='available' if row['has_log2FC'] else 'missing_log2FC',
        )
        if isinstance(row.get('gene_name'), str):
            attrs['gene_name'] = row['gene_name']
        for column in numeric_columns:
            value = pd.to_numeric(row[column], errors='coerce')
            if pd.notna(value) and np.isfinite(value):
                attrs[column] = float(value)
        if 'log2FoldChange' in attrs:
            fc = attrs['log2FoldChange']
            attrs['direction'] = 'up' if fc > 0 else 'down' if fc < 0 else 'unchanged'
            if 'padj' in attrs and 0 <= attrs['padj'] <= 1:
                attrs['is_DE'] = bool(
                    attrs['padj'] < spec['padj'] and abs(fc) >= spec['min_abs_log2fc']
                )
                attrs['signed_score'] = float(
                    fc * -np.log10(max(attrs['padj'], 1e-300))
                )
                attrs['ranking_score'] = float(
                    ranking_score(fc, attrs['padj'], spec['rank_mode'])
                )
    for row in final_edges.to_dict('records'):
        a = row.pop('stringId_A')
        b = row.pop('stringId_B')
        graph.add_edge(a, b, **row)
    isolated = list(nx.isolates(graph))
    save_table(
        pd.DataFrame(
            [{'stringId': node, **graph.nodes[node]} for node in isolated],
            columns=None if isolated else ['stringId', 'name', 'role'],
        ),
        folder / 'removed_isolated_nodes.tsv',
    )
    graph.remove_nodes_from(isolated)
    for node in graph:
        graph.nodes[node].update(
            degree=int(graph.degree(node)),
            seed_connections=sum(
                (neighbor in seed_ids for neighbor in graph.neighbors(node))
            ),
        )
    client.progress('Laying out network')
    positions = (
        nx.forceatlas2_layout(graph, max_iter=60, seed=42)
        if len(graph) > 1
        else {node: np.zeros(2) for node in graph}
    )
    nodes = []
    for node, attrs in graph.nodes(data=True):
        nodes.append(
            {
                'stringId': node,
                **attrs,
                'x': float(positions[node][0]),
                'y': float(positions[node][1]),
            }
        )
    scores = [row['signed_score'] for row in nodes if 'signed_score' in row]
    limits = np.percentile(scores, [5, 95]).tolist() if scores else [0.0, 0.0]
    save_table(
        pd.DataFrame(
            nodes, columns=None if nodes else ['stringId', 'name', 'role', 'degree']
        ),
        folder / 'nodes.tsv',
    )
    save_table(final_edges, folder / 'edges.tsv')
    save_table(
        expression.loc[expression.stringId.isin(graph)],
        folder / 'node_expression_mapping.tsv',
    )
    nx.write_graphml(graph, folder / 'network.graphml')
    summary.update(
        unmapped_genes=int(mapped.stringId.isna().sum()),
        seed_proteins=len(seed_ids),
        expansion_budget=budget,
        added_proteins=len(added_ids),
        retained_seed_proteins=len(seed_ids & set(graph)),
        added_with_log2FC=sum(
            ('log2FoldChange' in graph.nodes[node] for node in added_ids & set(graph))
        ),
        removed_isolates=len(isolated),
        nodes=len(graph),
        edges=graph.number_of_edges(),
        status=(
            'no_DE_genes'
            if top.empty
            else 'no_mapped_seeds' if not seed_ids else 'complete'
        ),
    )
    value = {
        'nodes': nodes,
        'edges': final_edges.to_dict('records'),
        'summary': summary,
        'color_limits': limits,
        'settings': spec,
    }
    encoded = json.dumps(value, allow_nan=False)
    if len(encoded.encode()) > 64 * 1024 * 1024:
        raise ValueError(
            'The network is too large for the interactive viewer. Raise confidence or reduce seeds.'
        )
    (folder / 'network.json').write_text(encoded)
    (folder / 'summary.json').write_text(json.dumps(summary, indent=2))
    return summary


'Shared POSIX disk cache for STRING 12.5, including overlapping query batches.\n\nA single advisory lock serializes cache misses, HTTP calls, retries and publication\nacross all jobs sharing this directory. The upstream client waits at least one\nsecond before every uncached HTTP attempt. No lock is held outside API/cache work.\nThe filesystem must support cross-host flock (also required on HPC shared storage).\n'
from contextlib import contextmanager
import fcntl
import hashlib
import json
import os
from pathlib import Path
import tempfile


@contextmanager
def locked(path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('a') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, name = tempfile.mkstemp(dir=path.parent, prefix='.', suffix='.tmp')
    try:
        with os.fdopen(handle, 'w') as stream:
            json.dump(value, stream, allow_nan=False)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(name, path)
    finally:
        if os.path.exists(name):
            os.unlink(name)


class CachedStringClient(StringClient):

    def __init__(self, cache, species, progress=print):
        super().__init__(Path(cache) / 'requests', species, progress)
        self.shared = Path(cache)
        self.hits = 0
        self.misses = 0

    def entity_path(self, method, identifier, parameters):
        key = hashlib.sha256(
            json.dumps(
                [URL, self.species, method, identifier, parameters], sort_keys=True
            ).encode()
        ).hexdigest()
        return self.shared / 'entities' / key[:2] / (key + '.json')

    def read_entity(self, path):
        try:
            if path.stat().st_size > 64 * 1024 * 1024:
                return None
            value = json.loads(path.read_text())
            return value if isinstance(value, list) else None
        except (OSError, ValueError):
            return None

    def request(self, method, identifiers, **parameters):
        if not identifiers:
            return []
        if method not in ['get_string_ids', 'interaction_partners']:
            raise ValueError('Unknown STRING method.')
        if len(set(identifiers)) != len(identifiers):
            raise ValueError('STRING queries must contain distinct identifiers.')
        with locked(self.shared / 'https.lock'):
            entries = {
                gene: self.read_entity(self.entity_path(method, gene, parameters))
                for gene in identifiers
            }
            missing = [gene for gene, value in entries.items() if value is None]
            self.hits += len(identifiers) - len(missing)
            self.misses += len(missing)
            if missing:
                # Canonical batches reuse exact requests after interrupted jobs.
                missing.sort()
                for start in range(0, len(missing), 500):
                    batch = missing[start : start + 500]
                    records = super().request(method, batch, **parameters)
                    if len(records) > MAX_EDGES:
                        raise ValueError(
                            'Too many STRING response records; raise the confidence threshold.'
                        )
                    grouped = {gene: [] for gene in batch}
                    if method == 'get_string_ids':
                        seen = set()
                        for record in records:
                            try:
                                index = int(record['queryIndex'])
                                taxon = int(record['ncbiTaxonId'])
                                if (
                                    str(index) != str(record['queryIndex'])
                                    or not 0 <= index < len(batch)
                                    or index in seen
                                    or (taxon != self.species)
                                    or (
                                        not record['stringId'].startswith(
                                            str(self.species) + '.'
                                        )
                                    )
                                ):
                                    raise ValueError()
                                seen.add(index)
                                grouped[batch[index]].append(dict(record, queryIndex=0))
                            except (KeyError, ValueError, TypeError, AttributeError):
                                raise ValueError(
                                    'Invalid or ambiguous STRING mapping response.'
                                ) from None
                    else:
                        # Validate records before publishing them to the shared cache.
                        filter_edges(
                            records,
                            dict(
                                DEFAULTS,
                                species=self.species,
                                min_combined_score=0,
                                evidence_min={},
                            ),
                        )
                        for record in records:
                            endpoints = {
                                record['stringId_A'],
                                record['stringId_B'],
                            } & grouped.keys()
                            if not endpoints:
                                raise ValueError(
                                    'STRING returned an interaction unrelated to the query.'
                                )
                            for gene in endpoints:
                                grouped[gene].append(record)
                    for gene, value in grouped.items():
                        atomic_json(self.entity_path(method, gene, parameters), value)
                        entries[gene] = value
            result = []
            seen = set()
            for index, gene in enumerate(identifiers):
                for record in entries[gene]:
                    if method == 'get_string_ids':
                        result.append(dict(record, queryIndex=index))
                    else:
                        key = json.dumps(record, sort_keys=True)
                        if key not in seen:
                            seen.add(key)
                            result.append(record)
                if len(result) > MAX_EDGES:
                    raise ValueError(
                        'Too many STRING records; raise confidence or reduce seeds.'
                    )
            return result


"Run the app's immutable scientific workflow with an HPC-safe shared cache."
import hashlib
import json
from pathlib import Path
import platform
import shutil
import tempfile
import time
import networkx
import numpy
import pandas as pd

REFERENCE_DIR = Path(__file__).resolve().parent / 'reference'


def checksum(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        while data := handle.read(1024 * 1024):
            digest.update(data)
    return digest.hexdigest()


def run(config, client_factory=CachedStringClient):
    spec = validate({}, config['metadata'])
    local_client = None
    database_source = {'backend': 'api', 'url': URL}
    if config.get('STRINGDB'):
        import sys
        sys.path.insert(0, str(Path(config['module_file']).parent))
        from utils.local_stringdb import LocalStringClient
        local_client = LocalStringClient(config['STRINGDB'], spec['species'],
                                         network_type=spec['network_type'])
        database_source = {'backend': 'local', 'path': str(local_client.folder),
                           'manifest': local_client.manifest,
                           'adapter_sha256': checksum(Path(config['module_file']).parent / 'utils' / 'local_stringdb.py')}
    source = Path(config['input'])
    samples = Path(config['samples'])
    output = Path(config['output'])
    prefix = config['contrast'] + '_'
    signatures = {
        name: checksum(path)
        for name, path in {
            'input': source,
            'samples': samples,
            'module': Path(config['module_file']),
            'species': REFERENCE_DIR / 'species.v12.5.tsv.gz',
        }.items()
    }
    identity = {
        'database_source': database_source,
        'files': signatures,
        'settings': spec,
        'metadata': config['metadata'],
        'pipeline_revision': config['pipeline_revision'],
        'module_revision': config.get('module_revision', 'unrecorded'),
        'container': config['container'],
        'packages': {
            'networkx': networkx.__version__,
            'numpy': numpy.__version__,
            'pandas': pd.__version__,
        },
        'python': platform.python_version(),
    }
    fingerprint = hashlib.sha256(
        json.dumps(identity, sort_keys=True).encode()
    ).hexdigest()
    output.parent.mkdir(parents=True, exist_ok=True)
    with locked(output.parent / ('.' + output.name + '.lock')):
        previous = output.parent / ('.' + output.name + '.previous')
        if previous.exists() and (not output.exists()):
            previous.rename(output)
        try:
            old = json.loads((output / (prefix + 'manifest.json')).read_text())
            artifacts = old.get('artifacts_sha256', {})
            if (
                old.get('request_fingerprint') == fingerprint
                and artifacts
                and all(
                    (
                        (output / name).is_file() and checksum(output / name) == digest
                        for name, digest in artifacts.items()
                    )
                )
            ):
                print(
                    'Completed STRING network unchanged; reusing ' + str(output),
                    flush=True,
                )
                return old
        except (ValueError, OSError):
            pass
        temporary = Path(
            tempfile.mkdtemp(prefix='.' + output.name + '.', dir=output.parent)
        )
        try:
            created = time.time()
            shutil.copyfile(source, temporary / 'inputs.tsv')
            shutil.copyfile(samples, temporary / 'samples.tsv')
            frame = pd.read_csv(
                temporary / 'inputs.tsv',
                sep='\t',
                comment='#',
                dtype={'ensembl_gene_id': str},
            )
            annotations = pd.read_csv(
                temporary / 'samples.tsv', sep='\t', dtype=str, keep_default_na=False
            )
            annotations = annotations.rename(columns={annotations.columns[0]: 'sample'})
            columns = [sample for sample in annotations['sample'] if sample in frame]
            annotations = annotations[annotations['sample'].isin(columns)]
            client = local_client or client_factory(config['cache'], spec['species'])
            _, _, selection = rank_genes(frame, spec)
            summary = build_network(
                frame,
                spec,
                temporary,
                client,
                columns,
                {
                    'network_name': config['name'],
                    'dataset': config['dataset'],
                    'contrast': config['contrast'],
                    'input_sha256': signatures['input'],
                    'string_url': URL,
                },
            )
            shutil.copyfile(__file__, temporary / 'string_network.py')
            shutil.copyfile(config['module_file'], temporary / 'string_networks.py')
            if local_client:
                (temporary / 'utils').mkdir()
                for name in ('__init__.py', 'local_stringdb.py'):
                    shutil.copyfile(Path(config['module_file']).parent / 'utils' / name, temporary / 'utils' / name)
            shutil.copytree(REFERENCE_DIR, temporary / 'reference')
            manifest = {
                'schema': 1,
                'id': config['id'],
                'dataset': config['dataset'],
                'contrast': config['contrast'],
                'name': config['name'],
                'metadata': config['metadata'],
                'settings': spec,
                'selection': selection,
                'summary': summary,
                'created_at': created,
                'completed_at': time.time(),
                'creator': {},
                'string_url': URL,
                'database_source': database_source,
                'notebook_reference': 'ppis.ipynb',
                'notebook_sha256': NOTEBOOK_SHA256,
                'input_sha256': signatures['input'],
                'source_results_sha256': signatures['input'],
                'code_sha256': checksum(temporary / 'string_network.py'),
                'source_datasets': [config['dataset']],
                'sample_columns': columns,
                'source_analysis': {
                    'model': (
                        {'formula': config['model_formula']}
                        if config.get('model_formula')
                        else None
                    ),
                    'sample_annotations': annotations.to_dict('records'),
                    'pipeline': config.get('pipeline', 'jawm_string'),
                    'pipeline_revision': config['pipeline_revision'],
                    'string_module_revision': config.get(
                        'module_revision', 'unrecorded'
                    ),
                },
                'request_fingerprint': fingerprint,
                'runtime': identity,
                'cache_keys': client.keys,
                'cache_statistics': {
                    'entity_hits': getattr(client, 'hits', 0),
                    'entity_misses': getattr(client, 'misses', 0),
                },
                'disclosure': ('STRING bulk files were downloaded; gene mapping and network queries ran locally.' if local_client else 'Stable gene identifiers and species were sent to STRING by the pipeline; sample expression was retained locally.'),
            }
            # Publish every top-level artifact with the source comparison name.
            # Reference data and importable utility files retain their own names.
            for path in temporary.iterdir():
                if path.is_file():
                    path.rename(path.with_name(prefix + path.name))
            manifest['artifacts_sha256'] = {
                str(path.relative_to(temporary)): checksum(path)
                for path in temporary.rglob('*')
                if path.is_file()
            }
            atomic_json(temporary / (prefix + 'manifest.json'), manifest)
            if previous.exists():
                shutil.rmtree(previous)
            if output.exists():
                output.rename(previous)
            temporary.rename(output)
            if previous.exists():
                shutil.rmtree(previous)
            print(
                'Completed STRING network: ' + str(output / (prefix + 'network.graphml')),
                flush=True,
            )
            return manifest
        finally:
            if temporary.exists():
                shutil.rmtree(temporary)


if __name__ == '__main__':
    source = Path(r"""{{map.string_input}}""")
    project = Path(r"""{{project_folder}}""")
    contrast = source.name[: -len('.results.tsv')]
    output = Path(r"""{{mk.string_output}}""") / contrast
    REFERENCE_DIR = Path(r"""{{map.string_module}}""") / 'reference'
    metadata = {
        'organism': r"""{{organism}}""",
        'release': r"""{{release}}""",
        'taxon_id': r"""{{taxon_id}}""",
        'ncbi_taxon_id': r"""{{ncbi_taxon_id}}""",
        'assembly': r"""{{assembly}}""",
        'bioproject': r"""{{bioproject}}""",
        'geo': r"""{{geo}}""",
    }
    run(
        {
            'id': r"""{{network_id}}""",
            'dataset': project.name,
            'contrast': contrast,
            'name': contrast + ' · STRING (default)',
            'input': str(source),
            'samples': r"""{{map.string_samples}}""",
            'output': str(output),
            'cache': r"""{{mk.string_cache}}""",
            'STRINGDB': r"""{{STRINGDB}}""",
            'module_file': str(REFERENCE_DIR.parent / 'string_networks.py'),
            'metadata': {key: value for key, value in metadata.items() if value},
            'pipeline': r"""{{pipeline}}""",
            'pipeline_revision': r"""{{pipeline_revision}}""",
            'module_revision': r"""{{module_revision}}""",
            'model_formula': r"""{{model_formula}}""" or None,
            'container': r"""{{jawm.Process.container}}""",
        }
    )
''',
    # Defaults and variable descriptions follow the other JAWM modules.
    var={
        "map.string_module": str(ROOT),
        "input_file": "",
        "STRINGDB": "",
        "organism": "",
        "release": "",
        "taxon_id": "",
        "ncbi_taxon_id": "",
        "assembly": "",
        "bioproject": "",
        "geo": "",
        "model_formula": "",
        "pipeline": "jawm_string",
        "pipeline_revision": "unrecorded",
        "module_revision": "unrecorded",
    },
    desc={
        "project_folder": "Dataset directory containing deseq2_output/annotated.",
        "input_file": "Optional full path to one annotated .results.tsv contrast.",
        "string_input": "Annotated results table for this process clone (map.).",
        "string_samples": "Contrast sample table or samples_MasterTable.txt (map.).",
        "string_module": "Module directory containing the taxonomy reference (map.).",
        "string_output": "Parent network output directory (mk.), shared by contrast clones.",
        "network_id": "Stable network identifier recorded in the manifest.",
        "string_cache": "Shared persistent cache with cross-host POSIX flock (mk.).",
        "STRINGDB": "Optional bulk database root; STRINGDB/12.5/ holds shared files; empty uses the API.",
        "organism": "Scientific organism name, for example caenorhabditis_elegans.",
        "release": "Source Ensembl release, recorded in the network manifest.",
        "taxon_id": "Optional NCBI taxonomy ID; takes precedence over organism.",
        "ncbi_taxon_id": "Alternative NCBI taxonomy ID variable.",
        "assembly": "Optional source genome assembly.",
        "bioproject": "Optional source BioProject accession.",
        "geo": "Optional source GEO accession.",
        "model_formula": "Source differential-expression model, if known.",
        "pipeline": "Calling pipeline name for provenance.",
        "pipeline_revision": "Calling pipeline Git revision for provenance.",
        "module_revision": "STRING module Git revision for provenance.",
    },
    container=STRING_IMAGE,
    environment="docker",
    manager_slurm={"-c": 1, "--mem": "4GB", "-t": "24:00:00"},
)



string_database = jawm.Process(
    name="string_database",
    script=r'''#!/usr/bin/env python3
import sys
from pathlib import Path
folder = Path(r"""{{map.string_module}}""")
sys.path.insert(0, str(folder))
from utils.local_stringdb import ensure_database, resolve_species
species = resolve_species(
    r"""{{taxon_id}}""" or r"""{{ncbi_taxon_id}}""",
    r"""{{organism}}""", folder / 'reference/species.v12.5.tsv.gz')
ensure_database(r"""{{mk.string_db}}""", species)
''',
    var={"map.string_module": str(ROOT), "mk.string_db": "",
         "taxon_id": "", "ncbi_taxon_id": "", "organism": ""},
    desc={"string_db": "Shared STRING database root (mk.); versioned files and download.lock live below it.",
          "string_module": "STRING module directory (map.).",
          "taxon_id": "NCBI taxon ID.", "ncbi_taxon_id": "Alternative taxon ID.",
          "organism": "Scientific organism name."},
    container=STRING_IMAGE,
    environment="docker",
    manager_slurm={"-c": 1, "--mem": "4GB", "-t": "48:00:00"},
)


def database_root(variables):
    value = variables.get("STRINGDB") or os.environ.get("STRINGDB", "")
    return str(Path(value).expanduser().resolve()) if value else ""


def prepare_database(variables):
    """Prepare the optional shared database and wait before any network starts.

    The workers also validate readiness, since a scheduler may reuse a completed
    preparation job after the database has been moved or removed.
    """
    resolved = dict(string_network.var)
    resolved.update(variables)
    root = database_root(resolved)
    if not root:
        return
    process = string_database.clone()
    process.var.update({"mk.string_db": root, "map.string_module": str(ROOT)})
    for key in ("organism", "taxon_id", "ncbi_taxon_id"):
        process.var[key] = resolved.get(key) or ""
    process.execute()
    jawm.Process.wait(process.hash)


string_tests = jawm.Process(
    name="string_tests",
    script=r'''#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path

folder = Path(r"""{{map.string_module}}""")
suite = unittest.defaultTestLoader.discover(
    str(folder / "test"), pattern="test_*.py"
)
result = unittest.TextTestRunner(verbosity=2).run(suite)
if not result.wasSuccessful():
    sys.exit(1)

output = Path(r"""{{mk.test_output}}""")
output.mkdir(parents=True, exist_ok=True)
(output / "string_networks.txt").write_text("STRING module tests passed\n")
''',
    var={
        "map.string_module": str(ROOT),
        "mk.test_output": str(ROOT / "test/test-output"),
    },
    desc={
        "string_module": "Module directory containing the synthetic tests (map.).",
        "test_output": "Directory for the deterministic test result (mk.).",
    },
    container=STRING_IMAGE,
    environment="docker",
    manager_slurm={"-c": 1, "--mem": "4GB", "-t": "00:30:00"},
)


def git_revision(folder):
    """Return a source revision when the module is in a Git checkout."""
    try:
        return subprocess.check_output(
            ["git", "-C", str(folder), "rev-parse", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unrecorded"


def configure(result_file, variables, pipeline_file=None):
    """Clone the network process and assign its ordinary JAWM variables.

    Kept as the pipeline integration helper used by jawm_rnaseq.
    """
    result = Path(result_file).expanduser().resolve()
    if not result.name.endswith(".results.tsv"):
        raise ValueError("Expected an annotated .results.tsv file.")

    process = string_network.clone()
    # Read the resolved Process variables so process-scoped YAML and CLI
    # overrides work in addition to the caller's dataset variables.
    var = dict(process.var)
    var.update(variables)
    process.var["STRINGDB"] = database_root(var)
    if process.var["STRINGDB"]:
        process.var["mk.string_db"] = process.var["STRINGDB"]
    project = Path(var["project_folder"]).expanduser().resolve()
    contrast = result.name[: -len(".results.tsv")]
    identifier = uuid.uuid5(
        uuid.NAMESPACE_URL, "jawm-rnaseq:string:12.5:" + contrast
    ).hex
    samples = result.parent.parent / (contrast + ".input.tsv")
    if not samples.is_file():
        samples = result.parent.parent / "samples_MasterTable.txt"
    if not samples.is_file():
        raise ValueError("Missing contrast sample metadata: " + str(samples))

    for key in [
        "organism",
        "release",
        "taxon_id",
        "ncbi_taxon_id",
        "assembly",
        "bioproject",
        "geo",
    ]:
        process.var[key] = var.get(key) if var.get(key) is not None else ""
    cache = (
        Path(var.get("string_cache") or project.parent / ".jawm-string-cache")
        .expanduser()
        .resolve()
    )
    process.var.update(
        {
            "project_folder": str(project),
            "input_file": str(result),
            "map.string_input": str(result),
            "map.string_samples": str(samples),
            "map.string_module": str(ROOT),
            "mk.string_output": str(project / "string-networks"),
            "network_id": identifier,
            "mk.string_cache": str(cache),
            "pipeline": "jawm_rnaseq" if pipeline_file else "jawm_string",
            "pipeline_revision": git_revision(
                Path(pipeline_file).resolve().parent if pipeline_file else ROOT
            ),
            "module_revision": git_revision(ROOT),
            "model_formula": var.get("model_formula")
            or ("~ group" if pipeline_file else ""),
        }
    )
    return process


def report_files(project_folder):
    """Return network artifacts for the calling pipeline's report manifest."""
    folder = Path(project_folder) / "string-networks"
    files = sorted(folder.glob("*/*_network.graphml"))
    files += sorted(folder.glob("*/*_manifest.json"))
    return {"STRING networks": [str(path) for path in files]}


def validate_test_outputs(folders):
    """Check real network artifacts and write stable tables for jawm-test hashes.

    Layout coordinates, timestamps, paths and cache counters are intentionally
    excluded; node biology and STRING interaction scores are compared.
    """
    import csv
    import json
    import xml.etree.ElementTree as ET

    def read_table(path):
        with path.open(newline="") as stream:
            return list(csv.DictReader(stream, delimiter="\t"))

    def write_table(path, columns, rows):
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(
                stream, fieldnames=columns, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)

    for folder in folders:
        prefix = folder.name + "_"
        for name in [
            "network.graphml",
            "network.json",
            "nodes.tsv",
            "edges.tsv",
            "manifest.json",
            "inputs.tsv",
            "samples.tsv",
            "mapping.tsv",
        ]:
            if not (folder / (prefix + name)).is_file():
                raise ValueError("Missing real test artifact: " + str(folder / (prefix + name)))
        manifest = json.loads((folder / (prefix + "manifest.json")).read_text())
        nodes = read_table(folder / (prefix + "nodes.tsv"))
        edges = read_table(folder / (prefix + "edges.tsv"))
        if len(nodes) < 2 or not edges:
            raise ValueError(
                "The live STRING integration test produced an empty network."
            )
        if manifest["summary"]["nodes"] != len(nodes) or manifest["summary"][
            "edges"
        ] != len(edges):
            raise ValueError("Manifest and network table sizes differ.")
        sample_columns = manifest["sample_columns"]
        if not sample_columns or not any(row.get(sample_columns[0]) for row in nodes):
            raise ValueError(
                "Measured sample expression is missing from network nodes."
            )
        snapshot = json.loads((folder / (prefix + "network.json")).read_text())
        graphml = ET.parse(folder / (prefix + "network.graphml"))
        namespace = {"g": "http://graphml.graphdrawing.org/xmlns"}
        if (
            len(snapshot["nodes"]) != len(nodes)
            or len(snapshot["edges"]) != len(edges)
            or len(graphml.findall(".//g:node", namespace)) != len(nodes)
            or len(graphml.findall(".//g:edge", namespace)) != len(edges)
        ):
            raise ValueError("GraphML, JSON and network tables disagree.")

        # Preserve measured node statistics and expression, but omit x/y layout.
        node_columns = [
            "stringId",
            "name",
            "role",
            "expression_gene_id",
            "log2FoldChange",
            "padj",
            *sorted(sample_columns),
        ]
        canonical_nodes = []
        for row in sorted(nodes, key=lambda item: item["stringId"]):
            record = {key: row.get(key, "") for key in node_columns}
            for key in ["log2FoldChange", "padj", *sample_columns]:
                if record[key]:
                    record[key] = format(float(record[key]), ".12g")
            canonical_nodes.append(record)
        write_table(folder / (prefix + "test_nodes.tsv"), node_columns, canonical_nodes)

        # An undirected edge has the same identity in either API orientation.
        score_columns = [
            "score",
            "nscore",
            "fscore",
            "pscore",
            "ascore",
            "escore",
            "dscore",
            "tscore",
        ]
        edge_columns = ["stringId_A", "stringId_B", *score_columns, "edge_role"]
        identifiers = {node["stringId"] for node in nodes}
        canonical_edges = []
        pairs = set()
        for row in edges:
            a, b = sorted([row["stringId_A"], row["stringId_B"]])
            if (
                a == b
                or a not in identifiers
                or b not in identifiers
                or (a, b) in pairs
            ):
                raise ValueError("Invalid or duplicate network edge.")
            if float(row["score"]) < manifest["settings"]["min_combined_score"]:
                raise ValueError("A network edge is below the requested confidence.")
            pairs.add((a, b))
            canonical_edges.append(
                {
                    "stringId_A": a,
                    "stringId_B": b,
                    "edge_role": row["edge_role"],
                    **{key: format(float(row[key]), ".6f") for key in score_columns},
                }
            )
        canonical_edges.sort(key=lambda row: (row["stringId_A"], row["stringId_B"]))
        write_table(folder / (prefix + "test_edges.tsv"), edge_columns, canonical_edges)
        logger.info(
            "Real STRING test: %s nodes, %s edges; artifacts in %s",
            len(nodes),
            len(edges),
            folder,
        )


if __name__ == "__main__":
    import sys
    from jawm.utils import workflow

    workflows, var, args, unknown_args = jawm.utils.parse_arguments(
        ["main", "string", "test", "unit_test"]
    )

    if workflow(["main", "string", "test"], workflows):
        # Use the Process's resolved variables, including process-scoped YAML.
        var = string_network.var
        project = Path(var["project_folder"]).expanduser().resolve()
        if var.get("input_file"):
            results = [Path(var["input_file"])]
        else:
            results = sorted(
                (project / "deseq2_output/annotated").glob("*.results.tsv")
            )
        if not results:
            raise ValueError("No annotated contrast results found in " + str(project))

        prepare_database(var)
        string_jobs = []
        string_outputs = []
        for result in results:
            string_network_ = configure(result, var)
            string_network_.execute()
            string_jobs.append(string_network_.hash)
            string_outputs.append(
                Path(string_network_.var["mk.string_output"])
                / result.name[: -len(".results.tsv")]
            )

        jawm.Process.wait(string_jobs)

    if workflow("test", workflows):
        validate_test_outputs(string_outputs)

    if workflow("unit_test", workflows):
        string_tests.execute()
        jawm.Process.wait(string_tests.hash)
        logger.info("Test completed.")

    sys.exit(0)
