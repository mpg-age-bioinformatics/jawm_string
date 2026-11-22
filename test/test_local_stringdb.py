import gzip
import importlib.util
import json
import multiprocessing
from pathlib import Path
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from utils import local_stringdb as local

INFO = '#string_protein_id\tpreferred_name\tprotein_size\tannotation\n6239.a\talpha\t10\tA\n6239.b\tbeta\t20\tB\n6239.c\tgamma\t30\tC\n'
ALIASES = '#string_protein_id\talias\tsource\n6239.a\tWBGene1\tEnsembl\n6239.b\tWBGene2\tEnsembl\n6239.b\tambiguous\tEnsembl\n6239.a\tambiguous\tEnsembl\n'
LINKS = 'protein1 protein2 neighborhood fusion cooccurence coexpression experimental database textmining combined_score\n6239.a 6239.b 100 200 300 400 500 600 700 900\n6239.b 6239.a 100 200 300 400 500 600 700 900\n6239.b 6239.c 0 0 0 100 200 300 400 500\n'


def fake_download(url, target):
    text = INFO if 'protein.info' in url else ALIASES if 'protein.aliases' in url else LINKS
    with gzip.open(target, 'wt') as stream:
        stream.write(text)


def concurrent_prepare(root, log):
    def download(url, target):
        with open(log, 'a') as stream:
            stream.write(url + '\n')
        time.sleep(.03)
        fake_download(url, target)
    with patch.object(local, '_download', download):
        local.ensure_database(root, 6239, poll=.01)


class LocalTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_mapping_edges_scores_and_reuse_without_network(self):
        with patch.object(local, '_download', side_effect=fake_download) as download:
            client = local.LocalStringClient(self.root, 6239)
            self.assertEqual(download.call_count, 3)
        with patch.object(local, '_download', side_effect=AssertionError('unexpected download')):
            client = local.LocalStringClient(self.root, 6239)
        mapped = client.request('get_string_ids', ['absent', 'WBGene1', 'ambiguous', 'beta'])
        self.assertEqual([(r['queryIndex'], r['stringId']) for r in mapped], [(1, '6239.a'), (2, '6239.a'), (3, '6239.b')])
        edges = client.request('interaction_partners', ['6239.a', '6239.b'], required_score=700)
        self.assertEqual(len(edges), 1)
        self.assertEqual(edges[0]['score'], .9)
        self.assertEqual(edges[0]['nscore'], .1)
        self.assertEqual(edges[0]['tscore'], .7)
        self.assertEqual(edges[0]['preferredName_A'], 'alpha')
        self.assertEqual(len(client.request('interaction_partners', ['6239.c'], required_score=400)), 1)
        self.assertFalse((self.root / '12.5/download.lock').exists())
        self.assertEqual(client.manifest['mapping_policy'], local.MAPPING_POLICY)

    def test_concurrent_processes_download_once(self):
        context = multiprocessing.get_context('fork')
        log = self.root / 'requests'
        processes = [context.Process(target=concurrent_prepare, args=(self.root, log)) for _ in range(3)]
        for process in processes:
            process.start()
        for process in processes:
            process.join(10)
            if process.is_alive():
                process.terminate()
                process.join()
                self.fail('Downloader did not finish')
            self.assertEqual(process.exitcode, 0)
        self.assertEqual(len(log.read_text().splitlines()), 3)

    def test_existing_lock_blocks_even_ready_database(self):
        with patch.object(local, '_download', side_effect=fake_download):
            folder, _ = local.ensure_database(self.root, 6239)
        lock = folder / 'download.lock'
        lock.write_text('external owner')
        with self.assertRaises(TimeoutError):
            local.ensure_database(self.root, 6239, timeout=.02, poll=.01)
        self.assertEqual(lock.read_text(), 'external owner')

    def test_download_failure_does_not_publish_ready_and_can_retry(self):
        with patch.object(local, '_download', side_effect=OSError('interrupted')):
            with self.assertRaises(OSError):
                local.ensure_database(self.root, 6239)
        self.assertFalse((self.root / '12.5/6239.functional.manifest.json').exists())
        self.assertFalse((self.root / '12.5/download.lock').exists())
        with patch.object(local, '_download', side_effect=fake_download):
            folder, manifest = local.ensure_database(self.root, 6239)
        (folder / manifest['database']).write_bytes(b'broken')
        with patch.object(local, '_download', side_effect=fake_download) as download:
            local.ensure_database(self.root, 6239)
            self.assertEqual(download.call_count, 3)

    def test_malformed_download_is_not_published(self):
        def invalid(url, target):
            fake_download(url, target)
            if 'links' in url:
                with gzip.open(target, 'wt') as stream:
                    stream.write('bad header\n')
        with patch.object(local, '_download', side_effect=invalid):
            with self.assertRaises(ValueError):
                local.ensure_database(self.root, 6239)
        self.assertFalse((self.root / '12.5/6239.functional.manifest.json').exists())

    def test_physical_links_have_only_three_evidence_channels(self):
        def physical(url, target):
            fake_download(url, target)
            if 'links' in url:
                with gzip.open(target, 'wt') as stream:
                    stream.write('protein1 protein2 experimental database textmining combined_score\n6239.a 6239.b 800 0 0 800\n')
        with patch.object(local, '_download', side_effect=physical):
            client = local.LocalStringClient(self.root, 6239, network_type='physical')
        rows = client.request('interaction_partners', ['6239.a'], network_type='physical', required_score=700)
        self.assertEqual(rows[0]['escore'], .8)
        self.assertEqual(rows[0]['nscore'], 0)
        with self.assertRaises(ValueError):
            client.request('interaction_partners', ['6239.a'], network_type='functional')

    def test_species_resolution(self):
        reference = ROOT / 'reference/species.v12.5.tsv.gz'
        self.assertEqual(local.resolve_species('', 'caenorhabditis_elegans', reference), 6239)
        self.assertEqual(local.resolve_species('9606', '', reference), 9606)
        with self.assertRaises(ValueError):
            local.resolve_species('', 'not_a_species', reference)


if __name__ == '__main__':
    unittest.main()
