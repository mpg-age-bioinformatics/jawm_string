import hashlib
import importlib.util
import json
import multiprocessing
from pathlib import Path
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
CHANNELS=['nscore','fscore','pscore','ascore','escore','dscore','tscore']
WORKER=None


def load_worker():
    # Extract the embedded script without importing JAWM on the job node.
    import ast
    tree=ast.parse((ROOT/'string_networks.py').read_text())
    process=next(node.value for node in tree.body if isinstance(node,ast.Assign) and any(isinstance(target,ast.Name) and target.id=='string_network' for target in node.targets))
    script=next(ast.literal_eval(keyword.value) for keyword in process.keywords if keyword.arg=='script')
    temporary=tempfile.TemporaryDirectory()
    path=Path(temporary.name)/'string_worker.py';path.write_text(script)
    spec=importlib.util.spec_from_file_location('string_worker',path)
    worker=importlib.util.module_from_spec(spec);spec.loader.exec_module(worker)
    worker.REFERENCE_DIR=ROOT/'reference'
    return worker,temporary


def edge(species,a,b):
    return dict(stringId_A=str(species)+'.'+a,stringId_B=str(species)+'.'+b,
        preferredName_A=a,preferredName_B=b,score=.9,**dict.fromkeys(CHANNELS,.8))


def response(client,method,ids,**parameters):
    if method=='get_string_ids':return [dict(queryIndex=i,ncbiTaxonId=client.species,stringId=str(client.species)+'.'+gene,preferredName=gene) for i,gene in enumerate(ids) if gene!='unmapped']
    edges=[edge(client.species,a,b) for a,b in [('a','b'),('a','x'),('b','x'),('c','d'),('d','x'),('d','y')]]
    return [item for item in edges if item['stringId_A'] in ids or item['stringId_B'] in ids]


def query_in_child(cache):
    WORKER.CachedStringClient(cache,6239).request('get_string_ids',['a','b'],echo_query=1)


class StringTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name).resolve()
        global WORKER
        self.worker,temporary=load_worker();self.addCleanup(temporary.cleanup);WORKER=self.worker
        self.cache=self.root/'cache';self.calls=[]
        def fake(client,method,ids,**params):
            self.calls.append((client.species,method,list(ids),params));return response(client,method,ids,**params)
        patcher=patch.object(self.worker.StringClient,'request',fake);patcher.start();self.addCleanup(patcher.stop)

    def test_overlapping_mapping_batches_share_positive_and_negative_results(self):
        client=self.worker.CachedStringClient(self.cache,6239)
        client.request('get_string_ids',['a','b','unmapped'],echo_query=1)
        result=self.worker.CachedStringClient(self.cache,6239).request('get_string_ids',['unmapped','b','c'],echo_query=1)
        self.assertEqual(self.calls[-1][2],['c'])
        self.assertEqual([(r['queryIndex'],r['stringId']) for r in result],[(1,'6239.b'),(2,'6239.c')])
        self.assertEqual(len(self.calls),2)
        self.worker.CachedStringClient(self.cache,6239).request('get_string_ids',['b','a'],echo_query=1)
        self.assertEqual(len(self.calls),2)
        self.worker.CachedStringClient(self.cache,9606).request('get_string_ids',['b','a'],echo_query=1)
        self.assertEqual(len(self.calls),3)

    def test_partner_cache_reuses_overlaps_and_separates_network_type_and_score(self):
        client=self.worker.CachedStringClient(self.cache,6239);params={'required_score':700,'network_type':'functional'}
        client.request('interaction_partners',['6239.a','6239.b'],**params)
        result=client.request('interaction_partners',['6239.b','6239.x'],**params)
        self.assertEqual(self.calls[-1][2],['6239.x'])
        self.assertEqual(len({json.dumps(row,sort_keys=True) for row in result}),len(result))
        client.request('interaction_partners',['6239.a','6239.x'],**params);self.assertEqual(len(self.calls),2)
        client.request('interaction_partners',['6239.a'],**dict(params,network_type='physical'))
        client.request('interaction_partners',['6239.a'],**dict(params,required_score=900));self.assertEqual(len(self.calls),4)

    def test_failures_are_not_cached_as_absent_genes(self):
        with patch.object(self.worker.StringClient,'request',side_effect=ValueError('unavailable')):
            with self.assertRaises(ValueError):self.worker.CachedStringClient(self.cache,6239).request('get_string_ids',['a'],echo_query=1)
        self.assertEqual(self.worker.CachedStringClient(self.cache,6239).request('get_string_ids',['a'],echo_query=1)[0]['stringId'],'6239.a')

    def test_parallel_jobs_share_a_single_cache_miss(self):
        log=self.root/'http.log'
        def fake(client,method,ids,**params):
            with log.open('a') as stream:stream.write('call\n')
            time.sleep(.15)
            return response(client,method,ids,**params)
        with patch.object(self.worker.StringClient,'request',fake):
            context=multiprocessing.get_context('fork')
            children=[context.Process(target=query_in_child,args=(self.cache,)) for _ in range(3)]
            for child in children:child.start()
            for child in children:child.join(10);self.assertEqual(child.exitcode,0)
        self.assertEqual(log.read_text().splitlines(),['call'])

    def config(self,contrast):
        source=self.root/(contrast+'.results.tsv');samples=self.root/'samples.tsv'
        self.worker.pd.DataFrame({'ensembl_gene_id':['a','b','c','d','x','unmapped'],
            'log2FoldChange':[-4,3,-2,1,1,0],'padj':[.001,.001,.001,.001,.8,.9],
            'sample.rep_1':[10,20,30,40,50,0]}).to_csv(source,sep='\t',index=False)
        samples.write_text('sample\tgroup\nsample.rep_1\tA\n')
        return dict(id=contrast,dataset='study',contrast=contrast,name=contrast+' STRING',input=str(source),samples=str(samples),
            output=str(self.root/'networks'/contrast),cache=str(self.cache),metadata={'organism':'caenorhabditis_elegans','release':'115'},
            pipeline_revision='test-commit',module_file=str(ROOT/'string_networks.py'),container='mpgagebioinformatics/networkx:3.5')

    def test_offline_network_never_calls_api_and_records_provenance(self):
        from test_local_stringdb import local, fake_download
        config = self.config('offline')
        config['STRINGDB'] = str(self.root / 'database')
        frame = self.worker.pd.read_csv(config['input'], sep='\t')
        frame['ensembl_gene_id'] = ['WBGene1', 'WBGene2', 'absent1', 'absent2', 'absent3', 'absent4']
        frame.to_csv(config['input'], sep='\t', index=False)
        with patch.object(local, '_download', side_effect=fake_download), patch.object(
                self.worker.CachedStringClient, 'request', side_effect=AssertionError('API called')):
            manifest = self.worker.run(config)
        self.assertEqual(manifest['database_source']['backend'], 'local')
        self.assertEqual(manifest['summary']['edges'], 1)
        self.assertFalse(self.calls)
        with patch.object(local, '_download', side_effect=AssertionError('download repeated')):
            self.assertEqual(self.worker.run(config), manifest)
        self.assertTrue((Path(config['output']) / 'utils' / 'local_stringdb.py').is_file())

    def test_online_404_nothing_found_list_and_dict_only(self):
        import io
        from urllib.error import HTTPError
        worker, temporary = load_worker()
        self.addCleanup(temporary.cleanup)
        for body in [[{'Error': 'nothing found'}], {'Error': 'nothing found'},
                     [{'Error': 'different failure'}], {'unexpected': 'response'}]:
            client = worker.StringClient(self.cache, 6239)
            error = HTTPError('https://example.invalid', 404, 'Not Found', {},
                              io.BytesIO(json.dumps(body).encode()))
            with patch.object(client.opener, 'open', side_effect=error), patch.object(worker.time, 'sleep'):
                if body in [[{'Error': 'nothing found'}], {'Error': 'nothing found'}]:
                    self.assertEqual(client.request('get_string_ids', ['missing']), [])
                else:
                    with self.assertRaisesRegex(ValueError, 'HTTP 404'):
                        client.request('get_string_ids', ['missing'])

    def test_one_network_per_contrast_reuses_queries_and_is_restartable(self):
        first=self.config('one');second=self.config('two')
        manifest=self.worker.run(first);self.assertEqual(len(self.calls),3)
        self.assertEqual(manifest['settings'],dict(self.worker.DEFAULTS,species=6239))
        other=self.worker.run(second);self.assertEqual(len(self.calls),3);self.assertGreater(other['cache_statistics']['entity_hits'],0)
        graph=json.loads((Path(first['output'])/'one_network.json').read_text())
        self.assertEqual(graph['summary']['seed_proteins'],4);self.assertEqual(graph['summary']['added_proteins'],1)
        self.assertEqual(next(n for n in graph['nodes'] if n['stringId']=='6239.x')['sample.rep_1'],50)
        self.assertEqual(manifest['sample_columns'],['sample.rep_1'])
        with patch.object(self.worker,'build_network',side_effect=AssertionError('must resume')):self.assertEqual(self.worker.run(first),manifest)
        source=Path(first['input']);source.write_text(source.read_text().replace('50','75'))
        updated=self.worker.run(first);self.assertNotEqual(updated['request_fingerprint'],manifest['request_fingerprint']);self.assertEqual(len(self.calls),3)
        before=(Path(first['output'])/'one_manifest.json').read_bytes();source.write_text(source.read_text().replace('75','80'))
        with patch.object(self.worker,'build_network',side_effect=ValueError('failure')):
            with self.assertRaises(ValueError):self.worker.run(first)
        self.assertEqual((Path(first['output'])/'one_manifest.json').read_bytes(),before)

    def test_published_files_use_comparison_prefix_and_checksums(self):
        config = self.config('daf12_vs_WT')
        manifest = self.worker.run(config)
        folder = Path(config['output'])
        files = [path for path in folder.iterdir() if path.is_file()]
        self.assertTrue(files)
        self.assertTrue(all(path.name.startswith('daf12_vs_WT_') for path in files))
        self.assertTrue((folder / 'reference/species.v12.5.tsv.gz').is_file())
        for name, digest in manifest['artifacts_sha256'].items():
            self.assertEqual(self.worker.checksum(folder / name), digest)
        # A damaged published artifact must trigger a rebuild using cached queries.
        calls = len(self.calls)
        (folder / 'daf12_vs_WT_nodes.tsv').write_text('damaged')
        with patch.object(self.worker, 'build_network', wraps=self.worker.build_network) as build:
            self.worker.run(config)
            build.assert_called_once()
        self.assertEqual(len(self.calls), calls)

    def test_default_network_graphml_and_empty_selection(self):
        config=self.config('graphml');manifest=self.worker.run(config)
        graph=self.worker.nx.read_graphml(Path(config['output'])/'graphml_network.graphml')
        self.assertEqual(graph.number_of_nodes(),5);self.assertEqual(graph.number_of_edges(),5)
        self.assertEqual(graph.nodes['6239.x']['sample.rep_1'],50)
        self.assertEqual(graph.nodes['6239.x']['role'],'added');self.assertFalse(graph.nodes['6239.x']['is_DE'])
        self.assertEqual(manifest['code_sha256'],self.worker.checksum(Path(config['output'])/'graphml_string_network.py'))
        config=self.config('empty');source=Path(config['input'])
        frame=self.worker.pd.read_csv(source,sep='\t');frame['padj']=1;frame.to_csv(source,sep='\t',index=False)
        count=len(self.calls);result=self.worker.run(config)
        self.assertEqual(result['summary']['nodes'],0);self.assertEqual(len(self.calls),count)

    def test_taxonomy_validation_and_negative_fold_change_ranking(self):
        with self.assertRaises(ValueError):self.worker.validate({}, {'organism':'unrecognized species'})
        with self.assertRaises(ValueError):self.worker.validate({'species':9606},{'organism':'caenorhabditis_elegans'})
        config=self.config('ranking');frame=self.worker.pd.read_csv(config['input'],sep='\t')
        frame['log2FoldChange']=frame.log2FoldChange.astype(float)
        frame.loc[0,'padj']=0;frame.loc[2,'padj']=-1;frame.loc[3,'log2FoldChange']=float('inf')
        _,selected,summary=self.worker.rank_genes(frame,self.worker.validate({},config['metadata']))
        self.assertEqual(selected.iloc[0].ensembl_gene_id,'a');self.assertEqual(selected.iloc[0].signed_score,-1200)
        self.assertEqual(summary['invalid_rows'],2)

    def test_real_output_validation_hashes_biology_but_ignores_layout(self):
        import ast
        import csv
        import logging

        # Exercise the integration validator without importing the controller.
        tree = ast.parse((ROOT / 'string_networks.py').read_text())
        definition = next(node for node in tree.body
                          if isinstance(node, ast.FunctionDef)
                          and node.name == 'validate_test_outputs')
        namespace = {'logger': logging.getLogger('test_string_validator')}
        exec(compile(ast.Module(body=[definition], type_ignores=[]),
                     'validate_test_outputs', 'exec'), namespace)
        validate = namespace['validate_test_outputs']
        config = self.config('validation')
        self.worker.run(config)
        folder = Path(config['output'])
        validate([folder])
        before = [(folder / ('validation_' + name)).read_bytes()
                  for name in ['test_nodes.tsv', 'test_edges.tsv']]

        def edit_table(name, change):
            path = folder / ('validation_' + name)
            with path.open(newline='') as stream:
                reader = csv.DictReader(stream, delimiter='\t')
                columns = reader.fieldnames
                rows = list(reader)
            change(rows)
            with path.open('w', newline='') as stream:
                writer = csv.DictWriter(stream, fieldnames=columns,
                                        delimiter='\t', lineterminator='\n')
                writer.writeheader()
                writer.writerows(rows)

        edit_table('nodes.tsv', lambda rows: rows[0].update(x='12345', y='67890'))
        validate([folder])
        self.assertEqual(before, [(folder / ('validation_' + name)).read_bytes()
                                 for name in ['test_nodes.tsv', 'test_edges.tsv']])
        edit_table('edges.tsv', lambda rows: rows[0].update(score='0.95'))
        validate([folder])
        self.assertNotEqual(before[1], (folder / 'validation_test_edges.tsv').read_bytes())
        edit_table('edges.tsv', lambda rows: rows.clear())
        with self.assertRaisesRegex(ValueError, 'empty network'):
            validate([folder])

    @unittest.skipUnless(importlib.util.find_spec('jawm'),'JAWM controller integration is tested separately from the scientific runtime.')
    def test_repository_loader_needs_no_scientific_imports(self):
        from jawm.utils import load_modules
        with patch.dict(sys.modules,{'numpy':None,'pandas':None,'networkx':None}):
            import shutil
            modules=self.root/'modules';shutil.copytree(ROOT,modules/'jawm_string',ignore=shutil.ignore_patterns('.git','logs','__pycache__','.submodules','test-input','test-output'))
            loaded=load_modules(['jawm_string'],modules_root=modules)
        self.assertTrue('string_networks' in loaded or 'string_networks' in sys.modules)

    @unittest.skipUnless(importlib.util.find_spec('jawm'), 'JAWM controller integration requires JAWM.')
    def test_database_preparation_opt_in_and_rendered_script(self):
        import os
        import subprocess
        from test_local_stringdb import local, fake_download
        spec = importlib.util.spec_from_file_location('string_networks', ROOT / 'string_networks.py')
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with patch.dict(os.environ, {'STRINGDB': ''}), patch.object(module.string_database, 'clone') as clone:
            module.prepare_database({})
            clone.assert_not_called()
        root = self.root / 'shared database'
        with patch.object(local, '_download', side_effect=fake_download):
            local.ensure_database(root, 6239)
        process = module.string_database.clone()
        with patch.dict(os.environ, {'STRINGDB': str(root)}), patch.object(
                module.string_database, 'clone', return_value=process), patch.object(
                process, 'execute') as execute, patch.object(module.jawm.Process, 'wait') as wait:
            module.prepare_database({'organism': 'caenorhabditis_elegans'})
            execute.assert_called_once()
            wait.assert_called_once_with(process.hash)
        rendered = process._script_placeholders_and_mkdir(process.script)
        job = self.root / 'prepare.py'
        job.write_text(rendered)
        subprocess.run([sys.executable, str(job)], check=True, capture_output=True, text=True)
        self.assertFalse((root / '12.5/download.lock').exists())
        variables = {'STRINGDB': str(root), 'project_folder': str(self.root), 'organism': 'caenorhabditis_elegans'}
        result = self.root / 'deseq2_output/annotated/contrast.results.tsv'
        result.parent.mkdir(parents=True)
        (result.parent.parent / 'samples_MasterTable.txt').write_text('sample\tgroup\ns\tA\n')
        network = module.configure(result, variables)
        self.assertEqual(network.var['STRINGDB'], str(root))
        self.assertEqual(network.var['mk.string_db'], str(root))

    @unittest.skipUnless(importlib.util.find_spec('jawm'),'Controller integration requires JAWM; scientific/cache tests run without it.')
    def test_process_clones_and_python_string_name_does_not_collide(self):
        import string as standard_string
        spec=importlib.util.spec_from_file_location('string_networks',ROOT/'string_networks.py')
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        compile(module.string_tests.script,'string_tests','exec')
        project=self.root/"study with spaces and 'quotes'";folder=project/'deseq2_output/annotated';folder.mkdir(parents=True)
        (folder.parent/'samples_MasterTable.txt').write_text('sample\tgroup\ns\tA\n')
        processes=[module.configure(folder/(name+'.results.tsv'),{'project_folder':str(project),'organism':'homo_sapiens'}) for name in ['one','two']]
        self.assertNotEqual(processes[0].var['network_id'],processes[1].var['network_id'])
        self.assertEqual(processes[0].var['mk.string_output'], str(project/'string-networks'))
        self.assertEqual(processes[0].var['mk.string_output'],processes[1].var['mk.string_output'])
        self.assertEqual(processes[0].var['mk.string_cache'],processes[1].var['mk.string_cache'])
        self.assertEqual(processes[0].container,'mpgagebioinformatics/networkx:3.5')
        self.assertEqual(processes[0].var['organism'],'homo_sapiens')
        self.assertIn('map.string_module',processes[0].var);self.assertIn('ascii_letters',dir(standard_string))
        self.assertEqual(self.worker.species_hint({'organism':'Xenopus_tropicalis'}),8364)
        # Execute the actual rendered job script with no DE genes (and no HTTP).
        source=folder/'one.results.tsv'
        source.write_text('ensembl_gene_id\tlog2FoldChange\tpadj\ts\nA\t1\t1\t10\n')
        rendered=processes[0]._script_placeholders_and_mkdir(processes[0].script)
        job=self.root/'rendered_job.py';job.write_text(rendered)
        import subprocess
        subprocess.run([sys.executable,str(job)],check=True,capture_output=True,text=True)
        completed=json.loads((Path(processes[0].var['mk.string_output'])/'one'/'one_manifest.json').read_text())
        self.assertEqual(completed['summary']['status'],'no_DE_genes')
        self.assertEqual(completed['sample_columns'],['s'])
        outputs = project / 'string-networks/one'
        self.assertEqual(module.report_files(project), {'STRING networks': [
            str(outputs / 'one_network.graphml'), str(outputs / 'one_manifest.json')]})
        self.assertFalse((project / 'metadata/string-networks').exists())

if __name__=='__main__':unittest.main()
