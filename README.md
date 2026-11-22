# jawm_string

This is a jawm STRING module.

Installing jawm:
```
pip install git+ssh://git@github.com/mpg-age-bioinformatics/jawm.git
```
For more information on jawm please visit jawm's repo on [GitHub.com](https://github.com/mpg-age-bioinformatics/jawm/tree/main).

Example usage:
```
# clone this module
git clone git@github.com:mpg-age-bioinformatics/jawm_string.git

cd jawm_string

# download test data
jawm-test -r download

# docker
jawm string_networks.py string -p ./yaml/docker.yaml

# slurm & apptainer with multiple yaml files
jawm string_networks.py string -p ./yaml/docker.yaml ./yaml/hpc.yaml ./yaml/nexus.apptainer.yaml
```

Additional jawm workflows are available [here (GitHub.com)](https://github.com/mpg-age-bioinformatics?q=jawm_&type=all&language=&sort=).

## STRING-specific details

A JAWM module that creates one default STRING protein-interaction network per
annotated differential-expression contrast. It uses the RNA-seq app's selection,
expansion, filtering and visualization-data defaults, with a persistent shared
cache designed for datasets containing many comparisons.

`string_networks.py` contains the JAWM processes and their embedded worker scripts.
Importing the module on the controller requires JAWM and the Python standard
library; scientific dependencies run inside the job container.

As in the reference modules, `string_network` and `string_tests` are explicit
`jawm.Process` objects with inline `script`, `var`, `desc`, container and Slurm
settings. Inputs use `map.*` variables and writable directories use `mk.*`.
`configure()` clones `string_network` and assigns those named variables for one
contrast. The command-line entry point uses `parse_arguments`, workflow selection,
cloning, `execute()` and `Process.wait()`. Dataset variables can be supplied at
global or process scope in YAML; process variables are read from the resolved
Process configuration.

`mk.string_output` maps the parent `metadata/string-networks` directory, and
`network_id` selects a contrast subdirectory. Mapping the parent allows the worker
to atomically replace a completed contrast directory inside the container.

### Module organization

The layout follows the JAWM module conventions, including the shared test
workflow, parameter scopes and inline output-hash configuration.

```text
.github/workflows/test.yaml       # Shared JAWM module test workflow
string_networks.py                # Importable processes and workflow entry point
utils/
    __init__.py
    local_stringdb.py             # Shared downloads and local SQLite adapter
containers/string/
    Dockerfile
    requirements.txt
reference/
    species.v12.5.tsv.gz
    species.v12.5.source.json
    workflow-source.json
test/
    README.md
    data.txt                     # MD5/download entry; pending Figshare URL
    tests.txt                    # Semicolon-delimited jawm-test entry
    test_string_module.py        # Synthetic scientific/cache and controller tests
    yaml/
        test.yaml                # Real airway input and output hash scope
        unit_test.yaml           # Optional synthetic unit-test process
yaml/
    build.yaml                   # Optional /nexus build-machine overrides
    docker.yaml                  # Docker settings and dataset variables
    hpc.yaml                     # Slurm/Apptainer settings and /nexus bind
    nexus.apptainer.yaml          # Per-process site SIF paths
```

STRING requires the container recipe, taxonomy catalog, provenance and synthetic
tests in addition to the common module files. The Python filename remains
`string_networks.py` to preserve pipeline imports and avoid shadowing Python's
standard-library `string` module.

### Runtime image

Build the runtime image locally if needed:

```sh
docker build -t mpgagebioinformatics/networkx:3.5 -f containers/string/Dockerfile .
```

The runtime pins Python 3.13, NetworkX 3.5, NumPy 2.3.5 and pandas 2.3.3, matching
the app. The image contains dependencies; JAWM's job script contains the worker
code. JAWM mounts this module's directory for its taxonomy catalog and provenance,
as well as the result table, sample metadata, output directory and shared cache.
Publish the image before remote Docker workers use it. For Apptainer:

```sh
apptainer pull networkx-3.5.sif docker://mpgagebioinformatics/networkx:3.5
```

Override `string_network` and `string_tests` containers with that SIF, as illustrated
in `yaml/nexus.apptainer.yaml`. Set your site's partition/account in
`yaml/hpc.yaml`.
Use `automated_mount: false` only if explicit bind/PVC configuration already makes
all input, module, output and cache paths visible inside the jobs.

### Use from another JAWM pipeline

```python
from jawm.utils import load_modules

load_modules(["jawm_string"])  # exports the string_networks module

process = string_networks.configure(
    "/lake/STUDY/deseq2_output/annotated/group_B_vs_A.results.tsv",
    {"project_folder": "/lake/STUDY", "organism": "homo_sapiens", "release": "115"},
    pipeline_file=__file__,
)
process.execute(annotator.hash)
```

For all contrasts, loop over annotated `*.results.tsv` files and configure one
clone per file. JAWM controls local/Slurm execution and dependencies normally.
`jawm_rnaseq` already does this beside topGO after its annotator completes.
Pass `pipeline_file` to record the calling pipeline's Git revision independently
of the `jawm_string` revision. Production pipelines should load a published
`jawm_string@<commit>`; an unpinned name is useful while establishing the new repo.

### Run independently on an existing dataset

Edit `yaml/docker.yaml` to set `project_folder`, `organism` and optional `release`:

```sh
# Docker
jawm string_networks.py string -p ./yaml/docker.yaml

# Slurm and Apptainer with site-specific image paths
jawm string_networks.py string \
  -p ./yaml/docker.yaml ./yaml/hpc.yaml ./yaml/nexus.apptainer.yaml
```

The HPC profile follows the other modules' `/nexus:/nexus` bind and disables
automatic mounts. Place the module, input, output and cache paths under `/nexus`,
or adapt the bind settings to your site. `yaml/build.yaml` supplies the same
optional build-machine overrides used by `jawm_demo` and `jawm_template`:

```sh
jawm-test --python_versions system --jawm_versions local -y ./yaml/build.yaml
```

The standalone workflow discovers all
`<project_folder>/deseq2_output/annotated/*.results.tsv` files. Set a global
`input_file` to a single full path to process one contrast. Input tables require
`ensembl_gene_id`, `log2FoldChange` and `padj`; available expression columns are
carried into the network. An optional `model_formula` records a known source model;
standalone runs do not assume a model formula. Sample annotations are read from the contrast's
`deseq2_output/<contrast>.input.tsv`, falling back to `samples_MasterTable.txt`.

### Defaults and outputs

| Setting | Default |
| --- | --- |
| STRING database | 12.5, fixed version endpoint |
| Taxonomy | Dataset taxon ID, or exact scientific organism name in the bundled catalog |
| Differential-expression selection | padj < 0.05; minimum absolute log2FC 0 |
| Seed limit | 1,000 |
| Ranking | abs(log2FC × −log10(padj)) |
| Network | Functional associations |
| Combined confidence | ≥ 0.7 |
| Expansion | Up to floor(25% × unique mapped seeds) |
| Additional evidence gates | None |
| Saved layout | ForceAtlas2, 60 iterations, random seed 42 |

Unmapped top genes are reported rather than replaced. Expansion is one step,
ranked by distinct seed connections, then summed interaction confidence. Added
proteins retain available statistics from the full contrast; they need not be DE.
Isolated proteins are removed and reported. Limits, strict species checking and
API-error behavior match the app. No expression or fold-change data is sent to
STRING, only identifiers, taxonomy and network query parameters.

Each contrast writes to:

```text
<project_folder>/metadata/string-networks/<stable-network-id>/
    network.graphml               # Import into Cytoscape
    network.json                  # RNA-seq app viewer snapshot
    nodes.tsv / edges.tsv
    manifest.json                 # Contrast, settings, runtime, source hashes
    inputs.tsv / samples.tsv
    string_network.py             # Executed job script
    string_networks.py             # Module source
    reference/
    ... selection, mapping, expansion and isolate reports
```

The app discovers these completed networks after dataset upload or `.studio`
import. `report_files(project_folder)` returns GraphML and manifest paths for a
calling pipeline's report. Completed outputs are reused only if input, annotations,
module/catalog/runtime metadata and artifact checksums match. Changed inputs
rebuild with cached STRING responses; failed refreshes retain the last completed
network. Imported source files are provenance, never automatically executed.

### Shared cache and API coordination

The default cache is `<project_folder>/../.jawm-string-cache/`, outside dataset
exports and shared across neighboring datasets. Override it globally:

```yaml
- scope: global
  var:
    string_cache: /shared/reference-cache/string
```

All job nodes must see the same persistent directory, with working cross-host
POSIX `flock`. A global advisory lock covers cache lookup, missing HTTPS calls,
retries and atomic cache publication. Thus concurrent jobs cannot duplicate the
same miss or issue simultaneous requests; at least one second elapses before each
uncached attempt. Network calculation runs outside this shared lock.

Both exact queries and individual gene mappings / protein partner lists are
cached. Confirmed unmapped genes are cached as well; failures are not. Overlapping
comparisons reuse entities even when their batches differ. Keys include STRING
version, taxon, method and API parameters, keeping species, confidence thresholds
and network types separate. Misses are sent in batches of at most 500 identifiers.
Only misses require HTTPS access to `version-12-5.string-db.org`.

The species reference is unchanged STRING Consortium data, gzip compressed,
licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Source,
checksum and attribution are in `reference/species.v12.5.source.json`. Workflow
snapshot provenance is in `reference/workflow-source.json`.

### Tests

The standard test runs the production network process on real public airway
RNA-seq measurements and calls STRING 12.5. The fixture contains 80 genes
(64 DE and 16 non-DE), eight samples, and their measured DESeq2 statistics and
normalized expression. It is a subset of the
[public airway analysis](https://github.com/Bioinformatikci/rna-seq-workflow/tree/43eebf34e5ae6c088b03beba47b97c5502778a35).
Source revision, checksums, selection rule and license are included in the archive.

```sh
jawm-test --python_versions system --jawm_versions local
# The same test workflow directly:
jawm string_networks.py test -p ./test/yaml/test.yaml
```

Input is staged under `test/test-input/airway/`. The upload-ready archive is
`test/test-input/string_airway_test.tar.gz`; its commented MD5 entry is in
`test/data.txt`. Once it is uploaded to Figshare, insert the direct download URL
and uncomment that entry. JAWM will download, verify and extract it, following the
other modules. Until then, a new checkout needs the local archive extracted in
`test/test-input`.

The real GraphML, JSON, node/edge tables and manifest are written under
`test/test-output/metadata/string-networks/<network-id>/`. The test checks that the
network is nonempty and that these formats agree. `test_nodes.tsv` and
`test_edges.tsv` provide sorted output for reference-hash comparison, including
measured expression and STRING confidence/evidence scores. Layout coordinates,
timestamps and machine-specific provenance are excluded from this comparison.
No mocked responses are used by this integration test. Initial runs need HTTPS
access to `version-12-5.string-db.org`; subsequent runs reuse the persistent cache.
Public fixture identifiers and taxonomy are sent to STRING, while expression
values remain local.

The synthetic scientific/cache suite remains available separately:

```sh
jawm string_networks.py unit_test -p ./test/yaml/unit_test.yaml
# Or explicitly disable networking for the scientific unit tests:
docker run --rm --network none -v "$PWD:/module:ro" -w /tmp \
  mpgagebioinformatics/networkx:3.5 \
  python -m unittest discover -s /module/test -p test_string_module.py -v
```

See `test/README.md` for the fixture paths, download manifest and hash behavior.
GitHub Actions calls the shared `modules.webhook.yaml` on a self-hosted runner,
which runs the same real-input `jawm-test` entry. The runner needs the runtime
image and, until Figshare is configured, the staged input files.


## Optional shared local database

Set `STRINGDB` to opt into offline mapping and network queries. The same setting
works in the standalone STRING workflow and in the RNA-seq workflow:

```bash
export STRINGDB=/nexus/posix0/MAGE-flaski/service/databases/string
```

Alternatively, add it to the existing global `var` block in the project YAML:

```yaml
- scope: global
  var:
    STRINGDB: /nexus/posix0/MAGE-flaski/service/databases/string
```

A nonempty YAML/CLI variable takes precedence over the environment. The JAWM CLI
also accepts `--global.var.STRINGDB=/nexus/posix0/MAGE-flaski/service/databases/string`.
If neither supplies a path, the existing shared API cache is used. Unset the
variable in both places to return to online mode.

Before submitting contrast networks, `string_database` downloads species-specific
STRING **12.5** protein info, aliases, and detailed functional links, validates
the gzip files and schemas, and builds an indexed SQLite database. Files live in
`STRINGDB/12.5/`, prefixed with the taxon ID (e.g. `6239.functional.sqlite`). Other
species share the version directory without mixing their records. The backend
also supports the physical detailed-links format when requested by a caller.
Published manifests record source URLs, SHA-256 hashes, file sizes and the mapping
policy. Network provenance and restart fingerprints distinguish API and local data.

`STRINGDB/12.5/download.lock` is acquired exclusively across hosts. Contending
workflows wait, then validate the database before continuing. Downloads and indexes
are built in a temporary directory and the readiness manifest is published last.
An ordinary failure releases the lock and leaves no new ready database. After a
hard kill, the lock is deliberately retained: inspect its hostname and PID, confirm
that its owner has stopped, and only then remove the stale lock. Waiters time out
after 24 hours with an explanatory error; they never steal a live lock. Workers
also check readiness, including when a scheduler skips a previously completed
preparation step. The shared filesystem must support atomic exclusive creation
and rename. All jobs need write access for the lock, and the database root must be
visible inside containers (`mk.string_db` supplies automatic mounts; the Nexus
configuration explicitly binds `/nexus`).

Local mapping matches exact aliases. Ambiguities prefer an exact preferred-name
match, then the lexically first STRING protein ID. This reproducible rule may
differ from the API's best-match ranking. Missing IDs are omitted normally. Links
are deduplicated as undirected pairs and integer scores are divided by 1000;
selection, expansion, confidence and evidence filtering still use the network
worker's existing logic. No gene identifiers are sent to the API in local mode.

The site Apptainer configuration includes `string_database` alongside the network
and test jobs. When using these local edits, invoke this checkout's `rnaseq.py`
rather than fetching a fresh `jawm_rnaseq@latest-tag`; these changes have not been
published upstream.

Run synthetic tests (no bulk downloads) in the runtime image:

```bash
python -m unittest discover -s test -p 'test_*.py' -v
```
