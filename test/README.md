# jawm tests

The standard test runs the real STRING network workflow on a subset of the
[public airway RNA-seq analysis](https://github.com/Bioinformatikci/rna-seq-workflow/tree/43eebf34e5ae6c088b03beba47b97c5502778a35),
using 80 genes with measured statistics and normalized expression from eight
samples. The source analysis uses DESeq2 with `~ cell + dex` and apeglm shrinkage.
The subset is a software fixture, not a full biological analysis.

```sh
jawm-test --python_versions system --jawm_versions local
# Equivalent workflow:
jawm string_networks.py test -p ./test/yaml/test.yaml
```

Docker, `mpgagebioinformatics/networkx:3.5`, and HTTPS access to the fixed STRING
12.5 endpoint are required on the first run. There are no mocked STRING responses.
Later runs can reuse `test/.jawm-string-cache` and completed output when unchanged.
To build the image locally:

```sh
docker build -t mpgagebioinformatics/networkx:3.5 -f containers/string/Dockerfile .
```

## Input archive and Figshare

Upload this file to Figshare:

```text
test/test-input/string_airway_test.tar.gz
MD5: 9d1f162e0baab22acb6de79640c1c13d
```

The archive contains:

```text
airway/
    README.md
    LICENSE
    source.json
    deseq2_output/
        dexamethasone_vs_untreated.input.tsv
        annotated/
            dexamethasone_vs_untreated.results.tsv
```

`source.json` records the source revision, original file checksums, precise
selection rule and transformations. Numeric measurements are unchanged from the
source. The results table joins actual normalized expression for all eight
samples; the sample table records the four matched cell-line pairs.

The archive and extracted files are currently staged in `test/test-input`, which
is ignored by Git as in the other modules. Until the upload is configured, a new
checkout must receive and extract that archive manually:

```sh
tar -xzf test/test-input/string_airway_test.tar.gz -C test/test-input
```

After upload, replace `<FIGSHARE_DOWNLOAD_URL>` in `data.txt` with the direct file
download URL and uncomment the prepared entry:

```text
9d1f162e0baab22acb6de79640c1c13d  string_airway_test.tar.gz  <FIGSHARE_DOWNLOAD_URL>
```

`jawm-test` reads the conventional `<md5sum> <filename> <url>` columns, downloads
into `test/test-input`, extracts the archive there, and verifies the MD5. No URL
has been fabricated; the line remains commented until the file is uploaded.

## Real outputs and reference hashes

`test/tests.txt` runs the `test` workflow with `test/yaml/test.yaml`. Outputs are
under:

```text
test/test-output/metadata/string-networks/2c3292c1aed65bf69ea0f0cee2b74a4e/
    network.graphml
    network.json
    nodes.tsv
    edges.tsv
    manifest.json
    inputs.tsv
    samples.tsv
    ... selection, mapping, expansion and isolate reports
    test_nodes.tsv
    test_edges.tsv
```

The test rejects an empty network, missing artifacts, inconsistent GraphML/JSON/
table sizes, missing measured expression, duplicate edges, and scores below the
requested confidence. The hash scope compares `test_nodes.tsv` and
`test_edges.tsv`: sorted real node identities, roles, measured fold changes,
adjusted p-values, expression, and STRING confidence/evidence scores. Floats are
formatted consistently; coordinates, timestamps, local paths and cache statistics
are omitted. A network biology change therefore changes the reference hash.

`jawm-test` compares this hash with the recorded value in `tests.txt`. A blank
reference is populated after a successful first run. Use `jawm-test -o` only after
reviewing an intentional change. STRING versioned data can still change upstream;
review any resulting network difference before updating a reference.

## Separate synthetic unit tests

The prior synthetic scientific/cache tests remain available independently:

```sh
jawm string_networks.py unit_test -p ./test/yaml/unit_test.yaml
# Explicitly offline:
docker run --rm --network none -v "$PWD:/module:ro" -w /tmp \
  mpgagebioinformatics/networkx:3.5 \
  python -m unittest discover -s /module/test -p test_string_module.py -v
```

Two controller checks additionally verify repository loading, process cloning and
rendered inline-worker execution, including paths containing spaces and quotes.
They run when JAWM is installed alongside the pinned scientific dependencies:

```sh
python -m unittest discover -s test -p test_string_module.py -v
```

The shared GitHub module workflow runs the real `jawm-test` entry. Its runner needs
the image and access to STRING, plus locally staged inputs until the Figshare
entry has been enabled.
