# CarveMe Nextflow Pipeline

This repository runs [CarveMe](https://carveme.readthedocs.io/) on one or more protein FASTA files using Nextflow and an Apptainer or Singularity container.

Each input genome can be processed under multiple conditions in the same run:

- without gap filling;
- with a medium included in CarveMe, such as M9; or
- with a custom medium supplied in a media database.

The workflow currently builds metabolic models. Flux balance analysis (FBA) can be added as a separate module later, allowing the FBA tool to be changed without restructuring the CarveMe stage.

## Requirements

- Linux
- [Nextflow](https://www.nextflow.io/) 24 or later
- Java 17 or later, as required by your Nextflow version
- [Apptainer](https://apptainer.org/) or Singularity
- GNU Make, if using the supplied build and test commands

The input files must be protein FASTA files (`.faa`) containing translated genes. Raw nucleotide genome assemblies are not suitable CarveMe inputs.

## Repository layout

```text
.
├── assets/
│   └── media/                 Custom CarveMe media databases
├── containers/
│   ├── carveme.def            Container definition
│   └── carveme-1.6.6.sif      Built container, if present
├── modules/
│   └── local/
│       └── run_carveme.nf     CarveMe processes
├── test-params/               Parameter files used by tests
├── tests/                     Automated test scripts
├── main.nf                    Workflow entry point
├── nextflow.config            Defaults, profiles, and resources
└── Makefile                   Container and test targets
```

## Quick start

### 1. Obtain or build the container

If `containers/carveme-1.6.6.sif` is included with the repository, no build is needed. Otherwise, build it from the definition file:

```bash
apptainer build --fakeroot \
    containers/carveme-1.6.6.sif \
    containers/carveme.def
```

Confirm that CarveMe starts:

```bash
apptainer exec containers/carveme-1.6.6.sif carve --help
```

Messages about `squashfuse`, `fuse2fs`, or conversion to a temporary sandbox are informational if the command otherwise succeeds.

### 2. Create a parameter file

Create a YAML file such as `params.yml`:

```yaml
input: "genomes/*.faa"
outdir: "results"
solver: "scip"

media_databases:
  NO_GAPFILL: null
  M9: null
  SDMM: "assets/media/SDMM_media.tsv"
```

This example creates three models for every input `.faa` file:

1. a model without gap filling;
2. a model gap-filled using CarveMe's built-in M9 medium; and
3. a model gap-filled using the custom SDMM media database.

### 3. Run the pipeline

With Apptainer:

```bash
nextflow run main.nf \
    -profile apptainer \
    -params-file params.yml
```

With Singularity:

```bash
nextflow run main.nf \
    -profile singularity \
    -params-file params.yml
```

Nextflow resumes cached work when the same command is run with `-resume`:

```bash
nextflow run main.nf \
    -profile apptainer \
    -params-file params.yml \
    -resume
```

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `input` | `genomes/*.faa` | Glob selecting protein FASTA input files. Quote globs so the shell does not expand them. |
| `outdir` | `results` | Directory in which models, logs, and pipeline reports are published. |
| `solver` | `scip` | Optimization solver passed to CarveMe. The selected solver must be installed in the container. |
| `media_databases` | Repository defaults | Map of run-condition names to custom media database paths or `null`. See below. |
| `carveme_container` | `containers/carveme-1.6.6.sif` | Path to the local CarveMe SIF image. |

Defaults belong in `nextflow.config`. Use a YAML parameter file for settings that vary between runs. Command-line values can override either:

```bash
nextflow run main.nf \
    -profile apptainer \
    -params-file params.yml \
    --input "/data/proteins/*.faa" \
    --outdir results/experiment-01 \
    --solver scip
```

Nextflow parameters use two hyphens (`--input`), whereas Nextflow runtime options use one hyphen (`-profile`, `-params-file`, and `-resume`).

## Selecting media conditions

`media_databases` is a YAML mapping. Each key identifies a condition to run, and each value tells the workflow whether a custom database is required.

### No gap filling

Use the reserved `NO_GAPFILL` key:

```yaml
media_databases:
  NO_GAPFILL: null
```

### Built-in CarveMe medium

Set the value to `null` for a medium already known to CarveMe:

```yaml
media_databases:
  M9: null
```

The key must be a valid built-in CarveMe medium name.

### Custom medium

Set the value to the path of its media database:

```yaml
media_databases:
  SDMM: "assets/SDMM_media.tsv"
```

The mapping key must match the medium identifier defined in the database. Custom databases should be tab-separated and follow the format required by CarveMe, including BiGG metabolite identifiers.

An absolute path is safest when the parameter file is stored outside the repository:

```yaml
media_databases:
  SDMM: "/absolute/path/to/SDMM_media.tsv"
```

### Multiple conditions

List any combination of conditions in one mapping:

```yaml
media_databases:
  NO_GAPFILL: null
  M9: null
  SDMM: "assets/SDMM_media.tsv"
  ANOTHER_MEDIUM: "/data/media/another_medium.tsv"
```

The workflow creates one independent CarveMe task for every input genome and condition.

## Output

Published results are organized as:

```text
results/
├── models/
│   └── *.xml
├── logs/
│   └── *.log
└── pipeline_info/
    ├── dag.html
    ├── report.html
    ├── timeline.html
    └── trace.tsv
```

Models and logs may be published as symbolic links into Nextflow's `work/` directory. Do not delete `work/` while those links are needed. For durable archived results, change the relevant `publishDir` mode in `nextflow.config` from `symlink` to `copy`.

Use a new `outdir` for each run, or enable overwrite for the trace, report, timeline, and DAG files if reusing an output directory.

## Testing

Run the automated stub tests:

```bash
make test
```

Stub tests exercise workflow wiring and output publication without running the full CarveMe computation. They should cover no-gap-fill, built-in-medium, custom-medium, and combined configurations.

If the repository provides an integration-test target, run a real one-genome CarveMe test with:

```bash
make test-integration
```

You can also run an individual parameter file manually in stub mode:

```bash
nextflow run main.nf \
    -profile apptainer \
    -params-file test-params/no-gapfill.yml \
    -stub-run
```

## Troubleshooting

### Nextflow tries to pull `docker://${projectDir}/...`

The container path was interpreted as a remote image reference. Ensure `carveme_container` resolves to a local absolute path, for example in `nextflow.config`:

```groovy
carveme_container = projectDir.resolve('containers/carveme-1.6.6.sif').toString()
```

### A pipeline report or trace already exists

Use a different `outdir`, remove the old test output intentionally, or enable overwrite in `nextflow.config`:

```groovy
trace.overwrite = true
report.overwrite = true
timeline.overwrite = true
dag.overwrite = true
```

### Published model links exist but a test reports zero models

When outputs use `publishDir mode: 'symlink'`, tests must count symbolic links as well as regular files. Also confirm that the test is checking the configured `outdir`, usually its `models/` subdirectory.

### Custom media cannot be found

Check that:

- the database path is correct;
- the file is visible inside the container through Apptainer's bind mounts;
- the file is tab-separated;
- the mapping key matches the medium identifier in the database; and
- all compound identifiers are valid for CarveMe.

For additional runtime details, inspect `.nextflow.log`, the task log under `results/logs/`, and the corresponding task directory under `work/`.

## Extending the workflow

The CarveMe commands are isolated in `modules/local/run_carveme.nf`, while `main.nf` creates and connects channels. Additional modeling steps can be added as separate modules that consume the generated XML models.

Keeping FBA in its own module makes it possible to replace the FBA tool, container, resource requirements, or command without changing genome discovery or CarveMe model construction.
