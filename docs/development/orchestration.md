# Workflow orchestration

[Documentation](../README.md)

Make owns developer commands, environment setup, Docker builds, checks, and
frontend tools. Nextflow owns scientific dependencies, per-type and per-gene
tasks, task caching, and publication. Python modules own the calculations.

## Entry points

`main.nf` is the shared entry point. It selects named workflows under
`workflows/` through `--mode`:

| Mode | Make target | Stages |
|---|---|---|
| `harmonization` | `make run` | Study preparation, harmonization, pseudobulk/optional single-cell merges, QC and provenance. |
| `analysis` | `make run-analysis` | DE/diagnostics, top-gene profiles and cell-type analysis/reports from existing merges. |
| `full` (default) | `make run-all` | Harmonization, downstream analyses and integration benchmarking. |

Analysis mode requires both `--merged_input` and `--pseudobulk_input`.
Make supplies their usual locations under `OUTDIR`; it checks that they exist
and does not rebuild core processing. Full mode always creates the single-cell
merge required by its downstream branches, even when the optional merge is
disabled in configuration. Harmonization mode retains the existing merge setting.
Unknown modes and missing analysis inputs fail before task submission.

```mermaid
flowchart LR
    Core[Core harmonization and merges] --> DE[DE and diagnostics]
    Core --> Cell[Cell-type split and analysis]
    DE --> Genes[Profile model preparation and per-gene tasks]
    DE --> Reports[Cell-type reports]
    Cell --> Reports
    Core --> Clocks[Optional age clocks]
    Core --> Benchmark[Integration benchmark in full mode]
```

Full mode includes integration benchmarking by default. Its methods retain
`maxForks = 1`, so only one benchmark method runs at a time; other analysis
branches can still run concurrently. This does not impose a RAM limit or
reserve the whole host for benchmarking. For constrained production runs,
explicitly skip it and run the focused benchmark later:

```bash
make run-all RUN_BENCHMARK=false
make run-integration-benchmark
```

Analysis mode skips benchmarking by default; enable it with
`RUN_BENCHMARK=true`. Age clocks remain opt-in in either downstream mode with
`RUN_CLOCKS=true`. Harmonization mode invokes neither branch.
The corresponding Nextflow switches are `--run_benchmark true|false` and
`--run_clocks true|false`.

`main.nf --mode full` replaces the former `pipeline.nf` entry point;
`main.nf --mode analysis` replaces `pipeline.nf --analysis_only`.

The standalone `.nf` entry scripts use these same named workflows. They accept
existing inputs for focused work; targeted report commands never acquire an
implicit prerequisite on the core workflow. Production runs are resumable;
test aliases run fresh. Production/test Make targets share recipes and derive
input and output paths from `OUTDIR`, unless an input path is explicitly supplied.

## Configuration preflight

Each scientific entry point validates its consumed settings before task fan-out.
`VALIDATE_CONFIGURATION` reuses the Python configuration validators, produces an
unchanged copy of `pipeline.json`, and supplies that file to downstream tasks.
In full/harmonization mode, the preflight also checks the study and input-source
registries and schemas; even RDS conversion waits for it. Clock runs validate
method IDs, search grids, worker counts, external-study IDs and calibration options
before fitting. Top-gene counts must be positive integers before dispatch. Eligibility that depends on actual samples is still assessed by
the scientific task.

A preflight failure is fatal. Independent scientific/report task failures retain
the existing partial-run policy after validation succeeds. Preflight itself neither rebuilds core nor refits scientific inputs. Focused
report commands continue to reuse prepared artifacts when they remain valid.
`make workflow-test-preflight` uses the locked `.venv` and real validators to test
valid configuration, invalid analysis/clock settings and malformed core registries
without fitting or notebook execution. CI's unit and main-workflow stub tests cover
the validation rules and fatal dependency barrier separately.

The staged configuration path changes with this addition, so the first run after
migration can recompute tasks whose cache keys include that path. Configuration
bytes, scientific settings and published result paths are preserved.

## Dependencies and caching

Composed stages consume producer output channels. Gene ranking uses the current
run's DE results, and cell-type reports receive the current run's matching DE
files. Published output directories are not used to discover intermediate inputs
inside the composed graph. Standalone commands intentionally accept existing
published artifacts.

Gene profiles prepare one immutable PyDESeq2 model per cell type, then calculate
and render individual genes in separate tasks. A prepared model must match the
configuration, counts, samples, genes, design, and software version; a mismatch
fails the gene task rather than silently refitting. Per-type profile context is
collected independently of cross-type DE reporting, so an aggregation failure
does not prevent successful types from producing gene profiles.

Clock fitting runs per cell type, followed by independent method/comparison
reporting and a cross-type aggregation of this run's successful fits. Eligible
skips remain explicit records. Clock publication is confined to `ageing_clock/`;
it must not replace a whole cell-type directory shared with other analyses.

Moving a process into a named workflow changes its task name and can invalidate
old cache entries. The public focused entry points remain available, but the
first run after this orchestration migration may recompute tasks. Keep Nextflow's
cache and work directories for subsequent resumes. `make workflow-test-resume`
checks that two genes share one preparation task, identical runs reuse all tasks,
and changing only the report template reruns rendering alone.

## Failure isolation

Independent analysis and report task failures do not stop sibling tasks or
unrelated branches. Successful outputs remain published, including completed cell-type analysis
files when their report fails. Dependent tasks that
lack an upstream result cannot run; aggregators use completed streams and
explicit expected inputs to describe partial results, including the case where
all cell-type analyses fail. Eligibility skips are distinct from unexpected
clock failures. Essential input/configuration errors and core processing failures
remain fatal.

Every Make workflow invocation goes through `scripts/run_workflow.py`. It writes:

- `<outdir>/provenance/workflows/<target>-<run-id>/trace.tsv`: task status and exit codes.
- The adjacent `nextflow.log` and `summary.json`: `complete`, `partial`, or `failed`, with failed tasks
  and declared application-level partial summaries.
- `<outdir>/workflow_summary.json`: an atomic copy of the latest invocation summary.

An analysis run with ignored task failures exits successfully and records
`partial`. Verification checks require `complete` and fail after the workflow
has finished; old published artifacts cannot make a partial current run pass.

The runner prints failed task names, exit codes, error messages and log paths
after completion. It also lists failed DE models from current per-type metadata
and the collector manifest, including their cell type, target covariate, sample
count and fitting error. This exposes model failures inside otherwise successful
tasks. Summaries retain these details in the additive `task_diagnostics` and
`analysis_failures` fields, with bounded task-log excerpts.

Inspect the latest saved run without executing tasks or modifying outputs:

```bash
make workflow-failures OUTDIR=output
```

For an older run, pass its archived summary:

```bash
python3 -m scripts.show_workflow_failures --summary output/provenance/workflows/<run>/summary.json
```

Older summaries without saved diagnostics can still be inspected while their
work directories exist. Console output lists up to ten failed tasks and models;
new JSON summaries retain the full diagnostic lists. Nextflow also supports inspecting its execution
history directly, using the run name or session ID from the relevant invocation:

```bash
nextflow log
nextflow log <run-name> -f name,status,exit,workdir -F 'status == "FAILED" || status == "ABORTED"'
```

The [Nextflow log reference](https://docs.seqera.io/nextflow/reference/cli/log/)
also documents the `stderr` field. In a task's work directory, `.command.err`
contains standard error, `.command.out` standard output and `.command.sh` the
executed task script. Nextflow sees task exits; the repository's model summary
adds scientific failures caught and recorded inside a task.

To publish successful DE results and continue combined reporting after an
individual fit fails, opt into partial collection:

```bash
make run-all ALLOW_PARTIAL_RESULTS=true
```

The failed fit remains marked `failed`, and the DE and workflow summaries remain
`partial`. No model design, inclusion filter or numerical fitting rule changes.
The default remains strict; verification rejects partial results. A run with no
completed DE models still fails collection even when partial results are allowed.

`make workflow-test-failures` injects fitting failures, rendering failures, and
an all-failed cell-type branch into the production process definitions. It
asserts continued sibling work, retained scientific artifacts, and partial
summaries. `make workflow-test-main-failures` additionally invokes the actual
`main.nf` entry point with one failed DE fit, DE report, or benchmark method per
case. It checks that unrelated branches finish, cell reports consume current DE
fits despite failed DE rendering, and queued benchmark methods continue after
the first failure. CI runs these checks without scientific inputs or Docker.

Failure-injection cases announce `[EXPECTED FAILURE TEST]` and describe the
deliberate failure before starting. Their raw Nextflow output and detailed
failure summaries are saved to the case's printed `console.log` path. A
`[PASS]` line appears only after the harness checks the expected task exits,
failure isolation and required outputs. An unexpected result prints `[FAIL]`,
shows the diagnostic log tail and fails the check. Invalid-configuration cases
use the same labels and verify that validation blocks scientific tasks.
These labels describe test expectations; saved workflow states and exit codes
retain their actual failed or partial values. Production error reporting is
unchanged.

## Integration tests and finalization

`make workflow-test-modes` previews the production entry point with placeholder
inputs. It checks every mode, benchmark defaults and overrides, optional clocks,
single-cell merge requirements, and preflight failures without executing tasks.

`make run-all-test` exercises full mode with a fresh core run, then
checks its scientific artifacts. The current trace must contain the required
branches; selected gene requests must have matching fit/render tasks and published
artifacts, with one model preparation per type. These checks reuse the existing
run and do not launch additional calculations. It creates the synthetic pseudobulk once and
shares it across DE, clocks, top-gene profiles, and an additional planted-gene
profile. The benchmark and three cell-type reports use the real test merge.

`make run-analysis-test` uses analysis mode with benchmarking and clocks enabled
to run that same downstream graph against existing test
merges. Use it while changing orchestration without rebuilding the core.
`VERIFY_CELL_TYPES` selects its cell types; the verification clock defaults to
`elastic_net` with the small configured search grid.

The shared runner refreshes the workspace report index once after Nextflow
finishes, including partial and failed runs. This is deliberately an uncached
workspace inventory, so it can also show existing published reports. Use
`make render-report-index` to refresh it manually. Frontend release export and
builds remain separate explicit operations.

## Sequential clocks and resuming

Clock fitting uses one Nextflow task per cell type and `maxForks = 1`: fits do
not overlap, although their execution order is not guaranteed. Methods within a
fit use the configured internal worker count; report tasks and unrelated branches
may overlap fitting. This does not reserve the entire machine for clocks.

Production Make targets pass `-resume`; test runs are fresh. Matching successful
fit/render tasks are reused when their inputs, scripts, image and workflow identity
match and their cached outputs remain available. Failed tasks are retried on a new
resumed invocation. Preserve both `.nextflow/cache/` and the work directory;
published files alone do not recreate the task cache. A template-only change can
rerun rendering while retaining fitting when the other fit inputs remain unchanged.

Bare `-resume` selects the most recent run, which may be a different target. After
intervening runs, use the original `session_id` from its archived workflow summary:

```bash
make run-age-clocks NF_RESUME="-resume <session-id>"
```

Reusing the same entry point and input paths is the reliable route to reuse.
A focused command reading a published merge can have different cache keys from
a composed run reading a producer's work-directory artifact.
