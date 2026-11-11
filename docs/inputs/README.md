# Acquiring inputs

[Documentation](../README.md)

## Download inputs

> **Strictly opt-in and potentially large.** Only `make download-inputs`,
> `make download-gene-set-inputs`, and `make download-inputs-strict` retrieve data; `make install`, test targets,
> and production workflow targets never do. Before it starts, the downloader
> warns about pending artifacts and their recorded size where available. A
> selected study can require multi-gigabyte source files.

Specify the studies to acquire. The best-effort command downloads configured
public URLs without portal authentication or scraping and records outcomes in
ignored `input_data/download_manifest.json`. It returns a failure if a direct
download fails or fails its size/checksum check; manual-only inputs remain
reported for best-effort acquisition:

```bash
make download-inputs STUDIES=aifi
make download-inputs STUDIES=nehar_belaid26
make download-gene-set-inputs
```

Study selections include the small pipeline-wide GO enrichment inputs in
`input_data/go/`. Use `make download-gene-set-inputs` to retrieve those files and
the static atlas example pathways, without downloading study data. Study
selections do not download the atlas pathways. GO source files can change upstream, so the
registry intentionally has no fixed checksum or release requirement; the
download manifest and each report's enrichment provenance record the exact
SHA-256 values used, along with the GO ontology `data-version`.

The downloader prints each artifact's position in the run and its associated
study name. During a transfer it also shows bytes downloaded, transfer rate,
and percentage when the source reports a size or the manifest records one;
otherwise it reports downloaded bytes and rate.
Interrupted direct downloads retain a compatible `.part` file and transfer
metadata. The downloader makes up to three attempts with delay and resume
messages, and resumes with HTTP byte ranges only when the source confirms that
the saved partial is compatible. A final size and checksum check is still
required before publishing the requested file.

`make download-inputs-strict STUDIES=fachrul26` is the explicit acquisition
preflight. In addition to failed downloads and size/checksum mismatches, it
exits nonzero for a manual-only artifact that has not yet been placed at its
configured path. It is an acquisition preflight, not a routine test, and never
runs as part of analysis commands. Run it before production when download
integrity matters. Existing files are left unchanged unless
`DOWNLOAD_FORCE=true` is supplied deliberately.

This preflight can check only integrity values recorded for each artifact. A
manual-only source with no expected size or checksum is checked for presence
only; for example, the local Terekhova23 H5AD is a development subsample and
can pass that presence check. See [INPUT_FILES.md](files.md) for the
development-versus-production input distinction and per-source details.

The manifest records manual-only sources and their acquisition instructions.
Synapse expression objects require portal access, the AIDA publisher
attachment rejects unauthenticated programmatic requests, and the AIFI full
H5AD link currently redirects through Allen's Google-authenticated workspace.
The downloader does not implement these access flows. Entries explicitly marked
as allowing a local subset retain their published full-file size as download
documentation; a normal download check reports `local_subset` rather than a
size failure, while the strict preflight still rejects it. Confirm that acquired
files are suitable for production analysis.
See `config/input_sources.json` and [INPUT_FILES.md](files.md).
