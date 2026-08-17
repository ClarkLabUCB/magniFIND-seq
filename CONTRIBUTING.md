# Contributing

Contributions should preserve the three public analysis areas—preprocessing,
fusion analysis, and downstream analysis—avoid site-specific paths or scheduler
directives, and include a regression test for behavioral changes.

## Development setup

```bash
conda env create -f environment.yml
conda activate magnifind-seq-analysis
```

## Before opening a pull request

```bash
make test
```

Keep fixture data minimal and document whether it is synthetic or derived.
Never commit identifiable instrument read names, unrestricted human sequence,
credentials, absolute private paths, production BAMs, or production outputs.

Update `docs/methods.md` and the applicable analysis-area README when parameters,
filters, output semantics, or scientific interpretation change.
