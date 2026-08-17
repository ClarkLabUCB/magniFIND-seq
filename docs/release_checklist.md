# Repository release checklist

Unknown author-, manuscript-, and archive-controlled facts are kept separately
in `docs/information_to_complete.md`. This checklist contains only operations
that must be completed before calling the repository a full reproduction
package.

- Resolve every applicable entry in `docs/information_to_complete.md` without
  substituting guessed values or software defaults.
- Run `make test` and `make example` in the exact `linux-64` lock environment
  through GitHub Actions. Archive the successful workflow URL and commit SHA.
- Run every configured study workflow through
  `03_downstream_analysis/code/run_all.sh` using complete
  production inputs. Preserve configurations, normalized manifests, input and
  source SHA-256 records, effective parameters, and R session information.
- Compare regenerated study outputs with every manuscript number and Figure
  6B–G source value. Investigate discrepancies; do not edit expected values to
  make a failing comparison pass.
- Review all files proposed for tracking, especially sequence-derived examples,
  for data-use permission, personal identifiers, credentials, absolute local
  paths, and unexpectedly large files.
- Create the initial commit, add the intended remote, tag a reviewed version,
  publish the GitHub release, and archive that exact tag if a DOI is desired.
