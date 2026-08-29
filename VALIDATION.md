# Validation record

The repository package was validated on 28 August 2026 under Windows with R
4.5.x.

- The annotation workflow was rerun from the packaged compressed inputs.
- The downstream differential-expression, PLS, regression, GO-enrichment,
  figure, and supplementary-table workflow was rerun from the packaged fitted
  model checkpoints.
- The downstream CAFE summary, functional classification, node-wise GO,
  figure, and supplementary-table workflow was rerun from the packaged CAFE
  files.
- All three root entry scripts were tested with `rerun <- FALSE` from a working
  directory outside the repository and successfully loaded their compact
  workspaces.
- All packaged R scripts passed an R parser check.
- No packaged file is 100 MB or larger.

The computationally costly gene-wise and N13-HOG mixed-model refits are
provided as one-worker upstream scripts but were not repeated during this
final packaging check; their fitted checkpoints are included and were used by
the validated downstream workflow.

On 29 August 2026, the packaged input paths were migrated to the retained
`nextflow_runs/` directory. Byte-identical duplicate Salmon, EXCON GO,
OrthoFinder HOG and CAFE files were removed from the downstream input folders;
all R scripts were parsed again and the revised required-input paths were
checked against the repository contents.
