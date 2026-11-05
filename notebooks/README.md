# Course notebooks

Install hpgR from the repository root, then open a notebook with the R kernel.
Run cells in order. Student assignments contain the task and supplied code;
complete the remaining cells yourself.

Teaching examples and student assignments are listed in `catalog.json`.
`data/` contains their inputs. Models and estimators live in hpgR; figures
and plotting helpers live in the notebooks.

Commit notebooks and their inputs. HTML exports and ZIP downloads are local
generated files; the repository keeps the RDS inputs directly.

From the package root, check the notebooks with:

```sh
Rscript tools/sync_vignettes.R
```

This runs supplied R cells in fresh sessions using the current package.
The website uses a separate publication list. Answers and examination
materials are maintained in the instructor repository.
