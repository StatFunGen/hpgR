# Course notebooks

Install hpgR from the repository root, then open a notebook with the R kernel.
Run cells in order.

Teaching examples and data downloads are listed in `catalog.json`.
`data/` contains their inputs. Models and estimators live in hpgR; figures
and plotting helpers live in the notebooks.

Assignment requirements and function suggestions are in the lecture slides.
Students write their own reports using hpgR and the supplied data.

Commit notebooks and their inputs. HTML exports and ZIP downloads are local
generated files; the repository keeps the RDS inputs directly.

From the package root, check the notebooks with:

```sh
Rscript tools/sync_vignettes.R
```

This runs supplied R cells in fresh sessions using the current package.
The `website` field in `catalog.json` selects published notebooks; its entry
order sets the reading order within each lecture. Build and
preview the site with:

```sh
cd notebooks
pixi install
pixi run build
pixi run preview
```

The build verifies saved code and input records, then copies selected notebooks
and their inputs into an isolated directory. It checks page links and downloads.
Answers and examination materials are maintained in the instructor repository.

The Pages workflow builds this same publication list on pushes to `main`.
It sets `BASE_URL` from the repository's Pages configuration and publishes only
the verified HTML directory. Pull requests run the build checks.
