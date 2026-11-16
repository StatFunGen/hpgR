# Course notebooks

Install hpgR from the repository root, then open a notebook with the R kernel.
Run cells in order.

Teaching examples and data downloads are listed in `catalog.json`.
`data/` contains their inputs. Models and estimators live in hpgR; figures
and plotting helpers live in the notebooks.

## Lecture 9

1. [Population Size and History from Present-Day DNA](P8149_Lecture_9_Effective_Size_and_History.ipynb)
2. [Trace the Sample through Its Parents](P8149_Lecture_9_Two_DNA_Copies.ipynb)
3. [Mutations on Trees and Population Size History](P8149_Lecture_9_Mutations_on_a_Tree.ipynb)
4. [Simulate and Date Local Genealogies](P8149_Lecture_9_Ancestral_Graphs.ipynb) accompanies the material after the homework.

## Lecture 10

1. [Split, Drift, Migration, and Differentiation](P8149_Lecture_10_Split_Drift_Migration.ipynb)
2. [Manhattan Rats and a Barrier to Migration](P8149_Lecture_10_Manhattan_Rats_and_Migration_Barriers.ipynb)
3. [Archaic Gene Flow and Shared Alleles](P8149_Lecture_10_Archaic_Gene_Flow.ipynb)
4. [Local Genealogies and Introgressed Segments](P8149_Lecture_10_Ancestral_Graphs.ipynb)

The Lecture 9 panel and the local genealogy notebook use the same simulated
population under `data/population_history/`.

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
