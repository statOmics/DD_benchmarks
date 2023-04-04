# Differential detection benchmark project


This repository contains the code associated with our differential detection benchmark paper.

## Installing dependencies

The `package/` directory contains a helper package to provide code used throughout the project. This
can be installed from the command line with

```sh
R CMD INSTALL package
```

In addition, a [renv](https://rstudio.github.io/renv/articles/renv.html) environment is provided,
specifying all dependencies and their versions to run the analyses. To install the dependencies, run
this in an R console:

```r
renv::restore()
```

## Downloading raw data

### Lupus data

The Lupus data from [Perez *et al.* (2021)](https://doi.org/10.1126/science.abf1970) can be
downloaded from GEO using accession number
[GSE174188](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174188).

Download the
`GSE174188_CLUES1_adjusted.h5ad.gz ` file and place it in the `data-raw/` directory under `benchmarks/lupus/`.

## Running the benchmarks

The `benchmarks/` folder contains the following benchmarks:

* `lupus/`: benchmarks on the Lupus data

Each subfolder is accompanied by R scripts in a `scripts/` folder and R Markdown files in the
`analysis/` folder. `Makefile`s are provided to control the order in which the scripts should be
run.

A master `Makefile` in the root of this repository is also provided to run each of the benchmarks.
From the command line at the root of this repository, run

```sh
make lupus
```

to run each of the respective benchmarks. Or simply

```sh
make benchmarks
```

to run them all in one go.

**Note** that this will take **a considerable amount of time** and computational resources. It is
recommended to run this on an HPC instance.
