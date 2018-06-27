
# NMR Integrals

Version: 0.1.0


## Short Description

Compares specific metabolite levels in two NMR spectra.


## Description

This tool implements a comparison between the NMR spectra of a reference sample
and a test sample, based on selected metabolite signals.  The operation is part
of the QC procedure implemented in a number of biobanks, where samples are
periodically analyzed to detected and measure long-term quality decay.


The tool generates a table and a set of plots indicating the signal differences
and whether they are significant.


## Key features

- NMR Processing

## Instrument Data Types

- NMR


## Approaches
- Metabolomics / Targeted


## Tool Authors 

- Leonardo Tenori (CERM)
- Veronica Ghini (CERM)
- Antonio Rosato (CERM)
- Luca Pireddu (CRS4)


## Container Contributors

- Luca Pireddu (CRS4)

## Git Repository

- https://github.com/phnmnl/container-nmr-integrals.git

## Installation

For local individual installation:

```
docker pull docker-registry.phenomenal-h2020.eu/phnmnl/nmr-integrals
```

## Usage Instructions

For direct Docker usage:
```
docker run docker-registry.phenomenal-h2020.eu/phnmnl/nmr-integrals
```


```
Usage: ./integrals.R [options] METABOLITES_TABLE REFERENCE_SPECTRUM TEST_SPECTRUM


Options:
        --left=N
                Left

        --right=N
                Right

        --where=N
                Where

        --plotfile=FILE.PDF
                Plot file

        -h, --help
                Show this help message and exit

        METABOLITES_TABLE
                Metabolites table in tab-separated format
        REFERENCE_SPECTRUM
                Reference dataset (Bruker NMR)

        TEST_SPECTRUM
                Test dataset (Bruker NMR)
```
