
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


### Example metabolites table

Fields are tab-separated:

    totalarea 	1.3000	0.9000	1.2797	0	1	
    totalarea1 	4.2	1.3000	1.3246	0	1
    totalarea2	9.0	5.00	5.2431	0	1
    glucose		5.2500	5.2200	5.2367	1	1
    lactate		4.1137	4.0920	4.1085	1	1
    citrate		2.5600	2.5400	2.55	1	1
    pyruvate	2.3780	2.3650	2.373	1	1
    ethanol		1.2028	1.1926	1.1969	0	0
    ethanol1	1.1913	1.1798	1.1885	0	0
    glycine		3.5700	3.5583	3.563	1	1
    acetate		1.9231	1.9165	1.9189	1	1
    lypids		1.3150	1.2500	1.278	1	1
    glycerol	3.6959	3.6368	3.6674	0	0

