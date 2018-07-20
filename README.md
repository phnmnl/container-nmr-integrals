
# NMR Integrals

Version: 0.1.0


## Short Description

Compares specific metabolite levels in two NMR spectra of blood serum/plasma samples.


## Description

This tool implements a comparison between the NMR spectra of a reference sample and a test sample, based on selected metabolite signals.

Two aliquots of the same sample have the same metabolomic profiles, and thus have exactly the same NMR spectra. Based on this premise, the NMR spectrum can be used to check the quality of samples. 

This tool compares the NMR spectrum of a blood serum/plasma sample with the NMR spectrum of a different aliquot of the same sample acquired at a different time (e.g., spectrum acquired at the moment of the sample collection vs. spectrum acquired after long-term storage). The metabolites tested are those ones that were experimentally identified as "sensible" indicators of sample quality. The respective NMR signals are selected, and their relative concentrations inside each of the samples compared are measured. Spectral total area is also calculated to have information about the presence of contaminations. The operation is part of a QC procedure for biobanks, where samples are periodically analyzed to detected and measure long-term quality decay.

The tool generates a table and a set of plots indicating the signal differences and whether they are significant.

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
Usage: ./integrals.R [options] REFERENCE_SPECTRUM TEST_SPECTRUM


Options:
        --left=N
                Left

        --right=N
                Right

        --where=N
                Where

        --plotfile=FILE.PDF
                Plot file

        --metabolites
                Metabolites table in tab-separated format

        -h, --help
                Show this help message and exit

        REFERENCE_SPECTRUM
                Reference dataset (Bruker monodimensional 1H cpmg NMR spectrum)

        TEST_SPECTRUM
                Test dataset (Bruker monodimensional 1H cpmg NMR spectrum)
```


### Input


Left - Left margin of the integration range (ppm value)

Right - Right margin of the integration range (ppm value) (must be greater than left)

Where - ppm value (inside the integration rage) where we want to align the signals of the two spectra.

Reference spectrum - Bruker monodimensional 1H cpmg NMR spectrum

Test spectrum - Bruker monodimensional 1H cpmg NMR spectrum

Signal Metabolites - Optional table of signal metabolites, to override default behaviour.


#### Metabolites table


The metabolites table specifies, for each selected metabolite, an integration
window for one or more of its signals inside the spectra.  The program will by
default use a table (shown below) made specifically for the purpose of detecting
decay of blood serum/plasma samples.  You may however customize this table to
apply the technique to other types of samples.


Here is the table format.

Columns:

1st column - Metabolite name

2nd column - Right margin of the integration range (ppm value)

3rd column - Left margin of the integration range (ppm value)

4th column - ppm value (inside the integration rage) where we want to align the signals of the two spectra.

5th column - base line correction (0: no /1: yes)

6th column - signal alignment at where (0: no /1: yes)


The table below is the *default table*, which has been verified to return sensible
values:

    Metabolite    right ppm   left ppm   ppm align   baseline correct on/off   alignment on/off  

    totalarea       1.3000   0.9000      1.2797      0                         1                 
    totalarea1      4.2      1.3000      1.3246      0                         1                 
    totalarea2      9.0      5.00        5.2431      0                         1                 
    glucose         5.2500   5.2200      5.2367      1                         1                 
    lactate         4.1137   4.0920      4.1085      1                         1                 
    citrate         2.5600   2.5400      2.55        1                         1                 
    pyruvate        2.3780   2.3650      2.373       1                         1                 
    ethanol         1.2028   1.1926      1.1969      0                         0                 
    ethanol1        1.1913   1.1798      1.1885      0                         0                 
    glycine         3.5700   3.5583      3.563       1                         1                 
    acetate         1.9231   1.9165      1.9189      1                         1                 
    lypids          1.3150   1.2500      1.278       1                         1                 
    glycerol        3.6959   3.6368      3.6674      0                         0                 


#### Output

- Plots of ppm vs relative intensity, for each metabolite (PDF format).
- Table (tab-separated columns) with the difference in the relative
concentration of the selected metabolites. Threshold: +/-20%.

Table columns:

1st column - Metabolite name

2nd column - ok/bad

4th column - relative concetration difference (percentage)


The 3rd column is always "diff" and the 5th one is always a "%" sign.
