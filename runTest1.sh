#!/bin/bash

function error_trap() {
  echo "#### Error at line ${BASH_LINENO[0]} running command ${BASH_COMMAND} ####" >&2
}

trap error_trap ERR

set -o errexit
set -o nounset

cd /nmr-integrals

RefZip=reference_dataset.zip
TestZip=test_dataset.zip
MetabolitesFile=metabolites.txt
Expected=expected_output.tsv
PlotFile=somefile.pdf
OutputTable=output_table.tsv

echo "Downloading data..." >&2
wget -O "${RefZip}" 'https://drive.google.com/uc?export=download&id=1Ywd4IMsZZHsiDg4NUcEos9Q475ac9vBP'
wget -O "${TestZip}" 'https://drive.google.com/uc?export=download&id=1OV77pvk5aE4eiyARHfA_vf_DeX-B5HG3'

wget -O "${MetabolitesFile}"  'https://drive.google.com/uc?export=download&id=1B9SoJN0vJOFtWEllgM24rdZ8kTJU06Ei'
wget -O "${Expected}" 'https://drive.google.com/uc?export=download&id=1g29hI7qTSfMqUv3lByF3nwkBRrUYa5S-'

echo "##### Running test 1 with metabolites file #######" >&2

integrals_galaxy_wrapper.sh "${PWD}" 5.257 5.225 5.243 "${PlotFile}" "${RefZip}" "${TestZip}" "${OutputTable}" "${MetabolitesFile}"

# for ascii-based sorting.  Line order isn't important in our table comparison
export LC_ALL=C
if ! diff <(sort "${Expected}") <(sort "${OutputTable}") ; then
  echo "UNEXPECTED DIFFERENCES IN OUTPUT TABLE" >&2
  echo "Test failed" >&2
  exit 2
elif [[ ! -f "${PlotFile}" ]] ; then
  echo "Plot file ${PlotFile} wasn't generated" >&2
  echo "Test failed" >&2
  exit 3
fi

echo "##### Running test 2 with default metabolites #######" >&2

integrals_galaxy_wrapper.sh "${PWD}" 5.257 5.225 5.243 "${PlotFile}" "${RefZip}" "${TestZip}" "${OutputTable}"

# for ascii-based sorting.  Line order isn't important in our table comparison
export LC_ALL=C
if ! diff <(sort "${Expected}") <(sort "${OutputTable}") ; then
  echo "UNEXPECTED DIFFERENCES IN OUTPUT TABLE" >&2
  echo "Test failed" >&2
  exit 2
elif [[ ! -f "${PlotFile}" ]] ; then
  echo "Plot file ${PlotFile} wasn't generated" >&2
  echo "Test failed" >&2
  exit 3
fi

echo "Test passed" >&2
exit 0
