#!/bin/bash

set -o errexit
set -o nounset

IntegralsScript="integrals.R"

if [[ $# != 9 ]]; then
  echo "Usage error!  Got $# arguments" >&2
  exit $?
fi

function error_trap() {
  echo "#### Error at line ${BASH_LINENO[0]} running command ${BASH_COMMAND} ####" >&2
}

trap error_trap ERR


WorkDir="${1}"
Left="${2}"
Right="${3}"
Where="${4}"
OutputPlotfile="${5}"
SignalMetabolitesTable="${6}"
ReferenceSpectrumZip="${7}"
TestSpectrumZip="${8}"
OutputTable="${9}"


reference_dir="${WorkDir}/input/reference"
test_dir="${WorkDir}/input/test"

mkdir --parents "${reference_dir}" "${test_dir}"

echo "Unzipping input archives" >&2
unzip -q -d "${reference_dir}" "${ReferenceSpectrumZip}"
unzip -q -d "${test_dir}" "${TestSpectrumZip}"

echo "Analyzing" >&2

"${IntegralsScript}" --left=${Left} --right=${Right} --where=${Where} --plotfile "${OutputPlotfile}" "${SignalMetabolitesTable}" "${reference_dir}" "${test_dir}" > "${OutputTable}"

echo "Finished" >&2
exit 0
