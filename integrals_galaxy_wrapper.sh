#!/bin/bash

set -o errexit
set -o nounset

IntegralsScript="integrals.R"

if [[ $# != 8 && $# != 9 ]]; then
  echo "Usage error!  Got $# arguments" >&2
  exit 1
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
ReferenceSpectrumZip="${6}"
TestSpectrumZip="${7}"
OutputTable="${8}"
SignalMetabolitesTable="${9:-}"


reference_dir="${WorkDir}/input/reference"
test_dir="${WorkDir}/input/test"


function unzip_strip_components_1() {
  # unzip the dataset archive stripping the first path component if its the only one.
  # We assume that the user zipped the entire dataset tree.
  # If we find more than one path at the base of the archive tree then we move
  # everything to the input directory.
  local archive="${1}"
  local dest="${2}"

  local tmpdir=$(mktemp -d --tmpdir="${WorkDir}")
  unzip -d "${tmpdir}" "${archive}"

  mkdir --parents "${dest}"

  shopt -s dotglob
  local items=("${tmpdir}"/*)
  if [[ ${#items[@]} == 1 ]]; then
    mv "${tmpdir}"/*/* "${dest}"
  else
    mv "${tmpdir}"/* "${dest}"
  fi
  rmdir "${tmpdir}"/* "${tmpdir}"
  return 0
}

#### main ####
echo "Unzipping input archives" >&2

unzip_strip_components_1 "${ReferenceSpectrumZip}" "${reference_dir}"
unzip_strip_components_1 "${TestSpectrumZip}" "${test_dir}"

echo "Analyzing" >&2

if [[ -n "${SignalMetabolitesTable:-}" ]]; then
  "${IntegralsScript}" --left=${Left} --right=${Right} --where=${Where} --metabolites "${SignalMetabolitesTable}" --plotfile "${OutputPlotfile}" "${reference_dir}" "${test_dir}" > "${OutputTable}"
else
  "${IntegralsScript}" --left=${Left} --right=${Right} --where=${Where} --plotfile "${OutputPlotfile}" "${reference_dir}" "${test_dir}" > "${OutputTable}"
fi


echo "Finished" >&2

rm -rf "${reference_dir}" "${test_dir}"
exit 0
