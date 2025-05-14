#!/usr/bin/bash

# Ensure bash v >= 4. is used.
if (( BASH_VERSINFO[0] < 4 )); then
    echo "${__name__} requires bash v4 or greater"
    echo "Current bash version: ${BASH_VERSION}"
    exit 1
fi

# Get the root path of the project
ROOT=$(cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)/.." &> /dev/null && pwd)

TESTDIR=${ROOT}/tests/dir.template/atom_samples

# Create the testdir if it does not exists yet
if [ -d ${TESDIR} ]; then
    rm -rf ${TESTDIR}
fi
mkdir -p ${TESTDIR}

# Check if atomsk is available in the path
if ! [ -x "$(command -v atomsk)" ]; then
    echo "ATOMSK was not found in the PATH."
    exit 1
fi

# create samples
cd ${TESTDIR}

a=3.58
s=Cu
sgrp=fcc
rx=2  # set rx != ry != rz avoid cubic supercell
ry=3
rz=4

for tilt in 0 1; do
    if (( $tilt == 1 )); then
        tilt_suffix="_triclinic"
        skdef="-def yz 1%"
    else
        tilt_suffix=""
        skdef=""
    fi

    outname=${s}_${sgrp}_${rx}x${ry}x${rz}${tilt_suffix}.lmp
    atomsk_cmd="--create ${sgrp} ${a} ${s} -duplicate ${rx} ${ry} ${rz} ${skdef} -overwrite ${outname}"

    atomsk ${atomsk_cmd}
    gzip -fvk ${outname}
    bzip2 -fvk ${outname}
    xz -fvk ${outname}

done
