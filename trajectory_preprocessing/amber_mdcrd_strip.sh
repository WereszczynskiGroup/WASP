#!/bin/bash

prmtop=$1
incrd=$2
nbp=$3
startframe=$4
endframe=$5
stride=$6

cpptraj<<EOF
parm $prmtop
trajin $incrd $startframe $endframe $stride
strip !(@C1')
trajout stripped_C1.mdcrd nobox
run
EOF

cpptraj<<EOF
parm $prmtop
parmstrip !(@C1')
parmwrite out stripped_C1.prmtop nobox
EOF
