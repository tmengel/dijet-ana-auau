#!/bin/bash

INSTALLDIR="/sphenix/user/tmengel/hijing_scaling/HIJING_SCALING/install"
source /opt/sphenix/core/bin/sphenix_setup.sh -n new
source $OPT_SPHENIX/bin/setup_local.sh $INSTALLDIR


# echo "Running nominal mb hijing"
# outfile="CALO_TREE_nominalMB_hijing31_pass1-00000.root"
# root -l -q -b "Fun4All_UEScaling_Pass1.C(-1, 31, 0, \"${outfile}\", false)"
# echo $?

# echo "Running nominal mb hijing with waveform fit"
# outfile="CALO_TREE_noNoise_hijing31_pass1-00000.root"
# root -l -q -b "Fun4All_UEScaling_Pass1.C(-1, 31, 0, \"${outfile}\", true)"
# echo $?

echo "Making plots"
outdir="plots"
mkdir -p $outdir
root -l -q -b "plot_tower_energy_comp.C(\"CALO_TREE_nominalMB_hijing31_pass1-00000.root\", \"CALO_TREE_noNoise_hijing31_pass1-00000.root\", \"${outdir}\")"
echo $?

