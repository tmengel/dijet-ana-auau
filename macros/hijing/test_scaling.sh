#!/bin/bash

INSTALLDIR="/sphenix/user/tmengel/hijing_scaling/HIJING_SCALING/install"
source /opt/sphenix/core/bin/sphenix_setup.sh -n new
source $OPT_SPHENIX/bin/setup_local.sh $INSTALLDIR

echo "Running nominal mb hijing with waveform fit"
embfile="CALO_TREE_noNoise_hijing31_pass1-00000.root"
root -l -q -b "Fun4All_UEScaling_Pass1.C(10, 31, 0, \"${embfile}\", true)"
echo $?

jetid=10
echo "Running nominal mb hijing with waveform fit for jetid $jetid"
outfile="DST_SCALED_jet${jetid}_hijing31_pass2-00000.root"
root -l -q -b "Fun4All_UEScaling_Pass2.C(10, 31, 0, ${jetid}, \"${embfile}\", \"${outfile}\")"
echo $?

echo "Checking output DST"
dstoutfile="CALO_TREE_DST_CHECK_sHijing_0_20fm-00000031-00000.root"
root -l -q -b "Fun4All_UEScaling_CheckDst.C(10, \"${outfile}\", \"${dstoutfile}\")"
echo $?

echo "Checking output DST nominal"
dstoutfile_nominal="CALO_TREE_DST_CHECK_nominal_emb_0_20fm-00000031-00000.root"
root -l -q -b "Fun4All_UEScaling_CheckDst.C(10, \"${outfile}\", \"${dstoutfile_nominal}\", true)"
echo $?

echo "Making plots"
outdir="scaled_plots"
mkdir -p $outdir
root -l -q -b "plot_tower_energy_comp.C(\"${dstoutfile}\", \"${dstoutfile_nominal}\", \"${outdir}\")"
echo $?
