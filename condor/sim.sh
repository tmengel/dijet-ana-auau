#!/bin/bash

INSTALLDIR="/sphenix/user/tmengel/dijet-ana-auau/install"
source /opt/sphenix/core/bin/sphenix_setup.sh -n new
source $OPT_SPHENIX/bin/setup_local.sh $INSTALLDIR

echo "Running nominal mb hijing with waveform fit"

macrodir="/sphenix/user/tmengel/dijet-ana-auau/macros/hijing"
cd $macrodir || exit 1

SEGMENT="$1"
# versionTag="08_14_2026_v001"
versionTag="$2"

NEVENTS=-1
RUN_NUM=31


OUTPUTDIR="/sphenix/tg/tg01/jets/tmengel/ppg14/sim_scaling/$versionTag"
mkdir -p $OUTPUTDIR

treeoutdir="$OUTPUTDIR/trees"
mkdir -p $treeoutdir  

dstoutdir="$OUTPUTDIR/dsts"
mkdir -p $dstoutdir

embfile="${treeoutdir}/CALO_TREE_noNoise_hijing${RUN_NUM}_pass1-$(printf "%05d" "$SEGMENT").root"
root -l -q -b "Fun4All_UEScaling_Pass1.C(${NEVENTS}, ${RUN_NUM}, ${SEGMENT}, \"${embfile}\", true)"
EXITCODE=$?
if [ $EXITCODE -ne 0 ]; then
    echo "Error: ${EXITCODE} running Fun4All_UEScaling_Pass1.C"
    exit $EXITCODE
fi

for jetid in 10 20 30 ; do 

    if [ $jetid -eq -1 ]; then

        echo "Running nominal mb hijing with waveform fit for jetid $jetid"

        this_dstoutdir="${dstoutdir}/mb_hijing/scaled"
        mkdir -p $this_dstoutdir
        # outfile="DST_SCALED_mb_hijing${RUN_NUM}_pass2-$(printf "%05d" "$SEGMENT").root"
        outfile="${this_dstoutdir}/DST_SCALED_mb_hijing${RUN_NUM}_pass2-$(printf "%05d" "$SEGMENT").root"
        root -l -q -b "Fun4All_UEScaling_Pass2.C(${NEVENTS}, ${RUN_NUM}, ${SEGMENT}, ${jetid}, \"${outfile}\" ,\"${embfile}\")"
        EXITCODE=$?
        if [ $EXITCODE -ne 0 ]; then
            echo "Error: ${EXITCODE} running Fun4All_UEScaling_Pass2.C for jetid $jetid"
        fi

        this_dstoutdir="${dstoutdir}/mb_hijing/unsclaed"
        mkdir -p $this_dstoutdir
        outfile="${this_dstoutdir}/DST_UNSCALED_mb_hijing${RUN_NUM}_pass2-$(printf "%05d" "$SEGMENT").root"
        # outfile="DST_UNSCALED_mb_hijing${RUN_NUM}_pass2-$(printf "%05d" "$SEGMENT").root"
        blankfile=""
        root -l -q -b "Fun4All_UEScaling_Pass2.C(${NEVENTS}, ${RUN_NUM}, ${SEGMENT}, ${jetid}, \"${outfile}\" ,\"${blankfile}\")"
        EXITCODE=$?
        if [ $EXITCODE -ne 0 ]; then
            echo "Error: ${EXITCODE} running Fun4All_UEScaling_Pass2.C for jetid $jetid"
        fi

    else

        echo "Running nominal mb hijing with waveform fit for jetid $jetid"
        
        this_dstoutdir="${dstoutdir}/jet${jetid}_hijing/scaled"
        mkdir -p $this_dstoutdir
        outfile="${this_dstoutdir}/DST_SCALED_jet${jetid}_hijing${RUN_NUM}_pass2-$(printf "%05d" "$SEGMENT").root"
        # outfile="DST_SCALED_jet${jetid}_hijing${RUN_NUM}_pass2-$(printf "%05d" "$SEGMENT").root"
        root -l -q -b "Fun4All_UEScaling_Pass2.C(${NEVENTS}, ${RUN_NUM}, ${SEGMENT}, ${jetid}, \"${outfile}\", \"${embfile}\")"
        EXITCODE=$?
        if [ $EXITCODE -ne 0 ]; then
            echo "Error: ${EXITCODE} running Fun4All_UEScaling_Pass2.C for jetid $jetid"
        fi
    
        blankfile=""
        this_dstoutdir="${dstoutdir}/jet${jetid}_hijing/unsclaed"
        mkdir -p $this_dstoutdir
        outfile="${this_dstoutdir}/DST_UNSCALED_jet${jetid}_hijing${RUN_NUM}_pass2-$(printf "%05d" "$SEGMENT").root"
        # outfile="DST_UNSCALED_jet${jetid}_hijing${RUN_NUM}_pass2-$(printf "%05d" "$SEGMENT").root"
        root -l -q -b "Fun4All_UEScaling_Pass2.C(${NEVENTS}, ${RUN_NUM}, ${SEGMENT}, ${jetid}, \"${outfile}\", \"${blankfile}\")"
        EXITCODE=$?
        if [ $EXITCODE -ne 0 ]; then
            echo "Error: ${EXITCODE} running Fun4All_UEScaling_Pass2.C for jetid $jetid"
        fi

    fi
done
