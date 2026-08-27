#!/bin/bash

INSTALLDIR="/sphenix/user/tmengel/dijet-ana-auau/install"
source /opt/sphenix/core/bin/sphenix_setup.sh -n new
source $OPT_SPHENIX/bin/setup_local.sh $INSTALLDIR

macrodir="/sphenix/user/tmengel/dijet-ana-auau/macros/truth-matching"
cd $macrodir || exit 1

INLIST="$1"
versionTag="$2"
listbase="$(basename "$INLIST" .list)"
OUTPUTDIR="$3"
mkdir -p $OUTPUTDIR
OUTPUTFILE="${OUTPUTDIR}/${listbase}_${versionTag}.root"
hadd -f -k $OUTPUTFILE $(cat $INLIST)
# root -l -q -b "match.C(\"${INLIST}\", \"${OUTPUTFILE}\")"
EXITCODE=$?
if [ $EXITCODE -ne 0 ]; then
    echo "Error: ${EXITCODE} running match.C for input list ${INLIST}"
    exit $EXITCODE
fi
