#!/bin/bash
source /opt/sphenix/core/bin/sphenix_setup.sh -n new
export INSTALLDIR="/sphenix/user/tmengel/dijet-ana-auau/install"
mkdir -p "$INSTALLDIR"
source $OPT_SPHENIX/bin/setup_local.sh $INSTALLDIR
