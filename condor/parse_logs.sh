#!/bin/bash

LOGDIR="/sphenix/user/tmengel/dijet-ana-auau/condor/logs/08_14_2026_v001"

echo "Scanning ${LOGDIR}/*.out for ERROR..."
echo

FOUND=0

for LOGFILE in "$LOGDIR"/*.out; do
    [[ -e "$LOGFILE" ]] || continue

    if grep -q "ERROR" "$LOGFILE"; then
        FOUND=1

        BASENAME=$(basename "$LOGFILE")

        echo "========================================"
        echo "ERRORS in: $BASENAME"
        echo "========================================"

        grep -n "ERROR" "$LOGFILE"

        echo
    fi
done

if [ "$FOUND" -eq 0 ]; then
    echo "No ERROR messages found."
fi