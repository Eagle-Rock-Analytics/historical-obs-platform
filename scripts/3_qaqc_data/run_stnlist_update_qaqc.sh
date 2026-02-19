#!/bin/bash
# Run stnlist_update_qaqc.py for all networks.
# Usage: ./run_stnlist_update_qaqc.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
START_TIME=$SECONDS

NETWORKS=(
    ASOSAWOS CAHYDRO CIMIS CW3E CDEC CNRFC CRN CWOP HADS
    HNXWFO HOLFUY HPWREN LOXWFO MAP MTRWFO NCAWOS NOS-NWLON
    NOS-PORTS otherisd RAWS SGXWFO SHASAVAL VCAPCD MARITIME
    NDBC SCAN SNOTEL VCAPCD
)

TOTAL=${#NETWORKS[@]}
COUNT=0

for network in "${NETWORKS[@]}"; do
    COUNT=$((COUNT + 1))
    echo "[$COUNT/$TOTAL] Updating QAQC station list for $network..."
    python "$SCRIPT_DIR/stnlist_update_qaqc.py" "$network"
done

ELAPSED=$(( SECONDS - START_TIME ))
echo "Done. Total runtime: $((ELAPSED / 60))m $((ELAPSED % 60))s"
