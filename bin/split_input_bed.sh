#!/usr/bin/env bash

echo "===================================================================================="
echo "[nf-core/circrna]: circRNA annotation script                                        "
echo "[nf-core/circrna]: Author: Roman Reggiardo                                              "
echo "[nf-core/circrna]: Institution: Stanford University              "
echo "[nf-core/circrna]:                                                                  "
echo "[nf-core/circrna]: MIT license                                                      "
echo "===================================================================================="

mkdir -p splitBed

# $1 - input bed file
# $2 - ${task.ncpu}
# using GNU parallel, splits the bed line by line into seperate single-line files
$(cat $1 | parallel -j $2 --pipe -u --colsep '\s+' -kN1 'cat > splitBed/num{#}.bed')
