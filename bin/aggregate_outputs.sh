#!/usr/bin/env bash

echo "===================================================================================="
echo "[nf-core/circrna]: circRNA annotation script                                        "
echo "[nf-core/circrna]: Author: Barry Digby                                              "
echo "[nf-core/circrna]: Institution: National University of Ireland, Galway              "
echo "[nf-core/circrna]:                                                                  "
echo "[nf-core/circrna]: MIT license                                                      "
echo "===================================================================================="


# remove the trailing commas that appear on end of Exon Size, Exon Starts (col 11, 12)
cat bed12/*.bed12.bed > master_bed12.bed.tmp

awk 'BEGIN{FS=OFS="\t"} {gsub(/,$/,"",$11);gsub(/,$/,"",$12)} 1' master_bed12.bed.tmp > master_bed12.bed && rm master_bed12.bed.tmp

# remove the pre-processed bed lines
rm -rf splitBed/

echo "                   Thank you for using nf-core/circrna! - Barry                     "
echo "===================================================================================="
echo "===================================================================================="
