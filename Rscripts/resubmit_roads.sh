#!/bin/bash
for cid in 129 130 131 132 133 134 135 137 138 139 141 145 146 147 149 153 161 162 163 165 169 177 193 194 195 197 201 209 225; do
  cat ~/scratch/impact_assessment/Rscripts/12B_repredict_birds.sh | sbatch \
    --account=def-bayne --mem=500G \
    --export=ALL,COALITION_ID=$cid \
    --job-name=coal_${cid}
done
