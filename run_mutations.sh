#!/bin/bash
TZ="Europe/London" date
cd /cephfs/covid/bham/climb-covid19-wrightd/COG-UK-Dashboard
source /cephfs/covid/bham/climb-covid19-wrightd/miniconda3/etc/profile.d/conda.sh
conda activate r_env
Rscript mutations.R
