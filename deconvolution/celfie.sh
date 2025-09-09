#!/bin/bash

set -ue #stop on any error

source ~/miniconda3/bin/activate celfie_env

python synced_git_files/cfDNA/celfie/scripts/celfie.py celfie_trial1/atlas_uxm_02/input.txt celfie_trial1/atlas_uxm_02/output_withoutunknowns 45 -u 0 -r 30
python synced_git_files/cfDNA/celfie/scripts/celfie.py celfie_trial1/atlas_uxm_02/input.txt celfie_trial1/atlas_uxm_02/output_withunknowns 45 -u 2 -r 30
