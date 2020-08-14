#!/bin/bash

# source conda
source /home/prashant/anaconda3/etc/profile.d/conda.sh

# activate conda
conda activate confen

# cd to build directory
cd build

rundir="$2"
if [[ ! -d $rundir ]]; then
	mkdir $rundir
fi
cd $rundir

# data and solver path
spath="/home/prashant/work/python_works/BayesForSEIRD/"
dpath="$spath""Data/"
mpath="$spath""Data/mesh/"
fwdpath="$spath""Model/seird/"
icpath="$spath""Model/ic_result_15July2020/"
submnpath="$spath""Model/subdomains/"
meshfile="mesh_5h"
samplespath="$spath""Results/prediction_11August2020/"
samplesfile="param_samples_mcmc_10.npy"
opath="../../"

# run
echo "pwd: $(pwd)"
echo "Running sim $2 out of $1+1 sims"
echo " "
python  "$fwdpath""run_prediction_qoi_11August2020.py" \
				--num_sim=$1 \
				--sim_id=$2 \
				--sim_time="80" \
				--time_step="0.1" \
				--data_path="$dpath" \
				--mesh_path="$mpath" \
				--mesh_file="mesh_5h" \
				--ic_path="$icpath" \
				--generate_ic="0" \
				--subdmn_path="$submnpath" \
				--out_path="$opath" \
				--samples_path="$samplespath" \
				--samples_file="$samplesfile"
