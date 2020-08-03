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
spath="../../../../"
dpath="$spath""Data/"
mpath="$spath""Data/mesh/"
fwdpath="$spath""Model/seird/"
icpath="$spath""Model/ic_result/"
submnpath="$spath""Model/subdomains/"
meshfile="mesh_5h"
samplespath="$spath""Results/validation_results/"
samplesfile="validation_posterior_samples.npy"
opath="./"

# run
echo "pwd: $(pwd)"
echo "Running sim $2 out of $1+1 sims"
echo " "
python  "$fwdpath""run_prediction_qoi.py" \
				--num_sim=$1 \
				--sim_id=$2 \
				--sim_time="110" \
				--time_step="0.1" \
				--fwd_path="$fwdpath" \
				--data_path="$dpath" \
				--mesh_path="$mpath" \
				--mesh_file="mesh_5h" \
				--ic_path="$icpath" \
				--generate_ic="0" \
				--subdmn_path="$submnpath" \
				--out_path="$opath" \
				--samples_path="$samplespath" \
				--samples_file="$samplesfile"
