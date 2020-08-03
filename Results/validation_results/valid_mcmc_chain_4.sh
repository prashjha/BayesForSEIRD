#!/bin/bash

# source conda
source /home/prashant/anaconda3/etc/profile.d/conda.sh

# activate conda
conda activate confen

# cd to build directory
cd mcmc_4

# set param dim
pardim="7"

# create directory if doesnt exists
pardir="param_dim_$pardim"
if [[ ! -d $pardir ]]; then
	mkdir $pardir
fi
cd $pardir

if [[ ! -d "validation_result" ]]; then
  mkdir validation_result
fi

rundir="$2"
if [[ ! -d $rundir ]]; then
	mkdir $rundir
fi
cd $rundir

# hyper parameters
sim_time="30.0"
time_step="0.1"
noise_inf="0.08"
noise_dec="0.04"
pcn_s="0.3"

# data and solver path
spath="../../../../"
dpath="$spath""Data/"
mpath="$spath""Data/mesh/"
fwdpath="$spath""Model/seird/"
icpath="$spath""Model/ic_result/"
submnpath="$spath""Model/subdomains/"
meshfile="mesh_5h"
priorpath="$spath""Results/calibration_results/"
priorfile="calibration_posterior_samples.npy"
opath="../validation_result/"

# run
echo "pwd: $(pwd)"
echo "Running chain $2 out of $1 total chains"
echo " "
python  "$fwdpath""run_validation_multichain.py" \
				--nsamples=5000 \
				--chain_id=$2 \
				--sim_time="$sim_time" \
				--time_step="$time_step" \
				--data_path="$dpath" \
				--mesh_path="$mpath" \
				--mesh_file="$meshfile" \
				--ic_path="$icpath" \
				--generate_ic="0" \
				--subdmn_path="$submnpath" \
				--samples_path="$priorpath" \
				--samples_file="$priorfile" \
				--out_path="$opath" \
				--noise_inf="$noise_inf" \
				--noise_dec="$noise_dec" \
				--pcn_s="$pcn_s" \
