#!/bin/bash

# source conda
# source /Users/prashant/opt/anaconda3/etc/profile.d/conda.sh
source /home/prashant/anaconda3/etc/profile.d/conda.sh

# activate conda
conda activate confen

# cd to build directory
cd mcmc_10

# set param dim
pardim="7"

# create directory if doesnt exists
pardir="param_dim_$pardim"
if [[ ! -d $pardir ]]; then
	mkdir $pardir
fi
cd $pardir

if [[ ! -d "calibration_result" ]]; then
  mkdir calibration_result
fi

rundir="$2"
if [[ ! -d $rundir ]]; then
	mkdir $rundir
fi
cd $rundir

# hyper parameters
sim_time="25.0"
time_step="0.1"
noise_inf="0.12"
noise_dec="0.05"
pcn_s="0.2"

# data and solver path
spath="../../../../../"
dpath="$spath""Data/"
mpath="$spath""Data/mesh/"
fwdpath="$spath""Model/seird/"
icpath="$spath""Model/ic_result_15July2020/"
submnpath="$spath""Model/subdomains/"
opath="../calibration_result/"
meshfile="mesh_5h"

# run
echo "pwd: $(pwd)"
echo "Running chain $2 out of $1 total chains"
echo " "
python  "$fwdpath""run_calibration_multichain_11August2020.py" \
				--nsamples=4000 \
				--chain_id=$2 \
				--sim_time="$sim_time" \
				--time_step="$time_step" \
				--data_path="$dpath" \
				--mesh_path="$mpath" \
				--mesh_file="$meshfile" \
				--ic_path="$icpath" \
				--generate_ic="0" \
				--subdmn_path="$submnpath" \
				--out_path="$opath" \
				--noise_inf="$noise_inf" \
				--noise_dec="$noise_dec" \
				--pcn_s="$pcn_s" \
