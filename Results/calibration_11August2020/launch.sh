#!/bin/bash

# mcmc id (each id has different combination of noise and pcn s parameter)
mcmc_id="$1"
mcmc_dir="mcmc_$mcmc_id"

# create directories
if [[ ! -d "$mcmc_dir" ]]; then
	mkdir $mcmc_dir
fi

num_chains="4"
for i in $(seq 1 $num_chains); do screen -dmS mcmc-$mcmc_id-$i "./mcmc_$mcmc_id.sh" $num_chains $i; done

# for i in $(seq 1 $num_chains); do "./mcmc_$mcmc_id.sh" $num_chains $i; done
#"./mcmc_$mcmc_id.sh" 1 1 

