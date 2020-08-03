#!/bin/bash

# create directories
if [[ ! -d "qoi" ]]; then
	mkdir qoi
fi

if [[ ! -d "log" ]]; then
	mkdir log
fi

if [[ ! -d "build" ]]; then
	mkdir build
fi

num_sim="19"
for i in $(seq 0 $num_sim); do screen -dmS predqoi-cal-$i ./run.sh $num_sim $i; done

#for i in $(seq 0 $num_sim); do ./run.sh $num_sim $i; done
