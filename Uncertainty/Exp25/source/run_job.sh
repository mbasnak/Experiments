#!/bin/sh
exp=$1
file="/n/scratch3/users/m/mb491/data/Exp25/$exp/2p/sidtid.txt"
while IFS=',' read -r n1 n2
do
	sbatch run_continuous_PB_analysis.sh $exp $n1 $n2
	sleep 1
done<"$file"