#!/bin/bash
if [[ -z $1  ]]
then
	np=$(($(nproc)+1))
else
	np=$1
fi
if [[ -z $2  ]]
then
	run="run.py"
else
	run=$2
fi
mpirun -np $np python3 $run
