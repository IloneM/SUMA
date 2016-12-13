#!/bin/bash
if [[ -z $1  ]]
then
	np=$(($(nproc)+1))
else
	np=$1
fi
mpirun -np $np python3 run.py
