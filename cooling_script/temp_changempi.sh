#!/bin/bash
#$ -cwd           # Set the working directory for the job to the current directory
#$ -pe smp 4     # Request 1 core
#$ -l h_rt=24:0:0 # Request 24 hour runtime
#$ -l h_vmem=2G   # Request 1GB RAM
#$ -N BASE10cool

module load dl_poly/4.08-openmpi

MINTEMP=300
MAXTEMP=3000
COOLRATE=10 #K/ps

DIR=${MAXTEMP}
mkdir $DIR
cp init/* $DIR
cd $DIR
	newtemp="temperature            "$MAXTEMP
	sed -i "/temperature/c $newtemp" ./CONTROL
mpirun -np 4 DLPOLY.Z
cd ..

for ((t=$MAXTEMP-$COOLRATE;t>=$MINTEMP;t=t-$COOLRATE)) do
        echo "Running at temperature: " $t
	oldDIR=$DIR
        DIR=${t}
        mkdir $DIR
        cp $oldDIR/CONTROL $DIR
	cp $oldDIR/REVCON $DIR/CONFIG
	cp $oldDIR/FIELD $DIR
        cd $DIR
	newtemp="temperature            "$t
	sed -i "/temperature/c $newtemp" ./CONTROL
	mpirun -np 4 DLPOLY.Z

        cd ..
done
