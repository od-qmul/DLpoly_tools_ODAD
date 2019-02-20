#!/bin/bash

MINTEMP=300
MAXTEMP=3000
COOLRATE=100 #K/ps


for ((t=$MAXTEMP;t>=$MINTEMP;t=t-$COOLRATE)) do
        DIR=${t}
        cd $DIR
        python3 ~/bin/Sum_STATIS/read_statis.py 1 30 31 32 33 34 35
        cd ..
done

