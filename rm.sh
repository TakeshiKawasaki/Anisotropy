#!/bin/sh

#a=(0.1 0.2 0.3 0.4 0.5 0.6 0.8 1.0 2.0)
#a=(0.001 0.003 0.006 0.01 0.03 0.06 0.1)
#a=(0.001 0.003 0.006 0.01 0.03 0.06 0.1 0.3 0.6 1.0)
#a=(0.09 0.1 0.2 0.4 0.5)
a=(0 0.1 0.25 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 10.0 20.0 30.0 40.0 50.0 100.0)
#a=(4.0 5.0)
icc -O3 BD_eq.cpp -o BD_eq.out 

for ((i=0 ; i<19 ; i++))
#do qsub bd_eq.bat ${a[i]}
do rm ${a[i]}/*
done