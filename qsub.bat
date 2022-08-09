#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -M kawasaki@r.phys.nagoya-u.ac.jp
#$ -m ea
#$ -V
#
#$ -q all.q@banana

./a.out
