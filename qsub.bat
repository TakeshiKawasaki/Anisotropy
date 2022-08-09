#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -m ea
#$ -V
#
#$ -q all.q@banana
#$ -q all.q@nashi
#$ -q all.q@kuri
#$ -q all.q@lemon
#$ -q all.q@mango

./a.out
