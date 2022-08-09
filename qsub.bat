#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -M takeshi.kawasaki@univ-montp2.fr
#$ -m ea
#$ -V
#
#$ -q all.q@banana
#$ -q all.q@nashi
#$ -q all.q@kuri
#$ -q all.q@lemon
#$ -q all.q@mango




./Langevin_oscillation.out