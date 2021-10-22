#!/bin/zsh
#
#$ -cwd
#$ -o test_farm.o
#$ -e test_farm.e
#$ -V
#$ -l h_cpu=12:00:00
#$ -l h_rss=2G
#$ -P pitz
#$ -pe multicore 16
#$ -R y

python injector_optimization_demo.py