#!/bin/zsh
#
#$ -cwd
#$ -o test.o
#$ -e test.e
#$ -V
#$ -l h_cpu=12:00:00
#$ -l h_rss=2G
#$ -P pitz
##$ -pe multicore 32
##$ -R y

cd .
gencore test.in 2>&1 | tee test.log
