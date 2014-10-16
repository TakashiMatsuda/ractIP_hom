#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -v PATH
#$ -v LD_LIBRARY_PATH
#$ -e ./stderr.txt
#$ -o ./stdout.txt
#$ -N RactIP_homologous
#$ -M mazda.takasky@gmail.com

PARALLELLIMIT=64
NUM_PROCESS=0
mkdir tritest_hyb_t_u

for i in 0.05 0.1 0.15 0.2 0.3 0.35 0.4
do
    for j in 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4
    do
	for k in 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1
	do
	    # PARALLEL THREADS LIMITATION
	    if [ ${NUM_PROCESS} -gt ${PARALLELLIMIT} ] ; then
		wait
		NUM_PROCESS=0
	    fi
	    echo "ib-t == ${i}"
	    echo "ob-t == ${j}"
	    ${RACTIPHOM}/src/ractip_hom -t ${i} -u ${j} --hyb-mix-w=${k} RyhB.fa U00096_RyhB_sRNA_SodB_1.fa SodB.fa U00096_RyhB_sRNA_SodB_2.fa > tritest_hyb_t_u/ibt-${i}-obt-${j}-hyb-${k}.fa &
	    NUM_PROCESS=`expr ${NUM_PROCESS} + 1`
	done
    done
done
