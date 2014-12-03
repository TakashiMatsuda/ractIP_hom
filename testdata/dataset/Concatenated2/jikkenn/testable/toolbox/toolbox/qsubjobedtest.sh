#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -v PATH
#$ -v LD_LIBRARY_PATH
#$ -e ./stderr.txt
#$ -o ./stdout.txt
#$ -N RactIP_homologous
#$ -M mazda.takasky@gmail.com


${RACTIPHOM}/src/ractip_hom -a $1 GcvB.fa ok_AE006468_GcvB_sRNA_livK_1.fa livK.fa ok_AE006468_GcvB_sRNA_livK_2.fa > tritest_hyb_t_u/alpha-$1.fa
