#!/bin/sh
mkdir newres
for i in 0.05 0.1 0.15 0.2 
do
    for j in 0.05 0.1 0.12 0.14 0.16
    do
	echo "ib-t == ${i}"
	echo "ob-t == ${j}"
	time ${RACTIPHOM}/src/ractip_hom -t ${i} -u ${j} OxyS.fa U00096_OxyS_sRNA_fhlA_1.fa fhlA.fa U00096_OxyS_sRNA_fhlA_2.fa > newres/ibt-${i}-obt-${j}.fa
    done
done

