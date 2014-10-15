#!/bin/sh

for i in 0.05 0.1 0.15 0.2 
do
    for j in 0.05 0.1 0.12 0.14 0.16
    do
	echo "ib-t == ${i}"
	echo "ob-t == ${j}"
	time ${RACTIP_HOM}/src/ractip_hom -t ${i} -u ${j} MicA.fa easy1.fa ompA.fa easy2.fa > newres/ibt-${i}-obt-${j}.fa
    done
done

