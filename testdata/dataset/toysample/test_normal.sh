#!/bin/sh

mkdir res_normal
mkdir res_normal/bovine
for i in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
do
	ractip -a $i gca_bovine.fasta sample_nongap.fa gca_bovine2.fasta sample2_nongap.fa > res_normal/bovine/alpha_${i}.fa
	echo "alpha_${i}" >> res_normal/bovine/alpha_total.fa
	ractip -a $i gca_bovine.fasta sample_nongap.fa gca_bovine2.fasta sample2_nongap.fa >> res_normal/bovine/alpha_total.fa
done