#! /bin/sh

for i in 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1
do
	echo "hyb_mix_weight == ${i}"
	mkdir bp-matrix-wh-${i}
	$RACTIP_HOM/src/ractip_hom --hyb-mix-w=${i} --output-dir=bp-matrix-wh-${i} AE006468_1_1.fa AE006468_CyaR_sRNA_ompX_1.fa AE006468_1_2.fa AE006468_CyaR_sRNA_ompX_2.fa
done
