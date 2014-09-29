#!/bin/sh

for i in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
do
	echo "mix_weight == ${i}"
	/Users/takashi/cbrc/ractip_hom/src/ractip_hom --mix-weight=${i} -t 0.5 AE006468_1_1.fa AE006468_CyaR_sRNA_ompX_1.fa AE006468_1_2.fa AE006468_CyaR_sRNA_ompX_2.fa > "res/mix_w/mix_w_${i}.fa"
	# 最初の_を用いて切断して、第0成分+_1_1, 2の配列に対してractip_homを行う。
	# これRubyで書きなおした方が楽かも。
done

