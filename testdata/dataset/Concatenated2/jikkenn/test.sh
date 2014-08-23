#!/bin/sh

# つくった32このテストデータに対して、実際にractip_homを走らせてみる。まずはパラメータチューニング、次にractipとの性能比較。
# スコアはppvとsensitivityだけでよい。
# 計算時間を別に測っておいてもよい。が、それは必須ではないので後回し。
mkdir res
mkdir res/bovine
for i in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
do
	/Users/takashi/cbrc/ractip_hom/src/ractip_hom -a $i gca_bovine.fasta sample_nongap.fa gca_bovine2.fasta sample2_nongap.fa > res/bovine/alpha_${i}.fa
	echo "alpha_${i}\n" >> res/bovine/alpha_total.fa
	/Users/takashi/cbrc/ractip_hom/src/ractip_hom -a $i gca_bovine.fasta sample_nongap.fa gca_bovine2.fasta sample2_nongap.fa >> res/bovine/alpha_total.fa
done

	