#! /bin/bash

# integrated test process.
# すべてここに書く。このスクリプトを動かすことによってテストが行われる。

kdir res
mkdir res/bovine
for i in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
do
	/Users/takashi/cbrc/ractip_hom/src/ractip_hom -a $i 
	echo "alpha_${i}\n" >> res/bovine/alpha_total.fa
	/Users/takashi/cbrc/ractip_hom/src/ractip_hom -a $i gca_bovine.fasta sample_nongap.fa gca_bovine2.fasta sample2_nongap.fa >> res/bovine/alpha_total.fa
done
/Users/takashi/cbrc/ractip_hom/src/ractip_hom 