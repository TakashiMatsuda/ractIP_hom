#!/bin/sh

for i in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
do
	/Users/takashi/cbrc/ractip_hom/src/ractip_hom --acc-th=i gca_bovine.fasta sample.fasta gca_bovine2.fasta sample-2.fasta
done

