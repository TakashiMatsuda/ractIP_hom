#ifndef __INC_DUPLEXHOM_H__
#define __INC_DUPLEXHOM_H__

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include "centroidalifold/folding_engine.h"
#include "vvf.h"
#include "alignment.h"
#include "centroidalifold/util.h"

class RNAduplexHommodel
{
public:
	RNAduplexHommodel(const std::string& engine_a);
	VVF calculate_posterior(const TH& s1, const TH& s2);
private:
	std::string engine_a_;
	VVVF align_v(const TH& seq, double min_aln);
	VVF computeposterior(PROBCONS::Probcons& pc, const std::string& seq1, const std::string& seq2);
	VVVVF rnaduplex_hom(const TH& th1, const TH& th2);
	void aln_duplex(Alignment a1, Alignment a2, VVF& hp);
	void rnaduplex(const std::string& s1, const std::string& s2, VVF& hp) const;
};


#endif