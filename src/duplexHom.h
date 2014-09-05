#ifndef __INC_DUPLEXHOM_H__
#define __INC_DUPLEXHOM_H__

#include "centroidalifold/folding_engine.h"
#include "vvf.h"
#include "alignment.h"

class RNAduplexHommodel
{
public:
	RNAduplexHommodel(const std::string& engine_a);
	void calculate_posterior(const TH& s1, const TH& s2, double min_aln);
	VVF get_bp();
private:
	std::string engine_a_;

	std::vector<Alignment> align_v(const TH& seq, double min_aln);
	VVF aln_duplex(Alignment a1, Alignment a2, VVF& hp);
};


#endif