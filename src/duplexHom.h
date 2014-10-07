#ifndef __INC_DUPLEXHOM_H__
#define __INC_DUPLEXHOM_H__

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include <unistd.h>
#include <cstdlib>
#include <sys/time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <boost/algorithm/string.hpp>


#include "centroidalifold/folding_engine.h"
#include "centroidalifold/probconsRNA/probcons.h"
#include "vvf.h"
#include "alignment.h"
#include "centroidalifold/util.h"

class RNAduplexHommodel
{
public:
	~RNAduplexHommodel();
	RNAduplexHommodel(const std::string& engine_a="probcons");
	VVF calculate_posterior(const TH& s1, const TH& s2, double wh);

	// These functions are temporally publiciated for the test.
	VVF computeposterior(PROBCONS::Probcons& pc, const std::string& seq1, const std::string& seq2);
private:
	std::string engine_a_;
	VVVF align_v(const TH& seq, double min_aln);

	VVVVF rnaduplex_hom(const TH& th1, const TH& th2);
	void rnaduplex(const std::string& s1, const std::string& s2, VVF& hp) const;
};


#endif