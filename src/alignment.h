#ifndef __INC_ALIGNMENT_H__
#define __INC_ALIGNMENT_H__

#include "centroidalifold/folding_engine.h"
#include "vvf.h"

class Alignment : public std::pair<double, std::vector<std::string> >
{
public:
	Alignment(std::vector<std::string> seq)
	{
		this->second = seq;
	};
	double get_prob()
	{
		return this->first;
	};
	
	std::vector<std::string> get_seq()
	{
		return this->second;
	};
private:
};

#endif