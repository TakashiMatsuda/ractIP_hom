  /*
 * $Id$
 * 
 * Copyright (C) 2010 Kengo Sato
 *
 * This file is part of RactIP.
 *
 * RactIP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * RactIP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with RactIP.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "cmdline.h"
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
#include "fa.h"
#include "ip.h"
#include "duplexHom.h"
#include "centroidalifold/aln.h"
#include "centroidalifold/mea.h"
#include "centroidalifold/folding_engine.h"
#include "centroidalifold/engine/contrafold.h"
#include "centroidalifold/engine/contrafoldm.h"
#include "centroidalifold/engine/mccaskill.h"
#include "centroidalifold/engine/alifold.h"
#include "centroidalifold/engine/pfold.h"
#include "centroidalifold/engine/averaged.h"
#include "centroidalifold/engine/mixture.h"
#include "centroidalifold/engine/aux.h"
#include "centroidalifold/engine/contrafoldhom.h"
#include "centroidalifold/engine/mccaskillhom.h"

#define BP_OUTPUT 0


namespace Vienna {
extern "C" {
#include <ViennaRNA/fold.h>
#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/part_func_up.h>
#include <ViennaRNA/part_func_co.h>
#include <ViennaRNA/utils.h>
#ifndef HAVE_VIENNA20
#include <ViennaRNA/PS_dot.h>
  extern int eos_debug;
#endif
#include "pf_duplex.h"
  extern void read_parameter_file(const char fname[]);
};
};

namespace uShuffle {
extern "C" {
#include "ushuffle.h"
};
};

extern "C" {
#include "boltzmann_param.h"
};


#include "vvf.h"

typedef unsigned int uint;
#ifdef HAVE_VIENNA18
typedef Vienna::plist pair_info;
#else
typedef Vienna::pair_info pair_info;
#endif

class RactIP
{
public:
  RactIP()
    : alpha_(0.5),
      beta_(0.0),
      th_hy_(0.2),
      th_ss_(0.0),
      th_ac_(0.0),
      WH(0.6),
      max_w_(0),
      min_w_(0),
      enable_zscore_(0),
      num_shuffling_(1000),
      seed_(0),
      in_pk_(true),
      use_contrafold_(true),
      use_pf_duplex_(true),
      stacking_constraints_(true),
      show_energy_(false),
      allow_concat_(true),
      run_with_modena_(false),
      n_th_(1),
      rip_file_(),
      param_file_(),
      fa1_(),
      fa2_(),
      mix_w(),
      engine(),
      engine_a(),
      bp_dirname()

  {
  }
  
  ~RactIP()
  {
  }

  RactIP& parse_options(int& argc, char**& argv);
  int run();
  float solve(TH& s1, TH& s2, std::string& r1, std::string& r2, FoldingEngine<TH>* cf1, FoldingEngine<TH>* cf2);

  static void calculate_energy(const std::string s1, const std::string& s2,
                               const std::string r1, const std::string& r2,
                               float& e1, float& e2, float& e3);
  template <typename T> typename std::vector<T> change_l_v(std::list<T> list_seq) const;

private:
    void transBP_centroidfold_ractip(BPTable bp_centroidfold, VF& bp, VI& offset, const std::string& seq) const;
  int ComputeRowOffset(int i, int N, int w /*=0*/) const;
  void homfold(const TH& seq, VF& bp, VI& offset, VVF& up, FoldingEngine<TH>* cf) const;
  //void contrafold(const std::string& seq, VF& bp, VI& offset, VVF& up) const;
  void rnaduplex_aln(const Aln& a1, const Aln& a2, VVF& hp) const;
  void rnaduplex_hom(TH& a1, TH& a2, VVF& hp, double wh) const;
  //void contraduplex(const std::string& seq1, const std::string& seq2, VVF& hp) const;
  void rnafold(const std::string& seq, VF& bp, VI& offset) const;
  void rnafold(const std::string& seq, VF& bp, VI& offset, VVF& up, uint max_w) const;
  void rnaduplex(const std::string& seq1, const std::string& seq2, VVF& hp) const;
  void load_from_rip(const char* filename,
                     const std::string& s1, const std::string& s2,
                     VF& bp1, VI& offset1, VF& bp2, VI& offset2, VVF& hp) const;

private:
  // options
  float alpha_;                // weight for the hybridization score
  float beta_;                 // weight for unpaired bases
  float th_hy_;                // threshold for the hybridization probability
  float th_ss_;                // threshold for the base-pairing probability
  float th_ac_;                // threshold for the accessible probability
  double WH; 
  int max_w_;                  // maximum length of accessible regions
  int min_w_;                  // mimimum length of accessible regions
  int enable_zscore_;          // flag for calculating z-score
  int num_shuffling_;          // the number of shuffling for calculating z-score
  uint seed_;                  // seed for random()
  bool in_pk_;                 // allow internal pseudoknots or not
  bool use_contrafold_;        // use CONTRAfold model or not
  bool use_pf_duplex_;
  bool stacking_constraints_;
  bool show_energy_;
  bool allow_concat_;
  bool run_with_modena_;
  int n_th_;                   // the number of threads
  std::string rip_file_;
  std::string param_file_;
  std::string fa1_;
  std::string fa2_;
  std::string aln1_;
  std::string aln2_;
  std::string bp_dirname;
  

  /**
   * member for Alignment interaction
   **/
  std::vector<float> mix_w;   // mixture weights of inference engines
  std::vector<std::string> engine;
  std::vector<std::string> engine_a;


};


template 
<typename T>
std::vector<T>
RactIP::
change_l_v(std::list<T> list_seq) const{
  std::vector<T> vt_seq;
  typename std::list<T>::iterator itr_list;
  for(itr_list=list_seq.begin(); itr_list!=list_seq.end(); itr_list++){
    vt_seq.push_back(*itr_list);
  }
  return vt_seq;
}

void
RactIP::
transBP_centroidfold_ractip(BPTable bp_centroidfold, VF& bp, VI& offset, const std::string& seq) const
{
  int bpsize=bp_centroidfold.size();
  uint L = seq.size();
  //std::cout << "bpsize: " << std::endl;
  bp.clear();
  bp.resize((L+1)*(L+2)/2);
  offset.resize(L+1);
  for(uint i=0; i<=L; ++i)
    {
      offset[i]=i*((L+1)+(L+1)-i-1)/2;
    }
  for (int i=0; i<bpsize; i++)
    {
      int j=0;
      for (j=i; j<bpsize; j++)
	{
	  //bp.push_back(bp_centroidfold(i,j));
	  bp[offset[i]+j]=bp_centroidfold(i,j);
	}
    }
  
  //L = sstruct.G;
  //int max_bp_dist=0;
  //InferenceEngine<float> en(false);
  
  //  for (int i = 0; i <= bpsize; i++)
  //    {
  //      offset.push_back(ComputeRowOffset(i,bpsize+1,max_bp_dist));
      //allow_unpaired_position[i] = 1;
      //loss_unpaired_position[i] = RealT(0);
  //    }
}


int 
RactIP::ComputeRowOffset(int i, int N, int w /*=0*/) const
{
  //Assert(i >= 0 && i <= N, "Index out-of-bounds.");大文字Assertは見当たらない
  if (!(i>=0 && i<=N))
    {
      std::cout << "Index out-of-bounds." << std::endl;
    }

#define USE_EFFICIENT_WINDOW
#ifdef USE_EFFICIENT_WINDOW
    if (w==0)
    {
	// equivalent to:
	//   return N*(N+1)/2 - (N-i)*(N-i+1)/2 - i;
	return i*(N+N-i-1)/2;
    }
    else
    {
	return i*w - i;
    }
#else
    // equivalent to:
    //   return N*(N+1)/2 - (N-i)*(N-i+1)/2 - i;
    return i*(N+N-i-1)/2;
#endif
}


void
RactIP::
homfold(const TH& seq, VF& bp, VI& offset, VVF& up, FoldingEngine<TH>* cf) const
{
  // 塩基対確率を計算して取り出す関数を作っておいて、呼び出す。
  // basepair probabilityの計算
  //cf->stochastic_fold(seq, num_samples, *out)
  cf->calculate_posterior(seq);// error, ギャップ'-'を想定して改修したい
  BPTable bp_centroidfold;
  bp_centroidfold=cf->get_bp();

  // bpの変換
  transBP_centroidfold_ractip(bp_centroidfold, bp, offset, seq.first);

  const uint L=seq.size();
  up.resize(L, VF(1, 1.0));
  for (uint i=0; i!=L; ++i)
  {
    for (uint j=0; j<i; ++j)
      up[i][0] -= bp[offset[j+1]+(i+1)];
    for (uint j=i+1; j<L; ++j)
      up[i][0] -= bp[offset[i+1]+(j+1)];
    up[i][0] = std::max(0.0f, up[i][0]);
  }
}

/**
void
RactIP::
contrafold(const std::string& seq, VF& bp, VI& offset, VVF& up) const
{
  SStruct ss("unknown", seq);
  ParameterManager<float> pm;
  InferenceEngine<float> en(false);
  VF w = GetDefaultComplementaryValues<float>();
  bp.resize((seq.size()+1)*(seq.size()+2)/2);
  std::fill(bp.begin(), bp.end(), 0.0);
  en.RegisterParameters(pm);
  en.LoadValues(w);
  en.LoadSequence(ss);
  en.ComputeInside();
  en.ComputeOutside();
  en.ComputePosterior();
  en.GetPosterior(0, bp, offset);
  
  const uint L=seq.size();
  up.resize(L, VF(1, 1.0));
  for (uint i=0; i!=L; ++i)
  {
    for (uint j=0; j<i; ++j)
      up[i][0] -= bp[offset[j+1]+(i+1)];
    for (uint j=i+1; j<L; ++j)
      up[i][0] -= bp[offset[i+1]+(j+1)];
    up[i][0] = std::max(0.0f, up[i][0]);
  }
}
**/
/**
void
RactIP::
contraduplex(const std::string& seq1, const std::string& seq2, VVF& hp) const
{
  hp.resize(seq1.size()+1, VF(seq2.size()+1));
  SStruct ss1("unknown", seq1), ss2("unknown", seq2);
  ParameterManager<float> pm;
  DuplexEngine<float> en(false);
  VF w = GetDefaultComplementaryValues<float>();
  VF ip;
  en.RegisterParameters(pm);
  en.LoadValues(w);
  en.LoadSequence(ss1, ss2);
  en.ComputeInside();
  en.ComputeOutside();
  en.ComputePosterior();
  en.GetPosterior(th_hy_, ip);
  for (uint i=0; i!=seq1.size(); ++i)
    for (uint j=0; j!=seq2.size(); ++j)
      hp[i+1][j+1] = ip[en.GetOffset(i+1)+(j+1)];
}
**/

void
RactIP::
rnafold(const std::string& seq, VF& bp, VI& offset) const
{
  uint L=seq.size();
  bp.resize((L+1)*(L+2)/2);
  offset.resize(L+1);
  for (uint i=0; i<=L; ++i)
    offset[i] = i*((L+1)+(L+1)-i-1)/2;
#if 0
  std::string str(seq.size()+1, '.');
  float min_en = Vienna::fold(const_cast<char*>(seq.c_str()), &str[0]);
  float sfact = 1.07;
  float kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
  Vienna::pf_scale = exp(-(sfact*min_en)/kT/seq.size());
#else
  Vienna::pf_scale = -1;
#endif
#ifndef HAVE_VIENNA20
  Vienna::init_pf_fold(L);
#endif
  Vienna::pf_fold(const_cast<char*>(seq.c_str()), NULL);
#ifdef HAVE_VIENNA20
  FLT_OR_DBL* pr = Vienna::export_bppm();
  int* iindx = Vienna::get_iindx(seq.size());
#else
  FLT_OR_DBL* pr = Vienna::pr;
  int* iindx = Vienna::iindx;
#endif
  for (uint i=0; i!=L-1; ++i)
    for (uint j=i+1; j!=L; ++j)
      bp[offset[i+1]+(j+1)] = pr[iindx[i+1]-(j+1)];
  Vienna::free_pf_arrays();
}

void
RactIP::
rnafold(const std::string& seq, VF& bp, VI& offset, VVF& up, uint max_w) const
{
  uint L=seq.size();
  bp.resize((L+1)*(L+2)/2);
  offset.resize(L+1);
  for (uint i=0; i<=L; ++i)
    offset[i] = i*((L+1)+(L+1)-i-1)/2;
#if 0
  std::string str(seq.size()+1, '.');
  float min_en = Vienna::fold(const_cast<char*>(seq.c_str()), &str[0]);
  float sfact = 1.07;
  float kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
  Vienna::pf_scale = exp(-(sfact*min_en)/kT/seq.size());
#else
  Vienna::pf_scale = -1;
#endif
#ifndef HAVE_VIENNA20
  Vienna::init_pf_fold(L);
#endif
  Vienna::pf_fold(const_cast<char*>(seq.c_str()), NULL);
#ifdef HAVE_VIENNA20
  FLT_OR_DBL* pr = Vienna::export_bppm();
  int* iindx = Vienna::get_iindx(seq.size());
#else
  FLT_OR_DBL* pr = Vienna::pr;
  int* iindx = Vienna::iindx;
#endif
  for (uint i=0; i!=L-1; ++i)
    for (uint j=i+1; j!=L; ++j)
      bp[offset[i+1]+(j+1)] = pr[iindx[i+1]-(j+1)];

  up.resize(L, VF(max_w));
  Vienna::pu_contrib* pu = Vienna::pf_unstru(const_cast<char*>(seq.c_str()), max_w);
  assert((int)L==pu->length);
  for (uint i=0; i!=L; ++i)
    for (uint j=0; j!=max_w; ++j)
      up[i][j]=pu->H[i+1][j] + pu->I[i+1][j] + pu->M[i+1][j] + pu->E[i+1][j];
#ifdef HAVE_VIENNA20
  Vienna::free_pu_contrib_struct(pu);
#else
  Vienna::free_pu_contrib(pu);
#endif
  Vienna::free_pf_arrays();
}


// 
void
RactIP::
rnaduplex_hom(TH& s1, TH& s2, VVF& hp, double wh) const
{
  RNAduplexHommodel rdh;
  hp = rdh.calculate_posterior(s1, s2, wh);
}

void
RactIP::
rnaduplex_aln(const Aln& a1, const Aln& a2, VVF& hp) const
{
  // アラインメントに含まれるタグが同じことを前提とする。
  // 進化で構造が保存されるという考え方を利用するには、この条項が必要。
  // forループを書き換える。

  // hpに平均値を入れていく。
  // hybridization probability
  std::list<std::string> s1=a1.seq();
  std::list<std::string> s2=a2.seq();

  uint size1=a1.num_aln();
  uint size2=a2.num_aln();
  std::list<std::string>::iterator it1=s1.begin();
  std::list<std::string>::iterator it2=s2.begin();
  VVVVF vhp(size1);
  int i=0;
  ///////////
  it1++;
  //it1++;
  it2++;

  if (size1 == size2)
  {
    it1 = s1.begin();
    it2 = s2.begin();
    VVVF vtr_hp(size1);
    for (int i = 0; i < size1; i++)
    {
      rnaduplex(*it1, *it2, vtr_hp[i]);
      it1++;
      it2++;
    }
    VVVF::iterator itr_vtr_hp;
    int L = vtr_hp[0].size();
    int M = vtr_hp[0][0].size();
    hp.resize(L+1);// debug temporally
    for (int i=0; i<L; i++)
    {
      for (int j = 0; j < M; j++)
      {
        double sum = 0;
        for (itr_vtr_hp = vtr_hp.begin(); itr_vtr_hp != vtr_hp.end(); itr_vtr_hp++)
        {
          sum += (*itr_vtr_hp)[i][j];
        }
        hp[i][j] = sum / (double)size1;
      }
    }
  }
  else
  {
    for (it1=s1.begin(); it1 != s1.end(); it1++)
    {
      int j=0;
      vhp[i].resize(size2);
      for (it2=s2.begin(); it2!=s2.end(); it2++)
      {
       rnaduplex(*it1, *it2, vhp[i][j]);
       // correct probability for gap 
       int k=0;
       int l=0;
       if (false)
       {
        if (false)
          vhp[i][j][k][l]=0;
       }
       j++;
      }      
      i++;
    }
    // 各塩基配列対について平均をとる
    // 各要素doubleにしたほうがいいかも
    VVVVF::iterator it_vvhp;
    VVVF::iterator it_vhp;
    int L=vhp[0][0].size();
    int M=vhp[0][0][0].size();
    hp.resize(L+1);// debug temporally
    for(int i=0; i<L; i++)
    {
      hp[i].resize(M+1);
      for (int j=0; j<M; j++)
      {
       double sum=0;
       for (it_vvhp=vhp.begin(); it_vvhp!=vhp.end();it_vvhp++)
       {
         for(it_vhp=(*it_vvhp).begin(); it_vhp!=(*it_vvhp).end(); it_vhp++)
           sum+=(*it_vhp)[i][j];
       }
	     hp[i][j]=sum / (double)(size1*size2);
     }
   }
 }
}


void
RactIP::
rnaduplex(const std::string& s1, const std::string& s2, VVF& hp) const
{
  if (use_pf_duplex_)//true
  {
    Vienna::pf_scale = -1;
    hp.resize(s1.size()+1, VF(s2.size()+1));

    Vienna::pf_duplex(s1.c_str(), s2.c_str());// reading
    for (uint i=0; i!=s1.size(); ++i)
      for (uint j=0; j!=s2.size(); ++j)
        hp[i+1][j+1] = Vienna::pr_duplex[i+1][j+1];
    Vienna::free_pf_duplex();
  }
  else
  {
    hp.clear();
    hp.resize(s1.size()+1, VF(s2.size()+1, 0.0));
    std::string s=s1+s2;
    std::string c(s.size(), 'e');
    Vienna::pf_scale = -1;
    Vienna::cut_point = s1.size()+1;
    Vienna::co_pf_fold(const_cast<char*>(s.c_str()), const_cast<char*>(c.c_str()));// reading
    pair_info* pi = NULL;
#ifdef HAVE_VIENNA20
    Vienna::assign_plist_from_pr(&pi, Vienna::export_co_bppm(), s.size(), th_hy_);
#else
    pi = Vienna::get_plist((pair_info*)malloc(sizeof(*pi)*s.size()), s.size(), th_hy_);
#endif
    for (uint k=0; pi[k].i!=0; ++k)
      if (pi[k].i<Vienna::cut_point && pi[k].j>=Vienna::cut_point && pi[k].p>th_hy_)
        hp[pi[k].i][pi[k].j-Vienna::cut_point+1]=pi[k].p;
    if (pi) free(pi);
    Vienna::free_co_pf_arrays();
    Vienna::cut_point = -1;
  }
}

void
RactIP::
load_from_rip(const char* filename,
              const std::string& s1, const std::string& s2,
              VF& bp1, VI& offset1, VF& bp2, VI& offset2, VVF& hp) const
{
  enum { NONE, TABLE_R, TABLE_S, TABLE_I };
  uint st=NONE;

  uint L1=s1.size();
  bp1.resize((L1+1)*(L1+2)/2);
  offset1.resize(L1+1);
  for (uint i=0; i<=L1; ++i)
    offset1[i] = i*((L1+1)+(L1+1)-i-1)/2;

  uint L2=s2.size();
  bp2.resize((L2+1)*(L2+2)/2);
  offset2.resize(L2+1);
  for (uint i=0; i<=L2; ++i)
    offset2[i] = i*((L2+1)+(L2+1)-i-1)/2;

  hp.resize(L1+1, VF(L2+1));

  std::ifstream is(filename);
  std::string l;
  while (std::getline(is, l))
  {
    if (strncmp(l.c_str(), "Table R:", 8)==0) st=TABLE_R;//
    else if (strncmp(l.c_str(), "Table S:", 8)==0) st=TABLE_S;
    else if (strncmp(l.c_str(), "Table I:", 8)==0) st=TABLE_I;
    else if (st!=NONE && l[0]>='0' && l[0]<='9')
    {
      int i,j;
      double p;
      std::istringstream s(l.c_str());
      s >> i >> j >> p;
      switch (st)
      {
      case TABLE_R://reading
          bp1[offset1[i]+j] = p;
          break;
        case TABLE_S:
          bp2[offset2[L2-j+1]+L2-i+1] = p;
          break;
        case TABLE_I:
          hp[i][L2-j+1] = p;
          break;
        default:
          break;
      }
    }
    else st=NONE;
  }
}

float
RactIP::
solve(TH& s1, TH& s2, std::string& r1, std::string& r2, FoldingEngine<TH>* cf1, FoldingEngine<TH>* cf2)
{
  IP ip(IP::MAX, n_th_);
  VF bp1, bp2;
  VI offset1, offset2;
  VVF hp;
  VVF up1, up2;
  bool enable_accessibility = min_w_>1 && max_w_>=min_w_;


  homfold(s1, bp1, offset1, up1, cf1);
  homfold(s2, bp2, offset2, up2, cf2);
  rnaduplex_hom(s1, s2, hp, WH);
  

#if BP_OUTPUT
  std::string innerbp_name1 = bp_dirname+"/out_bp_1.csv";
  std::ofstream out_bp1(innerbp_name1.c_str());
  VF::iterator it_bp1 = bp1.begin();
  for (it_bp1 = bp1.begin(); it_bp1 < bp1.end(); it_bp1++)
    {
      out_bp1 << (*it_bp1) << ",";
    }
  std::string innerbp_name2 = bp_dirname+"/out_bp2_2.csv";
  std::ofstream out_bp2(innerbp_name2.c_str());
  VF::iterator it_bp2 = bp2.begin();
  for (it_bp2 = bp2.begin(); it_bp2 < bp2.end(); it_bp2++)
    {
      out_bp2 << (*it_bp2) << ",";
    }

  std::string outbp_name = bp_dirname+"/out_hp_2.csv";
  std::ofstream out_hp(outbp_name.c_str());
  VVF::iterator itit_hp = hp.begin();
//  int S = (*(hp.begin())).size();
//  out_hp << 0;
  //for (int i = 1; i < S; i++){
 //   out_hp << "," << i;
//  }
//  out_hp << std::endl;
//  int hp_count = 0;
  for (itit_hp = hp.begin(); itit_hp < hp.end(); itit_hp++)
    {
      //out_hp << hp_count;
      for (VF::iterator it_hp = (*itit_hp).begin(); it_hp < (*itit_hp).end(); it_hp++)
      {
        if (it_hp == (*itit_hp).begin()){
          out_hp << (*it_hp);
        }
        else{
          out_hp << "," << (*it_hp);
        }
     }
     out_hp << std::endl;
//     ++hp_count;
  }
  out_bp1.close();
  out_bp2.close();
  out_hp.close();
  return 0;
#endif


  /**
   * Takashi Matsuda commented out below in 2014
  else if (!rip_file_.empty())
  {
    load_from_rip(rip_file_.c_str(), s1, s2, bp1, offset1, bp2, offset2, hp);
  }
  else if (use_contrafold_)
  {
    contrafold(s1, bp1, offset1, up1);
    contrafold(s2, bp2, offset2, up2);
    //contraduplex(s1, s2, hp);
    rnaduplex(s1, s2, hp);
  }
  else
  {
#if 0
    if (enable_accessibility)
    {
      rnafold(s1, bp1, offset1, up1, max_w_);
      rnafold(s2, bp2, offset2, up2, max_w_);
    }
    else
    {
      rnafold(s1, bp1, offset1);
      rnafold(s2, bp2, offset2);
    }
#else
    rnafold(s1, bp1, offset1, up1, std::max(1, max_w_));
    rnafold(s2, bp2, offset2, up2, std::max(1, max_w_));
#endif
    rnaduplex(s1, s2, hp);
    }**/
  
  // make objective variables with their weights
  VVI x(s1.size(), VI(s1.size(), -1));
  VVI xx(s1.size());
  for (uint j=1; j!=s1.size(); ++j)
  {
    for (uint i=j-1; i!=-1u; --i)
    {
      const float& p=bp1[offset1[i+1]+(j+1)];
      if (p>th_ss_)
      {
        x[i][j] = ip.make_variable(p);
        xx[i].push_back(j);
      }
    }
  }

  VVI y(s2.size(), VI(s2.size(), -1));
  VVI yy(s2.size());
  for (uint j=1; j!=s2.size(); ++j)
  {
    for (uint i=j-1; i!=-1u; --i)
    {
      const float& p=bp2[offset2[i+1]+(j+1)];
      if (p>th_ss_)
      {
        y[i][j] = ip.make_variable(p);
        yy[i].push_back(j);// danger
      }
    }
  }

  VVI z(s1.size(), VI(s2.size(), -1));
  std::vector< std::vector<int> > zz(s1.size());
  for (uint i=0; i!=s1.size(); ++i)
  {
    for (uint j=0; j!=s2.size(); ++j)
    {
      const float& p=hp[i+1][j+1];
      if (p>th_hy_ && ((min_w_==1 && up1[i][0]>th_ac_ && up2[j][0]>th_ac_ )|| min_w_!=1))// error
      {
        z[i][j] = ip.make_variable(alpha_*p+beta_*(up1[i][0]+up2[j][0]));
        zz[i].push_back(j);
      }
    }
  }

  VI v, w;
  std::vector< std::pair<uint,uint> > vv, ww;
  if (enable_accessibility)
  {
    for (uint i=0; i!=up1.size(); ++i)
      for (uint j=min_w_-1; j<up1[i].size(); ++j)
        if (up1[i][j]>th_ac_)
        {
          v.push_back(ip.make_variable(0.0));
          vv.push_back(std::make_pair(i,i+j));
        }

    for (uint i=0; i!=up2.size(); ++i)
      for (uint j=min_w_-1; j<up2[i].size(); ++j)
        if (up2[i][j]>th_ac_)
        {
          w.push_back(ip.make_variable(0.0));
          ww.push_back(std::make_pair(i,i+j));
        }
  }

  ip.update();

  // constraint 1: each a_i is paired with at most one base
  for (uint i=0; i!=s1.size(); ++i)
  {
    int row = ip.make_constraint(IP::UP, 0, 1);
    if (enable_accessibility)
    {
      int row_ac = ip.make_constraint(IP::LO, 0, 0);
      for (uint j=0; j!=v.size(); ++j)
        if (vv[j].first<=i && i<=vv[j].second)
        {
          ip.add_constraint(row, v[j], 1);
          ip.add_constraint(row_ac, v[j], 1);
        }
      for (uint j=0; j!=s2.size(); ++j)
        if (z[i][j]>=0)
          ip.add_constraint(row_ac, z[i][j], -1);
    }
    else
    {
      for (uint j=0; j!=s2.size(); ++j)
        if (z[i][j]>=0)
          ip.add_constraint(row, z[i][j], 1);
    }

    for (uint j=0; j<i; ++j)
      if (x[j][i]>=0)
        ip.add_constraint(row, x[j][i], 1);
    for (uint j=i+1; j<s1.size(); ++j)
      if (x[i][j]>=0)
        ip.add_constraint(row, x[i][j], 1);
  }

  if (enable_accessibility && !allow_concat_)
  {
    for (uint i=0; i!=v.size(); ++i)
      for (uint j=i+1; j!=v.size(); ++j)
        if (vv[i].second+1==vv[j].first || vv[j].second+1==vv[i].first)
        {
          int row = ip.make_constraint(IP::UP, 0, 1);
          ip.add_constraint(row, v[i], 1);
          ip.add_constraint(row, v[j], 1);
        }
  }

  // constraint 2: each b_i is paired with at most one
  for (uint i=0; i!=s2.size(); ++i)
  {
    int row = ip.make_constraint(IP::UP, 0, 1);
    if (enable_accessibility)
    {
      int row_ac = ip.make_constraint(IP::LO, 0, 0);
      for (uint j=0; j!=w.size(); ++j)
        if (ww[j].first<=i && i<=ww[j].second)
        {
          ip.add_constraint(row, w[j], 1);
          ip.add_constraint(row_ac, w[j], 1);
        }
      for (uint j=0; j!=s1.size(); ++j)
        if (z[j][i]>=0)
          ip.add_constraint(row_ac, z[j][i], -1);
    }
    else
    {
      for (uint j=0; j!=s1.size(); ++j)
        if (z[j][i]>=0)
          ip.add_constraint(row, z[j][i], 1);
    }

    for (uint j=0; j<i; ++j)
      if (y[j][i]>=0)
        ip.add_constraint(row, y[j][i], 1);
    for (uint j=i+1; j<s2.size(); ++j)
      if (y[i][j]>=0)
        ip.add_constraint(row, y[i][j], 1);
  }

  if (enable_accessibility && !allow_concat_)
  {
    for (uint i=0; i!=w.size(); ++i)
      for (uint j=i+1; j!=w.size(); ++j)
        if (ww[i].second+1==ww[j].first || ww[j].second+1==ww[i].first)
        {
          int row = ip.make_constraint(IP::UP, 0, 1);
          ip.add_constraint(row, w[i], 1);
          ip.add_constraint(row, w[j], 1);
        }
  }

  // constraint 3: disallow external pseudoknots
  for (uint i=0; i<zz.size(); ++i)
  {
    for (uint k=i+1; k<zz.size(); ++k)
    {
      for (uint p=0; p<zz[i].size(); ++p)
      { 
        uint j=zz[i][p];
        for (uint q=0; q<zz[k].size(); ++q)
        {
          uint l=zz[k][q];
          if (j<l)
          {
            int row = ip.make_constraint(IP::UP, 0, 1);
            ip.add_constraint(row, z[i][j], 1);
            ip.add_constraint(row, z[k][l], 1);
          }
        }
      }
    }
  }

  if (in_pk_)
  {
    // constraint 4: disallow internal pseudoknots in a
    for (uint i=0; i<xx.size(); ++i)
    {
      for (uint p=0; p<xx[i].size(); ++p)
      {
        uint j=xx[i][p];
        for (uint k=i+1; k<j; ++k)
        {
          for (uint q=0; q<xx[k].size(); ++q)
          {
            uint l=xx[k][q];
            if (j<l)
            {
              int row = ip.make_constraint(IP::UP, 0, 1);
              ip.add_constraint(row, x[i][j], 1);
              ip.add_constraint(row, x[k][l], 1);
            }
          }
        }
      }
    }

    // constraint 5: disallow internal pseudoknots in b
    for (uint i=0; i<yy.size(); ++i)
    {
      for (uint p=0; p<yy[i].size(); ++p)
      {
        uint j=yy[i][p];
        for (uint k=i+1; k<j; ++k)
        {
          for (uint q=0; q<yy[k].size(); ++q)
          {
            uint l=yy[k][q];
            if (j<l)
            {
              int row = ip.make_constraint(IP::UP, 0, 1);
              ip.add_constraint(row, y[i][j], 1);
              ip.add_constraint(row, y[k][l], 1);
            }
          }
        }
      }
    }
  }

  if (stacking_constraints_)
  {
    // upstream of s1
    for (uint i=0; i<s1.size(); ++i)
    {
      int row = ip.make_constraint(IP::LO, 0, 0);
      for (uint j=0; j<i; ++j)
        if (x[j][i]>=0)
          ip.add_constraint(row, x[j][i], -1);
      if (i>0)
        for (uint j=0; j<i-1; ++j)
          if (x[j][i-1]>=0)
            ip.add_constraint(row, x[j][i-1], 1);
      if (i+1<s1.size())
        for (uint j=0; j<i+1; ++j)
          if (x[j][i+1]>=0)
            ip.add_constraint(row, x[j][i+1], 1);
    }

    // downstream of s1
    for (uint i=0; i<s1.size(); ++i)
    {
      int row = ip.make_constraint(IP::LO, 0, 0);
      for (uint j=i+1; j<s1.size(); ++j)
        if (x[i][j]>=0)
          ip.add_constraint(row, x[i][j], -1);
      if (i>0)
        for (uint j=i; j<s1.size(); ++j)
          if (x[i-1][j]>=0)
            ip.add_constraint(row, x[i-1][j], 1);
      if (i+1<s1.size())
        for (uint j=i+2; j<s1.size(); ++j)
          if (x[i+1][j]>=0)
            ip.add_constraint(row, x[i+1][j], 1);
    }

    // upstream of s2
    for (uint i=0; i<s2.size(); ++i)
    {
      int row = ip.make_constraint(IP::LO, 0, 0);
      for (uint j=0; j<i; ++j)
        if (y[j][i]>=0)
          ip.add_constraint(row, y[j][i], -1);
      if (i>0)
        for (uint j=0; j<i-1; ++j)
          if (y[j][i-1]>=0)
            ip.add_constraint(row, y[j][i-1], 1);
      if (i+1<s2.size())
        for (uint j=0; j<i+1; ++j)
          if (y[j][i+1]>=0)
            ip.add_constraint(row, y[j][i+1], 1);
    }

    // downstream of s1
    for (uint i=0; i<s2.size(); ++i)
    {
      int row = ip.make_constraint(IP::LO, 0, 0);
      for (uint j=i+1; j<s2.size(); ++j)
        if (y[i][j]>=0)
          ip.add_constraint(row, y[i][j], -1);
      if (i>0)
        for (uint j=i; j<s2.size(); ++j)
          if (y[i-1][j]>=0)
            ip.add_constraint(row, y[i-1][j], 1);
      if (i+1<s2.size())
        for (uint j=i+2; j<s2.size(); ++j)
          if (y[i+1][j]>=0)
            ip.add_constraint(row, y[i+1][j], 1);
    }

    // for s2
    for (uint i=0; i<s2.size(); ++i)
    {
      int row = ip.make_constraint(IP::LO, 0, 0);
      for (uint j=0; j<s1.size(); ++j)
        if (z[j][i]>=0)
          ip.add_constraint(row, z[j][i], -1);
      if (i>0)
        for (uint j=0; j<s1.size(); ++j)
          if (z[j][i-1]>=0)
            ip.add_constraint(row, z[j][i-1], 1);
      if (i+1<s2.size())
        for (uint j=0; j<s1.size(); ++j)
          if (z[j][i+1]>=0)
            ip.add_constraint(row, z[j][i+1], 1);
    }

    // for s1
    for (uint i=0; i<s1.size(); ++i)
    {
      int row = ip.make_constraint(IP::LO, 0, 0);
      for (uint j=0; j<s2.size(); ++j)
        if (z[i][j]>=0)
          ip.add_constraint(row, z[i][j], -1);
      if (i>0)
        for (uint j=0; j<s2.size(); ++j)
          if (z[i-1][j]>=0)
            ip.add_constraint(row, z[i-1][j], 1);
      if (i+1<s1.size())
        for (uint j=0; j<s2.size(); ++j)
          if (z[i+1][j]>=0)
            ip.add_constraint(row, z[i+1][j], 1);
    }
  }

  // execute optimization
  std::cout << "execute optimization" << std::endl;
  float ea = ip.solve();
  std::cout << "finished optimization" << std::endl;
  // build the resultant structure
  r1.resize(s1.size());
  r2.resize(s2.size());
  std::fill(r1.begin(), r1.end(), '.');
  std::fill(r2.begin(), r2.end(), '.');
  for (uint i=0; i!=s1.size(); ++i)
  {
    for (uint j=0; j!=s2.size(); ++j)
    {
      if (z[i][j]>=0 && ip.get_value(z[i][j])>0.2)
      {

        r1[i]='['; r2[j]=']';
      }
    }
  }
  if (in_pk_)
  {
    for (uint i=0; i<s1.size(); ++i)
    {
      for (uint j=i+1; j<s1.size(); ++j)
      {
        if (x[i][j]>=0 && ip.get_value(x[i][j])>0.2)
        {
          assert(r1[i]=='.'); assert(r1[j]=='.');
          r1[i]='('; r1[j]=')';
        }
      }
    }
    for (uint i=0; i<s2.size(); ++i)
    {
      for (uint j=i+1; j<s2.size(); ++j)
      {
        if (y[i][j]>=0 && ip.get_value(y[i][j])>0.2)
        {
          assert(r2[i]=='.'); assert(r2[j]=='.');
          r2[i]='('; r2[j]=')';
        }
      }
    }
  }

#if 0
  std::cout << "v: ";
  for (uint i=0; i!=v.size(); ++i)
  {
    if (ip.get_value(v[i])>0.5)
      std::cout << "(" << vv[i].first << "," << vv[i].second << ","
                << up1[vv[i].first][vv[i].second-vv[i].first]<< "), ";
  }
  std::cout << std::endl;

  std::cout << "w: ";
  for (uint i=0; i!=w.size(); ++i)
  {
    if (ip.get_value(w[i])>0.5)
      std::cout << "(" << ww[i].first << "," << ww[i].second << ","safa-
                << up2[ww[i].first][ww[i].second-ww[i].first]<< "), ";
  }
  std::cout << std::endl;
#endif
  
  return ea;
}

RactIP&
RactIP::
parse_options(int& argc, char**& argv)
{
  gengetopt_args_info args_info;// Gnu gengetopt
  if (cmdline_parser(argc, argv, &args_info)!=0) exit(1);
  // the last filepointer should be taken by centroid_alifold section.
  alpha_ = args_info.alpha_arg;
  beta_ = args_info.beta_arg;
  th_ss_ = args_info.fold_th_arg;
  th_hy_ = args_info.hybridize_th_arg;
  th_ac_ = args_info.acc_th_arg;
  mix_w.push_back(args_info.mix_weight_arg);//float mix_weight_arg
  max_w_ = args_info.max_w_arg;
  min_w_ = args_info.min_w_arg;
  WH = args_info.hyb_mix_w_arg;
  enable_zscore_ = args_info.zscore_arg;
  num_shuffling_ = args_info.num_shuffling_arg;
  seed_ = args_info.seed_arg;
  in_pk_ = args_info.no_pk_flag==0;
  use_contrafold_ = args_info.mccaskill_flag==0;
  if (args_info.engine_seq_arg != NULL) engine.push_back(args_info.engine_seq_arg);
  if (args_info.engine_aln_arg != NULL) engine_a.push_back(args_info.engine_aln_arg);
  //models = args_info.models_arg // not yet implemented
  //use_pf_duplex_ = args_info.pf_duplex_flag;
  stacking_constraints_ = args_info.allow_isolated_flag==0;
  //allow_concat_ = args_info.allow_concat_flag;
  //run_with_modena_ = args_info.modena_flag;
  n_th_ = 1; // args_info.n_th_arg;
  if (args_info.rip_given) rip_file_ = args_info.rip_arg;
  show_energy_ = args_info.show_energy_flag==1;
  if (args_info.param_file_given) param_file_ = args_info.param_file_arg;

  if (args_info.inputs_num==0 ||
      (min_w_!=0 && max_w_!=0 && (min_w_>max_w_ || use_contrafold_)))
  {
    cmdline_parser_print_help();
    cmdline_parser_free(&args_info);
    exit(1);
  }

  if (args_info.inputs_num>=1){
    fa1_ = args_info.inputs[0];

  }
  if (args_info.inputs_num>=2){
    aln1_ = args_info.inputs[1];

  }
  if (args_info.inputs_num>=3){
    fa2_ = args_info.inputs[2];

  }
  if (args_info.inputs_num>=4){
    aln2_ = args_info.inputs[3];

  }

  bp_dirname =  args_info.output_dir_arg;  

  cmdline_parser_free(&args_info);
  return *this;
}

// static
void
RactIP::
calculate_energy(const std::string s1, const std::string& s2,
                 const std::string r1, const std::string& r2,
                 float& e1, float& e2, float& e3)
{
#ifdef HAVE_VIENNA20
  e1=Vienna::energy_of_structure(s1.c_str(), r1.c_str(), -1);
  e2=Vienna::energy_of_structure(s2.c_str(), r2.c_str(), -1);
#else
  Vienna::eos_debug = -1;
  e1=Vienna::energy_of_struct(s1.c_str(), r1.c_str());
  e2=Vienna::energy_of_struct(s2.c_str(), r2.c_str());
#endif
  std::string ss(s1+"NNN"+s2);

  std::string r1_temp(r1);
  for (std::string::iterator x=r1_temp.begin(); x!=r1_temp.end(); ++x)
  {
    switch (*x)
    {
      case '(': case ')': *x='.'; break;
      case '[': *x='('; break;
      default: break;
    }
  }

  std::string r2_temp(r2);
  for (std::string::iterator x=r2_temp.begin(); x!=r2_temp.end(); ++x)
  {
    switch (*x)
    {
      case '(': case ')': *x='.'; break;
      case ']': *x=')'; break;
      default: break;
    }
  }

  std::string rr(r1_temp+"..."+r2_temp);
  //std::cout << ss << std::endl << rr << std::endl;
#ifdef HAVE_VIENNA20
  e3=Vienna::energy_of_structure(ss.c_str(), rr.c_str(), -1);
#else
  Vienna::eos_debug = -1;
  e3=Vienna::energy_of_struct(ss.c_str(), rr.c_str());
#endif
}

int
RactIP::
run()
{
  // set the energy parameters
  copy_boltzmann_parameters();
  if (!param_file_.empty())
    Vienna::read_parameter_file(param_file_.c_str());
  
  // Input Data
  Fasta fa1, fa2;
  std::vector<std::string> homs1;
  std::vector<std::string> homs2;
  if (!fa1_.empty() && !fa2_.empty() && !aln1_.empty() && !aln2_.empty())
  {
    //std::cout<<"both fa1 and fa2 are not empty"<<std::endl;
    std::list<Fasta> l_fa1, l_fa2;
    if (Fasta::load(l_fa1, fa1_.c_str())==0)// l_fa1にfa1_の中身を格納する。
      throw (fa1_+": Format error").c_str();
    if (Fasta::load(l_fa2, fa2_.c_str())==0)
      throw (fa2_+": Format error").c_str();
    fa1=l_fa1.front();
    fa2=l_fa2.front();
    

    if (aln1_ != "") {
      BOOST_SPIRIT_CLASSIC_NS::file_iterator<> fi1(aln1_.c_str());
      if (!fi1) {
       perror(aln1_.c_str());
       return 1;
     }
     while (1) {
       Fasta fa;
       if (fa.load(fi1)) {
         homs1.push_back (fa.seq());
       } else break;
     }
   }
    if (aln2_ != "") {
      BOOST_SPIRIT_CLASSIC_NS::file_iterator<> fi2(aln2_.c_str());
      if (!fi2) {
       perror(aln2_.c_str());
       return 1;
     }
     while (1) {
       Fasta fa;
       if (fa.load(fi2)) {
         homs2.push_back (fa.seq());
       } else break;
     }
   }


    std::replace(fa1.seq().begin(), fa1.seq().end(), 't', 'u');
    std::replace(fa1.seq().begin(), fa1.seq().end(), 'T', 'U');
    std::vector<std::string>::iterator it_homs1;
    for (it_homs1 = homs1.begin(); it_homs1 != homs1.end(); it_homs1++){
      std::replace((*it_homs1).begin(), (*it_homs1).end(), '-', '.');
      std::replace((*it_homs1).begin(), (*it_homs1).end(), 't', 'u');
      std::replace((*it_homs1).begin(), (*it_homs1).end(), 'T', 'U');
    }
    std::replace(fa2.seq().begin(), fa2.seq().end(), 't', 'u');
    std::replace(fa2.seq().begin(), fa2.seq().end(), 'T', 'U');
    std::vector<std::string>::iterator it_homs2;
    for (it_homs2 = homs2.begin(); it_homs2 != homs2.end(); it_homs2++){
      std::replace((*it_homs2).begin(), (*it_homs2).end(), '-', '.');
      std::replace((*it_homs2).begin(), (*it_homs2).end(), 't', 'u');
      std::replace((*it_homs2).begin(), (*it_homs2).end(), 'T', 'U');
    }
  }


  /**
  else if (!fa1_.empty())
  {
    std::list<Aln> l1;
    if (Aln::load(l1, fa1_.c_str())<2)
      throw (fa1_+": Format error").c_str();
    std::list<Aln>::const_iterator x=l1.begin();
    fa1=*(x++);
    fa2=*(x++);
  }

  else { throw "unreachable"; }
  **/
  
  std::vector<FoldingEngine<Aln>*> cf_list(engine.size(),NULL);
  if (engine.empty())
    {// Default setting for Inference Engine
#ifdef HAVE_LIB
      engine.push_back("McCaskill");
#else
      engine.push_back("CONTRAfold");
#endif
    }
  /**
  std::vector<std::string> homs;
  if (hom_seqs != "") {
    BOOST_SPIRIT_CLASSIC_NS::file_iterator<> fi2(hom_seqs.c_str());
    if (!fi2) {
      perror(hom_seqs.c_str());
      return 1;
    }
    while (1) {
      Fasta fa;
      if (fa.load(fi2)) {
	homs.push_back (fa.seq());
      } else break;
    }
  }
  **/
  FoldingEngine<TH>* cf1=NULL;
  FoldingEngine<TH>* cf2=NULL;
  std::vector<FoldingEngine<TH>*> cf_list_1(engine.size(), NULL);
  std::vector<FoldingEngine<TH>*> cf_list_2(engine.size(), NULL);

  // tmp code
  uint max_bp_dist = 1000;
  char* param_tmp = NULL;
  std::string param_tmp2;
  
  //cf_list_1[0] = new McCaskillHomModel(engine_a[0], false, max_bp_dist);
  //cf_list_2[0] = new McCaskillHomModel(engine_a[0], false, max_bp_dist);
  //cf1 = new McCaskillHomModel(engine_a[0], false, max_bp_dist);
  //cf2 = new McCaskillHomModel(engine_a[0], false, max_bp_dist);
  cf1 = new CONTRAfoldHomModel(param_tmp2, engine_a[0], true, max_bp_dist);// canonical=trueにするべきか
  cf2 = new CONTRAfoldHomModel(param_tmp2, engine_a[0], true, max_bp_dist);

  // とりあえずMcCaskillモデルで動かしてみる。以下はcentroidhomfoldからのコピーコード。
  /**
  for (uint i=0; i!=engine.size(); ++i)
  {
    if (engine[i]=="CONTRAfold")
    {
      if (gamma.empty()) gamma.push_back(vm.count("mea") ? 6.0 : 2.0);
      cf_list[i] = new CONTRAfoldHomModel(param, engine_a[0], !vm.count("noncanonical"), max_bp_dist, seed, vm.count("mea"));
    }
#ifdef HAVE_LIBRNA
    else if (engine[i]=="McCaskill")
    {
      if (gamma.empty()) gamma.push_back(vm.count("mea") ? 6.0 : 1.0);
      cf_list[i] = new McCaskillHomModel(engine_a[0], !vm.count("noncanonical"), max_bp_dist,
					 param.empty() ? NULL : param.c_str(),
                                         seed, vm.count("mea"));
    }
#endif
    else
    {
      std::cerr << engine[i] << std::endl;
      throw std::logic_error("unsupported inference engine");
    }
  }

  
  if (engine.size()==1)
    cf=cf_list[0];
  **/
  //std::vector<FoldingEngine<Aln>*> cf_list(engine.size(), NULL);
  //std::vector<FoldingEngine<std::string>*> src_list(engine.size(), NULL);
   /**
  for (uint i=0; i!=engine.size(); ++i)
  {
    if (engine[i]=="CONTRAfold")
    {
      if (gamma.empty()) gamma.push_back(vm.count("mea") ? 6.0 : 4.0);
      src_list[i] = new CONTRAfoldModel(param, !vm.count("noncanonical"), max_bp_dist, seed, vm.count("mea"));
      cf_list[i] = new AveragedModel(src_list[i], max_bp_dist, vm.count("mea"));
    }
    else if (engine[i]=="CONTRAfoldM")
    {
      if (gamma.empty()) gamma.push_back(vm.count("mea") ? 6.0 : 4.0);
      cf_list[i] = new CONTRAfoldMultiModel(param, !vm.count("noncanonical"), max_bp_dist, seed, vm.count("mea"));
    }
#ifdef HAVE_LIBRNA
    else if (engine[i]=="McCaskill")
    {
      if (gamma.empty()) gamma.push_back(vm.count("mea") ? 6.0 : 2.0);
      src_list[i] = new McCaskillModel(!vm.count("noncanonical"), max_bp_dist,
                                       param.empty() ? NULL : param.c_str(),
                                       seed, vm.count("mea"));
      cf_list[i] = new AveragedModel(src_list[i], max_bp_dist, vm.count("mea"));
    }
    else if (engine[i]=="Alifold")
    {
      if (gamma.empty()) gamma.push_back(vm.count("mea") ? 6.0 : 2.0);
      cf_list[i] = new AliFoldModel(!vm.count("noncanonical"), max_bp_dist,
                                    param.empty() ? NULL : param.c_str(),
                                    seed, vm.count("mea"));
    }
#endif
    else if (engine[i]=="pfold")
    {
      if (gamma.empty()) gamma.push_back(vm.count("mea") ? 6.0 : 1.0);
      std::string pfold_bin_dir(getenv("PFOLD_BIN_DIR") ? getenv("PFOLD_BIN_DIR") : ".");
      std::string awk_bin(getenv("AWK_BIN") ? getenv("AWK_BIN") : "mawk");
      std::string sed_bin(getenv("SED_BIN") ? getenv("SED_BIN") : "sed");
      cf_list[i] = new PfoldModel<Aln>(pfold_bin_dir, awk_bin, sed_bin, vm.count("mea"));
    }
    else if (engine[i]=="AUX")
    {
      if (gamma.empty()) gamma.push_back(vm.count("mea") ? 6.0 : 1.0);
      src_list[i] = new AuxModel(model, vm.count("mea"));
      cf_list[i] = new AveragedModel(src_list[i], 0, vm.count("mea"));
    }
    else
    {
      throw std::logic_error("unsupported engine");
    }
  }
  **/
  //  std::vector<std::pair<FoldingEngine<Aln>*,float> > models;


  /**
  // ここでモデルの設定を行う。
  // ゆくゆくは上で条件分岐できるように、コマンドラインを設計する。FUTURE work
  std::vector<FoldingEngine<Aln>*> cf_list(1, NULL);
  bool mea_bool = false;
  uint max_bp_dist = 0;
  uint seed_model = 0;
  
  cf_list[0] = new McCaskillHomModel(engine_a[0], !vm.count("noncanonical"), max_bp_dist,
					 param.empty() ? NULL : param.c_str(),
                                         seed, vm.count("mea"));
  uint enginenumber=2;// tmp
  {
    // if (gamma.empty()) gamma.push_back(vm.count("mea") ? 6.0 : 1.0);
    for (uint i=0; i!=enginenumber ; ++i)
    {
      if (engine.size()!=mix_w.size())
        models.push_back(std::make_pair(cf_list[i], 1.0));
      else
        models.push_back(std::make_pair(cf_list[i], mix_w[i]));
    }
    //cf = new MixtureModel<Aln>(models, vm.count("mea"));
    
  }
  **/

  // predict the interation
  // template std::vector<std::string> change_l_v<std::string> (std::list<std::string>&);
  std::string r1, r2;
  const std::string fa1_seq=fa1.seq();
  const std::string fa2_seq=fa2.seq();
  TH th1_ = TH(fa1_seq, homs1);
  TH th2_ = TH(fa2_seq, homs2);
  float ea = solve(th1_, th2_, r1, r2, cf1, cf2);

  // diplay the result
  const std::string fa1_name=fa1.name();
  std::cout << ">" << fa1_name << std::endl
	    << fa1_seq << std::endl;

  std::cout<< r1 << std::endl;
  
  const std::string fa2_name=fa2.name();
  std::cout << ">" << fa2_name << std::endl
		<< fa2_seq << std::endl;
 
  std::cout<< r2 << std::endl;
  
  //////////////////////// by Takashi Matsuda ///////////////////////////////
  //// future work: separate above code ///////
  //// future work: Revive the code below
  
  /**
  // show energy of the joint structure
  if (show_energy_ || enable_zscore_==1 || enable_zscore_==2 || enable_zscore_==12)
  {
    float e1, e2, e3;
    calculate_energy(fa1.seq(), fa2.seq(), r1, r2, e1, e2, e3);
    if (show_energy_)
      std::cout << "(E: S1=" << e1 << ", "
                << "S2=" << e2 << ", "
                << "H=" << e3 << ", "
                << "JS=" << e1+e2+e3 << ", "
                << "EA=" << ea << ")" << std::endl;

    if (enable_zscore_==1 || enable_zscore_==2 || enable_zscore_==12)
    {
      float sum=0.0, sum2=0.0;
      std::string s1(fa1.seq());
      std::string s2(fa2.seq());
      if (seed_==0)
      {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	seed_ = (unsigned long) tv.tv_sec;
      }
      srandom(seed_);
      uShuffle::set_randfunc((uShuffle::randfunc_t) random);
      for (uint i=0; i!=num_shuffling_; ++i)
      {
        if (enable_zscore_==1 || enable_zscore_==12)
          uShuffle::shuffle(fa1.seq().c_str(), &s1[0], fa1.seq().size(), 2);
        if (enable_zscore_==2 || enable_zscore_==12)
          uShuffle::shuffle(fa2.seq().c_str(), &s2[0], fa2.seq().size(), 2);
        solve(s1, s2, r1, r2);//
        float ee1, ee2, ee3;
        calculate_energy(s1, s2, r1, r2, ee1, ee2, ee3);
        float ee=ee1+ee2+ee3;
        sum += ee;
        sum2 += ee*ee;
      }
      float m=sum/num_shuffling_;
      float v=sum2/num_shuffling_-m*m;
      if (v<0.0) v=0.0;
      if (v!=0)
        std::cout << "z-score: " << (e1+e2+e3-m)/sqrt(v) << std::endl;
      else
        std::cout << "z-score: " << 0.0 << std::endl;
    }
  }
  **/
  return 0;
}

int
main(int argc, char* argv[])
{
  try
  {
    RactIP ractip;
    //std::cout<< "start" << std::endl;
    return ractip.parse_options(argc, argv).run();
  }
  catch (const char* msg)
  {
    std::cout << msg << std::endl;
  }
  catch (std::logic_error err)
  {
    std::cout << err.what() << std::endl;
  }
  return 1;
}
