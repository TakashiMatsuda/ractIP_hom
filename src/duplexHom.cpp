/* 
*
* Copyright (C) 2014 Takashi Matsuda
*
* This file is part of RactIP_with_homologous_sequences.
* GPL licence version 3.
*
*/


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

#include "duplexHom.h"
// include guard maybe necessary
namespace Vienna {
  extern "C" {
  #include <ViennaRNA/fold.h>
  #include <ViennaRNA/fold_vars.h>
  #include <ViennaRNA/part_func.h>
  #include <ViennaRNA/part_func_up.h>
  #include <ViennaRNA/part_func_co.h>
  #include <ViennaRNA/utils.h>
  #include "pf_duplex.h"
  };
};

#include "centroidalifold/util.h"
#include "centroidalifold/probconsRNA/probcons.h"

#include <omp.h>

RNAduplexHommodel::
RNAduplexHommodel(const std::string& engine_a){
	engine_a_ = engine_a;
}

RNAduplexHommodel::
~RNAduplexHommodel(){
}


// すべてのアラインメントパターンを計算して、それぞれについて
// 二本の間のアラインメント間の座標対塩基対確率を計算し、
// アラインメントパターンについて総和をとる。
//　結果として、各二次元座標の格子点ごとに、座標対が存在する確率が計算される。
VVF
RNAduplexHommodel::
calculate_posterior(const TH &th1, const TH &th2, double wh){
  const std::vector<std::string>& hom1 = th1.second;
  const std::vector<std::string>& hom2 = th2.second;

  VVF res;
  res.resize(th1.first.length()+1);
  for (VVF::iterator itr = res.begin(); itr != res.end(); ++itr){
    (*itr).resize(th2.first.length()+1);
  }

  VVVVF dup_matrix;
  dup_matrix.resize(hom1.size());
  for (VVVVF::iterator itr = dup_matrix.begin(); itr != dup_matrix.end(); ++itr){
    (*itr).resize(hom2.size());
  }
  VVVF align_matrix1;
  VVVF align_matrix2;

  dup_matrix=rnaduplex_hom(th1, th2);
  align_matrix1=align_v(th1, 0.000001);
  align_matrix2=align_v(th2, 0.000001);

  VVF no_hom_hp;
  rnaduplex(th1.first, th2.first, no_hom_hp);

  long double tmp = 0;
#if _OPENMP
  std::cout << "thread_num= " << omp_get_max_threads() << std::endl;
#endif
#pragma omp parallel
{
  #pragma omp for
  for (int i = 0; i < res.size()-1; ++i){
    std::cout << "HP calculation: " << i << " / " << res.size()-1 << std::endl;
    for (int k = 0; k < res[i].size()-1; ++k){
      tmp = 0;
      for (uint xi = 0; xi < hom1.size(); ++xi){
        for (uint eta = 0; eta < hom2.size(); ++eta){
          for (uint j = 0; j < hom1[xi].size(); ++j){
            for (uint l = 0; l < hom2[eta].size(); ++l){
              res[i+1][k+1] += (double)dup_matrix[xi][eta][j+1][l+1]*(double)align_matrix1[xi][i][j]*(double)align_matrix2[eta][k][l];
            }
          }
        }
      }
      res[i+1][k+1] = res[i+1][k+1] / ((double)hom1.size() * (double)hom2.size());
      //res[i+1][k+1] = tmp;
      // whの比にしたがって、no_hom_hp(hom未考慮のhp)とres(考慮したhp)を混合
      res[i+1][k+1] = wh * (double) no_hom_hp[i+1][k+1] + (1 - wh) * (double) res[i+1][k+1];
    }
  }
}

  return res;
}


// alignment probability
VVVF
RNAduplexHommodel::
align_v(const TH& th, double min_aln){
  const std::string& seq = th.first;
  const std::vector<std::string>& hom = th.second;
  VVVF res;
  res.resize(hom.size());
  PROBCONS::Probcons* pc = new PROBCONS::Probcons();
  for (uint n=0; n<hom.size(); ++n) {
    //pc->ComputePosterior(seq, hom[n], min_aln);
    // 型はfloat* 
    // float* Probcson::Impl::ComputePosterior(...)を読んで、なぜ塩基ごとに確率が計算できているかを確認するところから。
    // 確認できたら、それを参考にして、塩基ごとの確率を計算するコードをここに実装する。
    res[n] = computeposterior(*pc, seq, hom[n]);// probcons alignment prob
  }
  return res;
}

// outputtting alignment probability by base
VVF
RNAduplexHommodel::
computeposterior(PROBCONS::Probcons& pc, const std::string& seq1, const std::string& seq2){
  std::vector<float> ap;
  double th = 0.001;
  pc.ComputePosterior(seq1, seq2, ap, th);

  // transform from float* to Vector<Vector<float>>
  int L = seq1.length();
  int M = seq2.length();
  VVF res;
  res.resize(seq1.length());
  for (int i = 0; i < L; ++i){
    res[i].resize(M);
    for (int j = 0; j < M; ++j){
      res[i][j] = ap[(M+1)*(i+1) + j+1];
    }
  }
  return res;
}


// caution::secondの中にfirstも入れておくこと
VVVVF
RNAduplexHommodel::
rnaduplex_hom(const TH& th1, const TH& th2){
  const std::vector<std::string>& hom1 = th1.second;
  const std::vector<std::string>& hom2 = th2.second;

  VVVVF dup_matrix;
  dup_matrix.resize(hom1.size());
  for (uint xi = 0; xi < hom1.size(); ++xi){
    dup_matrix[xi].resize(hom2.size());
    for (uint eta = 0; eta < hom2.size(); ++eta){
      rnaduplex(hom1[xi], hom2[eta], dup_matrix[xi][eta]);
    }
  }
  return dup_matrix;
}



void
RNAduplexHommodel::
rnaduplex(const std::string& s1, const std::string& s2, VVF& hp) const
{
	Vienna::pf_scale = -1;
	hp.resize(s1.size()+1, VF(s2.size()+1));
	Vienna::pf_duplex(s1.c_str(), s2.c_str());
    for (uint i=0; i!=s1.size(); ++i)
      for (uint j=0; j!=s2.size(); ++j)
        hp[i+1][j+1] = Vienna::pr_duplex[i+1][j+1];
    Vienna::free_pf_duplex();
}