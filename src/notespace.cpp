{
    if (gamma.empty()) gamma.push_back(vm.count("mea") ? 6.0 : 1.0);
    std::vector<std::pair<FoldingEngine<Aln>*,float> > models;//watching
    for (uint i=0; i!=engine.size(); ++i)
    {
      if (engine.size()!=mix_w.size())
        models.push_back(std::make_pair(cf_list[i], 1.0)e);
      else
        models.push_back(std::make_pair(cf_list[i], mix_w[i]));
    }
    cf = new MixtureModel<Aln>(models, vm.count("mea"));
}

// 上が注目しているcentroid_alifoldのbpmatrix計算部分

// 下がAlifold.


void
RactIP::
alifold(const std::string& seq, VF& bp, VI& offset, VVF& up, engine, mix_w, cf_list, vm) const
{
  std::vector<std::string> engine;
  if (gamma.empty()) gamma.push_back(vm.count("mea") ? 6.0 : 1.0);
  std::vector<std::pair<FoldingEngine<Aln>*,float> > models;
  for (uint i=0; i!=engine.size(); ++i)
    {
      if (engine.size()!=mix_w.size())
	models.push_back(std::make_pair(cf_list[i], 1.0)e);
      else
	models.push_back(std::make_pair(cf_list[i],mix_w[i]));
    }
  cf=new MixtureModel<Aln>(models,vm.count("mea");

			   /**
			    * engine.size()
			    **/

  // upの計算部分、そのまま
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

// 下が注目しているcontrafoldのbpmatrix 投げられている
// reading
void
RactIP::
contrafold(const std::string& seq, VF& bp, VI& offset, VVF& up) const
{
  /**
   *この作り方は"contrafold/**Engine.hpp"にあわせているため
   **/
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
  en.GetPosterior(0, bp, offset);// alifold must implements this function type

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


// 下が注目しているRactIPのbpmatrix を投げてその結果を得る部分
    contrafold(s1, bp1, offset1, up1);//reading

/**
 * 
 **/

// Alifoldを用いる場合、次のように投げる。
/**
 *　上での計算結果と同じようにbp1,offset1,up1を返す。
 **/
  if(use_alifold){
  alifold(s1, bp1, offset1, up1);// define alifold later
  alifold(s2, bp2, offset2, up2);
  
  rnaduplex(s1,s2,hp)// maybe change this
  }

