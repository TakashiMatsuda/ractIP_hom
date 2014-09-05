#include "duplexHom.h"

RNAduplexHommodel::
RNAduplexHommodel(const std::string& engine_a){
	engine_a_ = engine_a;
}


// すべてのアラインメントパターンを計算して、それぞれについて
// 二本の間のアラインメント間の座標対塩基対確率を計算し、
// アラインメントパターンについて総和をとる。
//　結果として、各二次元座標の格子点ごとに、座標対が存在する確率が計算される。
// 詳しくは降機後、ノートを見ながら実装しよう。
void
RNAduplexHommodel::
calculate_posterior(const TH&s1, const TH&s2, double min_aln){
	std::vector<Alignment> a1_v = align_v(s1, min_aln);
	std::vector<Alignment> a2_v = align_v(s2, min_aln);
	VVF hp;
	for (std::vector<Alignment>::iterator itr_a1 = a1_v.begin(); itr_a1 != a1_v.end(); itr_a1++){
		for (std::vector<Alignment>::iterator itr_a2 = a2_v.begin(); itr_a2 != a2_v.end(); itr_a2++){
			VVF hp_tmp;
			aln_duplex(*itr_a1, *itr_a2, hp_tmp);
			for (int i = 0; i < hp.size(); i++){
				for (int j = 0; j < hp[0].size(); j++){
					hp[i][j] += hp_tmp[i][j] * (*itr_a1).get_prob() * (*itr_a2).get_prob();
				}
			}
		}
	}


}

std::vector<Alignment> 
RNAduplexHommodel::
align_v(const TH& seq, double min_aln){

}

VVF
RNAduplexHommodel::
aln_duplex(Alignment a1, Alignment a2, VVF& hp) {
	// アラインメントに含まれるタグが同じことを前提とする。
  // 進化で構造が保存されるという考え方を利用するには、この条項が必要。
  // forループを書き換える。

  // hpに平均値を入れていく。
  // hybridization probability
  std::vector<std::string> s1=a1.get_seq();
  std::vector<std::string> s2=a2.get_seq();

  uint size1=s1.size();
  uint size2=s2.size();
  std::vector<std::string>::iterator it1=s1.begin();
  std::vector<std::string>::iterator it2=s2.begin();
  VVVVF vhp(size1);
  int i=0;
  ///////////
  // gapの時の確率を見たい
  it1++;
  //it1++;
  it2++;
  //std::cout << "S1[2][15]: " << (*it1)[15] << std::cout;// 2-15の塩基を取り出すコード

  //////////


  /*** RNAアラインメントのギャップ'-'を考慮してメモリを損傷しないようなコードにしなければならない。***/

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
       uint ucount=0;
       for (it_vvhp=vhp.begin(); it_vvhp!=vhp.end();it_vvhp++)
       {
         for(it_vhp=(*it_vvhp).begin(); it_vhp!=(*it_vvhp).end(); it_vhp++)
           sum+=(*it_vhp)[i][j];
       }
	     hp[i][j]=sum / (double)(size1*size2);// 上から移動
     }
   }
 }
}
