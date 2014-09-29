#include "duplexHom.h"
#include "vvf.h"
#include "centroidalifold/folding_engine.h"
#include "centroidalifold/util.h"
#include "gtest/gtest.h"
#include "centroidalifold/probconsRNA/probcons.h"

#include "TestConfig.h"

#define HAVE_LIBRNA

#if 0
TEST(DuplexTest, SameSequence){
	PROBCONS::Probcons pc;
	std::vector<float> table;
	pc.ComputePosterior("AAA", "AAA", table, 0.001);
	// ComputePosteriorで出力されるtableの形式が知りたい
	std::vector<float> ans;
	ans.push_back(1);
	EXPECT_EQ(ans, table);
	std::cout << "size " << table.size() << std::endl;
}
# endif


#if 1
TEST(ProbabilityTest, AlignmentSameInput){
	float tmp_ans[2] = {1,0};
	VF tmp_ans2 (tmp_ans, tmp_ans + sizeof(tmp_ans) / sizeof(float));
	VVF mat;
	mat.push_back(tmp_ans2);
	tmp_ans[0] = 0;
	tmp_ans[1] = 1;
	tmp_ans2.clear();
	tmp_ans2.push_back(0);
	tmp_ans2.push_back(1);
	mat.push_back(tmp_ans2);

	RNAduplexHommodel ra;
	PROBCONS::Probcons pc;
	VVF res = ra.computeposterior(pc, "AAAAAA", "A");
	EXPECT_EQ(mat, res);
	EXPECT_EQ(6, res.size());
}
#endif

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}