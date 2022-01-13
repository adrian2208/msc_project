#pragma once
#include "SU3_mat.h"
#include <vector>
class SU3_gen{
public:
	SU3_gen();
	su3_mat& operator() (int i);
private:
	std::vector<su3_mat> T;
};
