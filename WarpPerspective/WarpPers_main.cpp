/*****************************************************************
Name :
Date : 2018/03/15
By   : CharlotteHonG
Final: 2018/03/16
*****************************************************************/
#include <iostream>
#include <vector>
#include <string>
#include "Raw2Img.hpp"
#include <timer.hpp>
using namespace std;

#include "WarpPers.hpp"
#include "WarpCyli.hpp"

//================================================================
int main(int argc, char const *argv[]) {
	Timer t1;
	const vector<double> HomogMat{
		 0.708484   ,  0.00428145 , 245.901,
		-0.103356   ,  0.888676   , 31.6815,
		-0.000390072, -1.61619e-05, 1
	};
	double ft = 672.673, Ax=752-533, Ay=500-496;

	/* 透視投影實現縫合 */
	test_WarpPers_Stitch();

	/* 圓柱投影實現縫合 */
	test_WarpCyli_AlphaBlend();
	return 0;
}
//================================================================
