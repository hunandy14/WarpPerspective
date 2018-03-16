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

//================================================================
int main(int argc, char const *argv[]) {
	Timer t1;
	const vector<double> HomogMat{
		 0.708484   ,  0.00428145 , 245.901,
		-0.103356   ,  0.888676   , 31.6815,
		-0.000390072, -1.61619e-05, 1
	};

	test1("kanna.bmp", HomogMat);
	test2("kanna.bmp", HomogMat);
	test3("kanna.bmp", HomogMat);
	return 0;
}
//================================================================
