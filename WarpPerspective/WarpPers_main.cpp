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

#include "imagedata.hpp"
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

	/* 透視投影 */
	//test1("sc03.bmp", HomogMat);
	//test2("kanna.bmp", HomogMat);
	//test3("kanna.bmp", HomogMat);

	/* 透視投影實現縫合 */
	//test_WarpPers_Stitch();

	/* 圓柱投影 */
	// 讀取影像
	basic_ImgData img1;
	Raw2Img::read_bmp(img1.raw_img, "sc02.bmp", &img1.width, &img1.height, &img1.bits);
	basic_ImgData img2, dst2, zero;
	Raw2Img::read_bmp(img2.raw_img, "sc03.bmp", &img2.width, &img2.height, &img2.bits);

	string name = "sc02.bmp";
	vector<unsigned char> raw_img;
	uint32_t weidth, heigh;
	uint16_t bits;
	Raw2Img::read_bmp(raw_img, name, &weidth, &heigh, &bits);

	Raw src(weidth, heigh), dst;
	src.RGB=raw_img;


	t1.start();
//	DealWithImgData(dst, src, ft);
	WarpCylindrical(dst2, img1, ft, 0, 0, 0.25);
	Raw2Img::raw2bmp("WarpCyli1.bmp", dst2.raw_img, dst2.width, dst2.height);

	WarpCylindrical(dst2, img2, ft, Ax, Ay, 0.25);
	t1.print("DealWithImgData");

	//Raw2Img::raw2bmp("WarpCyli.bmp", dst.RGB, dst.getCol(), dst.getRow());
	Raw2Img::raw2bmp("WarpCyli2.bmp", dst2.raw_img, dst2.width, dst2.height);



	return 0;
}
//================================================================
