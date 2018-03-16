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

void WarpPers_Stitch(basic_ImgData& matchImg, const basic_ImgData& imgL, const basic_ImgData& imgR, 
	const vector<double>& HomogMat)
{
	basic_ImgData warpImg;
	WarpPerspective(imgR, warpImg, HomogMat, 0);

	// 比例混合
	matchImg=warpImg;
	int i, j, start, end;
#pragma omp parallel for private(i, j, start, end)
	for(j = 0; j < imgL.height; j++) {
		start = imgL.width;
		end = imgL.width;
		for(i = 0; i <= (imgL.width-1); i++) {
			if(warpImg.raw_img[j*warpImg.width*3 + i*3+0] == 0 and 
				warpImg.raw_img[j*warpImg.width*3 + i*3+1] == 0 and
				warpImg.raw_img[j*warpImg.width*3 + i*3+2] == 0)
			{
				// 這裡要補原圖 L 的.
				matchImg.raw_img[j*matchImg.width*3 + i*3+0] = imgL.raw_img[j*imgL.width*3 + i*3+0];
				matchImg.raw_img[j*matchImg.width*3 + i*3+1] = imgL.raw_img[j*imgL.width*3 + i*3+1];
				matchImg.raw_img[j*matchImg.width*3 + i*3+2] = imgL.raw_img[j*imgL.width*3 + i*3+2];
			} else {
				// 這裡是重疊處.
				// 比例.
				if(start==end) {
					start=i; // 紀錄起頭
				}
				if(start<end) {
					float len = end-start;
					float ratioR = (i-start)/len;
					float ratioL = 1.0 - ratioR;
					matchImg.raw_img[j*matchImg.width*3 + i*3+0] = 
						imgL.raw_img[j*imgL.width*3 + i*3+0]*ratioL + warpImg.raw_img[j*warpImg.width*3 + i*3+0]*ratioR;
					matchImg.raw_img[j*matchImg.width*3 + i*3+1] = 
						imgL.raw_img[j*imgL.width*3 + i*3+1]*ratioL + warpImg.raw_img[j*warpImg.width*3 + i*3+1]*ratioR;
					matchImg.raw_img[j*matchImg.width*3 + i*3+2] = 
						imgL.raw_img[j*imgL.width*3 + i*3+2]*ratioL + warpImg.raw_img[j*warpImg.width*3 + i*3+2]*ratioR;
				}
			}
		}
	}

}
//================================================================
int main(int argc, char const *argv[]) {
	Timer t1;
	const vector<double> HomogMat{
		 0.708484   ,  0.00428145 , 245.901,
		-0.103356   ,  0.888676   , 31.6815,
		-0.000390072, -1.61619e-05, 1
	};

	//test1("sc03.bmp", HomogMat);
	//test2("kanna.bmp", HomogMat);
	//test3("kanna.bmp", HomogMat);

	// 讀取影像
	basic_ImgData img1;
	Raw2Img::read_bmp(img1.raw_img, "sc02.bmp", &img1.width, &img1.height, &img1.bits);
	basic_ImgData img2;
	Raw2Img::read_bmp(img2.raw_img, "sc03.bmp", &img2.width, &img2.height, &img2.bits);
	// 縫合影像
	basic_ImgData matchImg;
	t1.start();
	WarpPers_Stitch(matchImg, img1, img2, HomogMat);
	t1.print(" WarpPers_Stitch");
	Raw2Img::raw2bmp("WarpPers_Stitch.bmp", matchImg.raw_img, matchImg.width, matchImg.height, matchImg.bits);
	return 0;
}
//================================================================
