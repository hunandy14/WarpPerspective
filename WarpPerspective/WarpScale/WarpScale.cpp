/*****************************************************************
Name :
Date : 2018/03/15
By   : CharlotteHonG
Final: 2018/03/16
*****************************************************************/
#include <iostream>
#include <vector>
#include <string>
#include <timer.hpp>
using namespace std;

#include "Raw2Img.hpp"
#include "Sharelib.hpp"
#include "WarpScale.hpp"

void WarpScale(const basic_ImgData &src, basic_ImgData &dst, double Ratio){
	int newH = (int)(src.height * Ratio);
	int newW = (int)(src.width  * Ratio);
	// 初始化 dst
	dst.raw_img.resize(newW * newH * src.bits/8.0);
	dst.width  = newW;
	dst.height = newH;
	dst.bits   = src.bits;
	// 跑新圖座標

	int i, j;
#pragma omp parallel for private(i, j)
	for (j = 0; j < newH; ++j) {
		for (i = 0; i < newW; ++i) {
			// 調整對齊
			double srcY, srcX;
			if (Ratio < 1) {
				srcY = ((j+0.5f)/Ratio) - 0.5;
				srcX = ((i+0.5f)/Ratio) - 0.5;
			} else {
				srcY = j*(src.height-1.f) / (newH-1.f);
				srcX = i*(src.width -1.f) / (newW-1.f);
			}
			// 獲取插補值
			unsigned char* p = &dst.raw_img[(j*newW + i) *3];
			fast_Bilinear_rgb(p, src, srcY, srcX);
		}
	}
}

void test_WarpScale() {
	Timer t1;
	// 讀取影像
	basic_ImgData img1, img2;
	Raw2Img::read_bmp(img1.raw_img, "sc02.bmp", &img1.width, &img1.height, &img1.bits);
	
	t1.start();
	WarpScale(img1, img2, 2);
	t1.print(" WarpScale");

	Raw2Img::raw2bmp("WarpScale.bmp", img2.raw_img, img2.width, img2.height, img2.bits);
}