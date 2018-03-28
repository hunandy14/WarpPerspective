/*****************************************************************
Name :
Date : 2018/03/28
By   : CharlotteHonG
Final: 2018/03/28
*****************************************************************/
#include <iostream>
#include <vector>
#include <string>
#include <timer.hpp>
using namespace std;

#include "Raw2Img.hpp"

#include <opencv2/opencv.hpp>
using namespace cv;

#include "Sharelib.hpp"
#include "WarpScale.hpp"
#include "GauBlur.hpp"
#include "Pyramid.hpp"


static void WarpScale(const basic_ImgData &src, basic_ImgData &dst, double Ratio){
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
			fast_NearestNeighbor_rgb(p, src, srcY, srcX);
		}
	}
}

void pyraUp(const basic_ImgData &src, basic_ImgData &dst) {
	int newH = (int)(src.height * 2.0);
	int newW = (int)(src.width  * 2.0);

	// 初始化 dst
	dst.raw_img.resize(newW * newH * src.bits/8.0);
	dst.width  = newW;
	dst.height = newH;
	dst.bits   = src.bits;

	basic_ImgData temp;
	WarpScale(src, temp, 2.0);
	GauBlur(temp, dst, 1.6, 4);
}

void pyraDown(const basic_ImgData &src, basic_ImgData &dst) {
	int newH = (int)(src.height * 0.5);
	int newW = (int)(src.width  * 0.5);

	// 初始化 dst
	dst.raw_img.clear();
	dst.raw_img.resize(newW * newH * src.bits/8.0);
	dst.width  = newW;
	dst.height = newH;
	dst.bits   = src.bits;

	basic_ImgData temp;
	WarpScale(src, temp, 0.5);
	GauBlur(temp, dst, 1.6, 4);
}

void buildPyramids(const basic_ImgData &src, vector<basic_ImgData> &pyr, int octvs=5) {
	pyr.clear();
	pyr.resize(octvs);
	pyr[0]=src;
	for(int i = 1; i < octvs; i++) {
		pyraDown(pyr[i-1], pyr[i]);
	}
}
void buildLaplacianPyramids(const basic_ImgData &src, vector<basic_ImgData> &pyr, int octvs=5) {
	pyr.clear();
	pyr.resize(octvs);
	pyr[0]=src;

	for(int i = 1; i < octvs; i++) {
		basic_ImgData expend;
		pyraDown(pyr[i-1], pyr[i]);
		WarpScale(pyr[i], expend, 2.0);
#pragma omp parallel for
		for(int k = 0; k < expend.width * expend.height*3; k++) {
			pyr[i-1].raw_img[k] = (pyr[i-1].raw_img[k] - expend.raw_img[k])+128;
		}
	}
}
void reLaplacianPyramids(vector<basic_ImgData> &pyr, basic_ImgData &dst, int octvs=5) {
	int newH = (int)(pyr[0].height);
	int newW = (int)(pyr[0].width);

	// 初始化 dst
	dst.raw_img.clear();
	dst.raw_img.resize(newW * newH * pyr[0].bits/8.0);
	dst.width  = newW;
	dst.height = newH;
	dst.bits   = pyr[0].bits;

	for(int i = octvs-1; i >= 1; i--) {
		basic_ImgData expend;
		WarpScale(pyr[i], expend, 2.0);
		//pyraUp(pyr[i], expend);
#pragma omp parallel for
		for(int k = 0; k < expend.width * expend.height*3; k++) {
			pyr[i-1].raw_img[k] = (pyr[i-1].raw_img[k] + (expend.raw_img[k])-128);
		}
	}
	dst = pyr[0];
}

void test_pyramids() {
	Timer t1;
	// 讀取影像
	basic_ImgData img1, img2;
	Raw2Img::read_bmp(img1.raw_img, "sc02.bmp", &img1.width, &img1.height, &img1.bits);

	// 拉普拉斯金字塔
	vector<basic_ImgData> lap;
	t1.start();
	buildLaplacianPyramids(img1, lap);
	t1.print(" buildLaplacianPyramids");
	for(int i = 0; i < lap.size(); i++) {
		basic_ImgData& img2=lap[i];
		Raw2Img::raw2bmp("Lap\\Lap"+ to_string(i) +".bmp", img2.raw_img, img2.width, img2.height, img2.bits);
	}

	// 還原拉普拉斯金字塔
	basic_ImgData res;
	t1.start();
	reLaplacianPyramids(lap, res);
	t1.print(" rebuildLaplacianPyramids");
	Raw2Img::raw2bmp("Pyramids.bmp", res.raw_img, res.width, res.height, res.bits);
}