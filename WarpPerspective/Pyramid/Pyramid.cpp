/*****************************************************************
Name :
Date : 2018/03/28
By   : CharlotteHonG
Final: 2018/03/28

## 遇到的問題
1. 金字要用到縮放，如果是非偶數會導致一縮一放之後不一樣
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

// 高斯2d矩陣
static vector<float> gau_matrix2d(float p, size_t mat_len) {
	vector<float> gau_mat2d;
	gau_mat2d.resize(mat_len*mat_len);
	vector<float> gau_mat1d = gau_matrix(p, mat_len);
	cout << gau_mat1d.size() << endl;


	// 做 X
	int i=0, j=0;
	for(size_t j = 0; j < mat_len; j++){
		for(size_t i = 0; i < mat_len; i++){
			gau_mat2d[j*mat_len + i] = gau_mat1d[i];
		}
	}
	// 做 Y
	for(size_t j = 0; j < mat_len; j++){
		for(size_t i = 0; i < mat_len; i++){
			gau_mat2d[i*mat_len + j] *= gau_mat1d[i];
		}
	}
	//gau_mat2d[a*mat_len + j] *= gau_mat1d[a];
	return gau_mat2d;
}
static vector<float> gau_matrix2d(float p, size_t w, size_t h) {
	vector<float> gau_mat2d(w*h);
	vector<float> gau_matW = gau_matrix(p, w);
	vector<float> gau_matH = gau_matrix(p, h);

	// 做 X
	int i, j;
//#pragma omp parallel for private(a, j)
	for(j = 0; j < h; j++){
		for(i = 0; i < w; i++){
			gau_mat2d[j*w + i] = gau_matW[i];
		}
	}
	// 做 Y
//#pragma omp parallel for private(a, j)
	for(j = 0; j < w; j++){
		for(i = 0; i < h; i++){
			gau_mat2d[i*w + j] *= gau_matH[i];
		}
	}
	return gau_mat2d;
}
static vector<float> getGauKer(int x){
	vector<float> kernel(x);
	double half = (x-1) / 2.f;

	constexpr double rlog5_2 = -0.721348; // 1 / (2.f*log(0.5f))
	double sigma = sqrt( -powf(x-1-half, 2.f) * rlog5_2 );
	double rSigma22 = 1.0/(2 * sigma * sigma);

//#pragma omp parallel for
	for(int i = 0; i < x; i++){
		float g;
		if(i <= (x - half)){
			g = exp( -(i*i*rSigma22) );
		} else{
			g = 1.0 - exp(-powf(x-i-1, 2.f) * rSigma22);
		}
		kernel[i] = g;
	}
	return kernel;
}

// 圖片縮放
static void WarpScale(const basic_ImgData &src, basic_ImgData &dst, double Ratio){
	int newH = (int)((src.height * Ratio) +0.5);
	int newW = (int)((src.width  * Ratio) +0.5);
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
			//fast_NearestNeighbor_rgb(p, src, srcY, srcX);
			fast_Bilinear_rgb(p, src, srcY, srcX);
		}
	}
}

// 金字塔層處理
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
void imgSub(basic_ImgData &src, const basic_ImgData &dst) {
	int i, j;
#pragma omp parallel for private(i, j)
	for (j = 0; j < src.height; j++) {
		for (i = 0; i < src.width; i++) {
			int srcIdx = (j*src.width + i) * 3;
			int dstIdx = (j*dst.width + i) * 3;

			int pixR = (int)src.raw_img[srcIdx+0] - (int)dst.raw_img[dstIdx+0] +128;
			int pixG = (int)src.raw_img[srcIdx+1] - (int)dst.raw_img[dstIdx+1] +128;
			int pixB = (int)src.raw_img[srcIdx+2] - (int)dst.raw_img[dstIdx+2] +128;

			pixR = pixR <0? 0: pixR;
			pixG = pixG <0? 0: pixG;
			pixB = pixB <0? 0: pixB;
			pixR = pixR >255? 255: pixR;
			pixG = pixG >255? 255: pixG;
			pixB = pixB >255? 255: pixB;

			src.raw_img[srcIdx+0] = pixR;
			src.raw_img[srcIdx+1] = pixG;
			src.raw_img[srcIdx+2] = pixB;
		}
	}
}
void imgAdd(basic_ImgData &src, const basic_ImgData &dst) {
	int i, j;
#pragma omp parallel for private(i, j)
	for (j = 0; j < src.height; j++) {
		for (i = 0; i < src.width; i++) {
			int srcIdx = (j*src.width + i) * 3;
			int dstIdx = (j*dst.width + i) * 3;

			int pixR = (int)src.raw_img[srcIdx+0] + (int)dst.raw_img[dstIdx+0] -128;
			int pixG = (int)src.raw_img[srcIdx+1] + (int)dst.raw_img[dstIdx+1] -128;
			int pixB = (int)src.raw_img[srcIdx+2] + (int)dst.raw_img[dstIdx+2] -128;

			pixR = pixR <0? 0: pixR;
			pixG = pixG <0? 0: pixG;
			pixB = pixB <0? 0: pixB;
			pixR = pixR >255? 255: pixR;
			pixG = pixG >255? 255: pixG;
			pixB = pixB >255? 255: pixB;

			src.raw_img[srcIdx+0] = pixR;
			src.raw_img[srcIdx+1] = pixG;
			src.raw_img[srcIdx+2] = pixB;
		}
	}
}

// 金字塔
using LapPyr = vector<basic_ImgData>;
void buildPyramids(const basic_ImgData &src, vector<basic_ImgData> &pyr, int octvs=5) {
	pyr.clear();
	pyr.resize(octvs);
	pyr[0]=src;
	for(int i = 1; i < octvs; i++) {
		pyraDown(pyr[i-1], pyr[i]);
	}
}
void buildLaplacianPyramids(const basic_ImgData &src, LapPyr &pyr, int octvs=5) {
	Timer t1;
	pyr.clear();
	pyr.resize(octvs);
	pyr[0]=src;

	for(int i = 1; i < octvs; i++) {
		basic_ImgData expend;
		pyraDown(pyr[i-1], pyr[i]);
		WarpScale(pyr[i], expend, 2.0);
		imgSub(pyr[i-1], expend);
		//basic_ImgData& dst_out = pyr[i - 1];
		//Raw2Img::raw2bmp("Lap\\Lap_Test"+ to_string(i) +".bmp", dst_out.raw_img, dst_out.width, dst_out.height, dst_out.bits);
	}
}
void reLaplacianPyramids(LapPyr &pyr, basic_ImgData &dst, int octvs=5) {
	Timer t1;
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
		imgAdd(pyr[i-1], expend);
	}
	dst = pyr[0];
}

// 混合拉普拉斯金字塔
void blendLaplacianPyramids(LapPyr& LS, const LapPyr& LA, const LapPyr& LB) {
	LS.resize(LA.size());
	
	Timer t1;
	// 高斯矩陣
	auto gausKernal = getGauKer(LA.back().width);

	// 混合圖片
	for(int a = 0; a < LS.size(); a++) {
		int newH = (int)(LA[a].height);
		int newW = (int)(LA[a].width);
		int center = (int)(LA[a].width *0.5);

		// 初始化
		basic_ImgData dst;
		dst.raw_img.resize(newW * newH * LA[a].bits);
		dst.width  = newW;
		dst.height = newH;
		dst.bits   = LA[a].bits;

		// 開始跑 LS[a] 圖
		for(int j = 0; j < newH; j++) {
			for(int i = 0; i < newW; i++) {
				for(int c = 0; c < 3; c++) {
					// 拉普拉斯彩色區 (L*高斯) + (R*(1-高斯))
					if(a == LS.size()-1) {
						dst.raw_img[(j * dst.width + i) * 3 + c] = 
							LA[a].raw_img[(j * LA[a].width + i) * 3 + c] * gausKernal[i] +
							LB[a].raw_img[(j * LB[a].width + i) * 3 + c] * (1.f - gausKernal[i]);
						// 拉普拉斯差值區 (左邊就放左邊差值，右邊放右邊差值，正中間放平均)
					} else {
						if(i == center) {
							dst.raw_img[(j * dst.width + i) * 3 + c] = 
								(LA[a].raw_img[(j * LA[a].width + i) * 3 + c] + 
									LB[a].raw_img[(j * LB[a].width + i) * 3 + c]) * 0.5;
							// 右半部
						} else if(i > center) {
							dst.raw_img[(j * dst.width + i) * 3 + c] = 
								LB[a].raw_img[(j * LB[a].width + i) * 3 + c];
							// 左半部
						} else {
							dst.raw_img[(j * dst.width + i) * 3 + c] = 
								LA[a].raw_img[(j * LA[a].width + i) * 3 + c];
						}
					}
				}

			}
		}
		LS[a] = std::move(dst);
	}
}

// 混合圖片
void test_pyramids() {
	Timer t1;
	basic_ImgData src1, src2, dst;
	// 讀取影像
	//string name1="LA.bmp", name2="LB.bmp";
	string name1="white.bmp", name2="apple.bmp";
	Raw2Img::read_bmp(src1.raw_img, name1, &src1.width, &src1.height, &src1.bits);
	Raw2Img::read_bmp(src2.raw_img, name2, &src2.width, &src2.height, &src2.bits);
	
	// 拉普拉斯金字塔 AB
	vector<basic_ImgData> LA, LB;
	t1.start();
	buildLaplacianPyramids(src1, LA);
	t1.print(" buildLapA");
	t1.start();
	buildLaplacianPyramids(src2, LB);
	t1.print(" buildLapB");

	// 混合金字塔
	LapPyr LS;
	t1.start();
	blendLaplacianPyramids(LS, LA, LB);
	t1.print(" blendImg");
	
	// 還原拉普拉斯金字塔
	t1.start();
	reLaplacianPyramids(LS, dst);
	t1.print(" rebuildLaplacianPyramids");
	Raw2Img::raw2bmp("_rePyrImg.bmp", dst.raw_img, dst.width, dst.height, dst.bits);
	
}