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



// 線性取值
static double atBilinear_rgb(const vector<unsigned char>& img, 
	size_t width, double y, double x, size_t rgb)
{
	// 獲取鄰點(不能用 1+)
	double x0 = floor(x);
	double x1 = ceil(x);
	double y0 = floor(y);
	double y1 = ceil(y);
	// 獲取比例(只能用 1-)
	double dx1 = x - x0;
	double dx2 = 1 - dx1;
	double dy1 = y - y0;
	double dy2 = 1 - dy1;
	// 獲取點
	const double&& A = (double)img[(y0*width+x0)*3 +rgb];
	const double&& B = (double)img[(y0*width+x1)*3 +rgb];
	const double&& C = (double)img[(y1*width+x0)*3 +rgb];
	const double&& D = (double)img[(y1*width+x1)*3 +rgb];
	// 乘出比例(要交叉)
	double AB = A*dx2 + B*dx1;
	double CD = C*dx2 + D*dx1;
	double X = AB*dy2 + CD*dy1;
	return X;
}

// 輸入 dst 座標, 反轉 scr 輸出.
static void WarpPerspective_CoorTranfer_Inve(const vector<double>& HomogMat, double& x, double& y) {
	const double* H = HomogMat.data();
	const double i=x, j=y;

	x = (H[2] - H[8]*i) * (H[4] - H[7]*j) - 
		(H[1] - H[7]*i) * (H[5] - H[8]*j);
	y = (H[0] - H[6]*i) * (H[5] - H[8]*j) - 
		(H[2] - H[8]*i) * (H[3] - H[6]*j);

	double z = (H[1] - H[7]*i) * (H[3] - H[6]*j) - 
		(H[0] - H[6]*i) * (H[4] - H[7]*j);

	x /= z;
	y /= z;
}
// 輸入 scr 座標, 轉換 dst 輸出.
static void WarpPerspective_CoorTranfer(const vector<double>& HomogMat, double& x, double& y) {
	const double* H = HomogMat.data();
	const double i=x, j=y;

	x = H[0]*i + H[1]*y +H[2];
	y = H[3]*i + H[4]*y +H[5];
	double z = H[6]*i + H[7]*y +H[8];

	x /= z;
	y /= z;

	//x=round(x);
	//y=round(y);
}

// 透視轉換角點 輸入(xy*4) 輸出(dx, dy, minx, miny, maxx, maxy)
static vector<double> WarpPerspective_Corner(
	const vector<double>& HomogMat, size_t srcW, size_t srcH)
{
	vector<double> cn = {
		0, 0,  (double)(srcW-1.0), 0,
		0, (double)(srcH-1.0),   ((double)srcW-1.0), ((double)srcH-1.0)
	};
	// 透視轉換
	for(size_t i = 0; i < 4; i++) {
		//cout << cn[i*2+0] << ", " << cn[i*2+1] << "----->";
		WarpPerspective_CoorTranfer(HomogMat, cn[i*2+0], cn[i*2+1]);
		cn[i*2+0] = round(cn[i*2+0]);
		cn[i*2+1] = round(cn[i*2+1]);
		//cout << cn[i*2+0] << ", " << cn[i*2+1] << endl;
	}
	int max, min;
	// 找 x 最大最小
	max=INT_MIN, min=INT_MAX;
	for(size_t i = 0; i < 4; i++) {
		if(cn[i*2] > max) max=cn[i*2+0];
		if(cn[i*2] < min) min=cn[i*2+0];
	} cn[0] = max-min, cn[2] = min, cn[4] = max;
	// 找 y 最大最小
	max=INT_MIN, min=INT_MAX;
	for(size_t i = 0; i < 4; i++) {
		if(cn[i*2 +1] > max) max=cn[i*2+1];
		if(cn[i*2 +1] < min) min=cn[i*2+1];
	} cn[1] = max-min, cn[3] = min, cn[5] = max;
	// 反代逼近尋找正確值
	double dstW, dstH, tarW=cn[4], tarH=cn[5];
	for(size_t i=0, bw=0, bh=0; i < 100; i++) {
		double w=tarW+i, h=tarH+i;
		WarpPerspective_CoorTranfer_Inve(HomogMat, w, h);
		if(w > srcW and bw == 0) {
			bw=1;
			dstW=cn[4]+i;
		}
		if(h > srcH and bh == 0) {
			bh=1;
			dstH=cn[5]+i;
		}
		if(bh==1 and bw==1) break;
	}
	//cout << "final=" << dstW << ", " << dstH << endl;
	cn[0] = dstW-cn[2], cn[1] = dstH-cn[3];

	return cn;
}


// 圖像透視轉換 ImgRaw_basic 物件
void WarpPerspective(const basic_ImgData &src, basic_ImgData &dst, 
	const vector<double> &H, bool clip=0)
{
	int srcW = src.width;
	int srcH = src.height;
	// 獲得轉換後最大邊點
	vector<double> cn = WarpPerspective_Corner(H, srcW, srcH);
	// 起點位置
	int miny=0, minx=0;
	if(clip==1) { miny=-cn[3], minx=-cn[2]; }
	// 終點位置
	int dstW = cn[0]+cn[2]+minx;
	int dstH = cn[1]+cn[3]+miny;
	// 初始化 dst
	dst.raw_img.resize(dstW * dstH * src.bits/8.0);
	dst.width = dstW;
	dst.height = dstH;
	dst.bits = src.bits;
	// 透視投影
	int j, i;
	double x, y; 
	
#pragma omp parallel for private(i, j, x, y)
	for (j = -miny; j < dstH-miny; ++j) {
		for (i = -minx; i < dstW-minx; ++i){
			x = i, y = j;
			WarpPerspective_CoorTranfer_Inve(H, x, y);
			if ((x <= (double)srcW-1.0 && x >= 0.0) and
				(y <= (double)srcH-1.0 && y >= 0.0))
			{
				dst.raw_img[((j+miny)*dstW + (i+minx))*3 + 0] = atBilinear_rgb(src.raw_img, srcW, y, x, 0);
				dst.raw_img[((j+miny)*dstW + (i+minx))*3 + 1] = atBilinear_rgb(src.raw_img, srcW, y, x, 1);
				dst.raw_img[((j+miny)*dstW + (i+minx))*3 + 2] = atBilinear_rgb(src.raw_img, srcW, y, x, 2);

			}
		}
	}
}
void test1(string name, const vector<double>& HomogMat) {
	Timer t1;

	basic_ImgData img1, img2;
	Raw2Img::read_bmp(img1.raw_img, name, &img1.width, &img1.height, &img1.bits);
	t1.start();
	WarpPerspective(img1, img2, HomogMat, 0);
	t1.print(" WarpPerspective");

	Raw2Img::raw2bmp("WarpPers1.bmp", img2.raw_img, img2.width, img2.height, img2.bits);
}

// 圖像透視轉換 Raw 物件
/*void WarpPerspective(const Raw &src, Raw &dst, 
	const vector<double> &H, bool clip=0)
{

	int srcW = src.getCol();
	int srcH = src.getRow();
	// 獲得轉換後最大邊點
	vector<double> cn = WarpPerspective_Corner(H, srcW, srcH);
	// 起點位置
	int miny=0, minx=0;
	if(clip==1) {miny=-cn[3], minx=-cn[2];}
	// 終點位置
	int dstW = cn[0]+cn[2]+minx;
	int dstH = cn[1]+cn[3]+miny;
	
	dst.resize(dstW, dstH);
	// 透視投影
	int j, i;
	double x, y;
#pragma omp parallel for private(i, j, x, y)
	for (j = -miny; j < dstH-miny; ++j) {
		for (i = -minx; i < dstW-minx; ++i){
			x = i, y = j;
			WarpCylindrical_CoorTranfer_Inve(H, x, y);
			if ((x <= (double)srcW-1.0 && x >= 0.0) and
				(y <= (double)srcH-1.0 && y >= 0.0))
			{
				dst.RGB[((j+miny)*dstW + (i+minx))*3 + 0] = atBilinear_rgb(src.RGB, srcW, y, x, 0);
				dst.RGB[((j+miny)*dstW + (i+minx))*3 + 1] = atBilinear_rgb(src.RGB, srcW, y, x, 1);
				dst.RGB[((j+miny)*dstW + (i+minx))*3 + 2] = atBilinear_rgb(src.RGB, srcW, y, x, 2);

			}
		}
	}
}
void test2(string name, const vector<double>& HomogMat) {
	Timer t1;

	vector<unsigned char> raw_img;
	uint32_t weidth, heigh;
	uint16_t bits;
	Raw2Img::read_bmp(raw_img, name, &weidth, &heigh, &bits);

	Raw src(weidth, heigh), dst;
	src.RGB=raw_img;

	t1.start();
	WarpPerspective(src, dst, HomogMat, 1);
	t1.print(" WarpPerspective");
	Raw2Img::raw2bmp("WarpPers2.bmp", dst.RGB, dst.getCol(), dst.getRow());
}*/

// 圖像透視轉換 non 物件
void WarpPerspective(
	const vector<unsigned char> &src, const uint32_t srcW, const uint32_t srcH,
	vector<unsigned char> &dst, uint32_t& dstW, uint32_t& dstH,
	const vector<double> &H, bool clip=0)
{
	// 獲得轉換後最大邊點
	vector<double> cn = WarpPerspective_Corner(H, srcW, srcH);
	// 起點位置
	int miny=0, minx=0;
	if(clip==1) { miny=-cn[3], minx=-cn[2]; }
	// 終點位置
	dstW = cn[0]+cn[2]+minx;
	dstH = cn[1]+cn[3]+miny;
	dst.resize(dstW*dstH*3);
	// 透視投影
	int j, i;
	double x, y;
#pragma omp parallel for private(i, j, x, y)
	for (j = -miny; j < dstH-miny; ++j) {
		for (i = -minx; i < dstW-minx; ++i){
			x = i, y = j;
			WarpPerspective_CoorTranfer_Inve(H, x, y);
			if ((x <= (double)srcW-1.0 && x >= 0.0) and
				(y <= (double)srcH-1.0 && y >= 0.0))
			{
				dst[((j+miny)*dstW + (i+minx))*3 + 0] = atBilinear_rgb(src, srcW, y, x, 0);
				dst[((j+miny)*dstW + (i+minx))*3 + 1] = atBilinear_rgb(src, srcW, y, x, 1);
				dst[((j+miny)*dstW + (i+minx))*3 + 2] = atBilinear_rgb(src, srcW, y, x, 2);

			}
		}
	}
}
void test3(string name, const vector<double>& HomogMat) {
	Timer t1;

	vector<unsigned char> src;
	uint32_t srcW, srcH;
	uint16_t srcBits;
	Raw2Img::read_bmp(src, name, &srcW, &srcH, &srcBits);

	vector<unsigned char> dst;
	uint32_t dstW, dstH;
	uint16_t dstBits=srcBits;

	t1.start();
	WarpPerspective(src, srcW, srcH, dst, dstW, dstH, HomogMat, 1);
	t1.print(" WarpPerspective");
	Raw2Img::raw2bmp("WarpPers3.bmp", dst, dstW, dstH, dstBits);
}


// 透視轉換縫合範例
void Alpha_Blend(basic_ImgData& matchImg, 
	const basic_ImgData& imgL, const basic_ImgData& imgR, 
	const vector<double>& HomogMat)
{
	basic_ImgData warpImg;
	// 透視轉換
	WarpPerspective(imgR, warpImg, HomogMat, 0);
	// R 圖先補上去
	matchImg=warpImg;
	// 比例混合
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
				matchImg.raw_img[j*matchImg.width*3 + i*3+0] = 
					imgL.raw_img[j*imgL.width*3 + i*3+0];
				matchImg.raw_img[j*matchImg.width*3 + i*3+1] = 
					imgL.raw_img[j*imgL.width*3 + i*3+1];
				matchImg.raw_img[j*matchImg.width*3 + i*3+2] = 
					imgL.raw_img[j*imgL.width*3 + i*3+2];
			} else {
				// 這裡是重疊處.
				if(start==end) {
					start=i; // 紀錄起頭
				}
				if(start<end) {
					float len = end-start;
					float ratioR = (i-start)/len;
					float ratioL = 1.0 - ratioR;
					matchImg.raw_img[j*matchImg.width*3 + i*3+0] = 
						imgL.raw_img[j*imgL.width*3 + i*3+0]*ratioL + 
						warpImg.raw_img[j*warpImg.width*3 + i*3+0]*ratioR;
					matchImg.raw_img[j*matchImg.width*3 + i*3+1] = 
						imgL.raw_img[j*imgL.width*3 + i*3+1]*ratioL + 
						warpImg.raw_img[j*warpImg.width*3 + i*3+1]*ratioR;
					matchImg.raw_img[j*matchImg.width*3 + i*3+2] = 
						imgL.raw_img[j*imgL.width*3 + i*3+2]*ratioL + 
						warpImg.raw_img[j*warpImg.width*3 + i*3+2]*ratioR;
				}
			}
		}
	}

}
void test_WarpPers_Stitch() {
	Timer t1;
	// 透視矩陣
	const vector<double> HomogMat{
		0.708484   ,  0.00428145 , 245.901,
		-0.103356   ,  0.888676   , 31.6815,
		-0.000390072, -1.61619e-05, 1
	};
	// 讀取影像
	basic_ImgData img1;
	Raw2Img::read_bmp(img1.raw_img, "sc02.bmp", &img1.width, &img1.height, &img1.bits);
	basic_ImgData img2;
	Raw2Img::read_bmp(img2.raw_img, "sc03.bmp", &img2.width, &img2.height, &img2.bits);
	// 縫合影像
	basic_ImgData matchImg;
	t1.start();
	Alpha_Blend(matchImg, img1, img2, HomogMat);
	t1.print(" Alpha_Blend");
	Raw2Img::raw2bmp("Alpha_Blend.bmp", matchImg.raw_img, matchImg.width, matchImg.height, matchImg.bits);
}