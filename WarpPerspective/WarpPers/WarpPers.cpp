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
static void WarpPerspective_Corner(const vector<double>& HomogMat, vector<double>& cn) {
	int srcW = cn[6], srcH=cn[7];
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
}

// 圖像透視轉換
void WarpPerspective(const ImgRaw_basic &src, ImgRaw_basic &dst, 
	const vector<double> &H, bool clip=0)
{
	int miny=0, minx=0;

	// 輸入座標
	vector<double> cn={
		0, 0,  (double)(src.width-1.0), 0,
		0, (double)(src.height-1.0),   ((double)src.width-1.0), ((double)src.height-1.0)
	};
	// 獲得轉換後最大邊點
	WarpPerspective_Corner(H, cn);
	if(clip==1) miny=-cn[3], minx=-cn[2];
	int newW= cn[0]+cn[2], newH=cn[1]+cn[3];
	dst.resize(newW+minx, newH+miny, src.bits);
	// 透視投影
	int j, i;
	double x, y;
	int dstW = dst.width;
	int dstH = dst.height;
	int srcW = src.width;
	int srcH = src.height;
#pragma omp parallel for private(i, j, x, y)
	for (j = -miny; j < dstH-miny; ++j) {
		for (i = -minx; i < dstW-minx; ++i){
			x = i, y = j;
			WarpPerspective_CoorTranfer_Inve(H, x, y);
			if ((x <= (double)srcW-1.0 && x >= 0.0) and
				(y <= (double)srcH-1.0 && y >= 0.0))
			{
				dst[((j+miny)*dstW + (i+minx))*3 + 0] = atBilinear_rgb(src.raw_img, srcW, y, x, 0);
				dst[((j+miny)*dstW + (i+minx))*3 + 1] = atBilinear_rgb(src.raw_img, srcW, y, x, 1);
				dst[((j+miny)*dstW + (i+minx))*3 + 2] = atBilinear_rgb(src.raw_img, srcW, y, x, 2);

			}
		}
	}
}
void test1(string name, const vector<double>& HomogMat) {
	Timer t1;

	ImgRaw_basic img1, img2;
	Raw2Img::read_bmp(img1, name, &img1.width, &img1.height, &img1.bits);
	t1.start();
	WarpPerspective(img1, img2, HomogMat, 1);
	t1.print(" WarpPerspective");

	Raw2Img::raw2bmp("WarpPers1.bmp", img2, img2.width, img2.height, img2.bits);
}

void WarpPerspective(const Raw &src, Raw &dst, 
	const vector<double> &H, bool clip=0)
{
	int miny=0, minx=0;

	// 輸入座標
	vector<double> cn={
		0, 0,  (double)(src.getCol()-1.0), 0,
		0, (double)(src.getRow()-1.0),   ((double)src.getCol()-1.0), ((double)src.getRow()-1.0)
	};
	// 獲得轉換後最大邊點
	WarpPerspective_Corner(H, cn);
	if(clip==1) miny=-cn[3], minx=-cn[2];
	int newW= cn[0]+cn[2], newH=cn[1]+cn[3];
	dst.resize(newW+minx, newH+miny);
	// 透視投影
	int j, i;
	double x, y;
	int dstW = dst.getCol();
	int dstH = dst.getRow();
	int srcW = src.getCol();
	int srcH = src.getRow();
#pragma omp parallel for private(i, j, x, y)
	for (j = -miny; j < dstH-miny; ++j) {
		for (i = -minx; i < dstW-minx; ++i){
			x = i, y = j;
			WarpPerspective_CoorTranfer_Inve(H, x, y);
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
}