/*****************************************************************
Name :
Date : 2018/03/15
By   : CharlotteHonG
Final: 2018/03/21
*****************************************************************/
// http://blog.csdn.net/weixinhum/article/details/50611750

#include <iostream>
#include <vector>
#include <algorithm>
#include <timer.hpp>
using namespace std;

#include "Raw2Img.hpp"
#include "WarpCyli.hpp"

#define M_PI 3.14159265358979323846

// 線性取值
inline static 
double atBilinear_rgb(const vector<unsigned char>& img, 
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
// 快速 線性插值 (不做任何檢查可能會超出邊界)
inline static 
void fast_Bilinear_rgb(unsigned char* p, 
	const basic_ImgData& src, double y, double x)
{
	// 起點
	int _x = (int)x;
	int _y = (int)y;
	// 左邊比值
	float l_x = x - (float)_x;
	float r_x = 1.f - l_x;
	float t_y = y - (float)_y;
	float b_y = 1.f - t_y;
	int srcW = src.width;
	int srcH = src.height;

	// 計算RGB
	float R = 0.f, G = 0.f, B = 0.f;
	R = (float)src.raw_img[((_y)* srcW + (_x)) * 3 + 0] * (r_x * b_y);
	G = (float)src.raw_img[((_y)* srcW + (_x)) * 3 + 1] * (r_x * b_y);
	B = (float)src.raw_img[((_y)* srcW + (_x)) * 3 + 2] * (r_x * b_y);

	R += (float)src.raw_img[((_y)* srcW + (_x + 1)) * 3 + 0] * (l_x * b_y);
	G += (float)src.raw_img[((_y)* srcW + (_x + 1)) * 3 + 1] * (l_x * b_y);
	B += (float)src.raw_img[((_y)* srcW + (_x + 1)) * 3 + 2] * (l_x * b_y);
						
	R += (float)src.raw_img[((_y + 1) * srcW + (_x)) * 3 + 0] * (r_x * t_y);
	G += (float)src.raw_img[((_y + 1) * srcW + (_x)) * 3 + 1] * (r_x * t_y);
	B += (float)src.raw_img[((_y + 1) * srcW + (_x)) * 3 + 2] * (r_x * t_y);
						
	R += (float)src.raw_img[((_y + 1) * srcW + (_x + 1)) * 3 + 0] * (l_x * t_y);
	G += (float)src.raw_img[((_y + 1) * srcW + (_x + 1)) * 3 + 1] * (l_x * t_y);
	B += (float)src.raw_img[((_y + 1) * srcW + (_x + 1)) * 3 + 2] * (l_x * t_y);

	*(p+0) = (unsigned char) R;
	*(p+1) = (unsigned char) G;
	*(p+2) = (unsigned char) B;
}
// 比例混合
inline static 
void AlphaBlend(basic_ImgData& matchImg, 
	const basic_ImgData& imgL, const basic_ImgData& imgR) {
	// R 圖先補上去
	matchImg=imgR;
	// 比例混合
	int i, j, start, end;
#pragma omp parallel for private(i, j, start, end)
	for(j = 0; j < imgL.height; j++) {
		start = imgL.width;
		end = imgL.width;
		for(i = 0; i <= (imgL.width-1); i++) {
			if( (imgR.raw_img[j*imgR.width*3 + i*3+0] == 0 and 
				imgR.raw_img[j*imgR.width*3 + i*3+1] == 0 and
				imgR.raw_img[j*imgR.width*3 + i*3+2] == 0)
				)
			{
				// 這裡要補原圖 L 的.
				matchImg.raw_img[j*matchImg.width*3 + i*3+0] = 
					imgL.raw_img[j*imgL.width*3 + i*3+0];
				matchImg.raw_img[j*matchImg.width*3 + i*3+1] = 
					imgL.raw_img[j*imgL.width*3 + i*3+1];
				matchImg.raw_img[j*matchImg.width*3 + i*3+2] = 
					imgL.raw_img[j*imgL.width*3 + i*3+2];
			} else {
				if(imgL.raw_img[j*imgL.width*3 + i*3+0] != 0 or
					imgL.raw_img[j*imgL.width*3 + i*3+1] != 0 or
					imgL.raw_img[j*imgL.width*3 + i*3+2] != 0)
				{
					// 這裡是重疊處.
					if(start==end) {
						start=i; // 紀錄起頭
					}
					if(start<end) {
						float len = end-start;
						float ratioR = (i-start)/len;
						float ratioL = 1.0 - ratioR;
						matchImg.raw_img[j*matchImg.width*3 + i*3+0] = //100;
							imgL.raw_img[j*imgL.width*3 + i*3+0]*ratioL + 
							imgR.raw_img[j*imgR.width*3 + i*3+0]*ratioR;
						matchImg.raw_img[j*matchImg.width*3 + i*3+1] = //0;
							imgL.raw_img[j*imgL.width*3 + i*3+1]*ratioL + 
							imgR.raw_img[j*imgR.width*3 + i*3+1]*ratioR;
						matchImg.raw_img[j*matchImg.width*3 + i*3+2] = //0;
							imgL.raw_img[j*imgL.width*3 + i*3+2]*ratioL + 
							imgR.raw_img[j*imgR.width*3 + i*3+2]*ratioR;
					}
				}
			}
		}
	}
}




// 圓柱投影座標反轉換
inline static  void WarpCylindrical_CoorTranfer_Inve(double R,
	size_t width, size_t height, 
	double& x, double& y)
{
	double r2 = (x - width*.5);
	double k = sqrt(R*R + r2*r2) / R;
	x = (x - width *.5)*k + width *.5;
	y = (y - height*.5)*k + height*.5;
}
// 圓柱投影 basic_ImgData
void WarpCylindrical(basic_ImgData &dst, const basic_ImgData &src, 
	double R ,int mx, int my, double edge)
{
	int w = src.width;
	int h = src.height;
	int moveH = (h*edge) + my;
	unsigned int moveW = mx;

	dst.raw_img.clear();
	dst.raw_img.resize((w+moveW)*3 *h *(1+edge*2));
	dst.width = w+moveW;
	dst.height = h * (1+edge*2);

	int j, i;
	// 圓柱投影
#pragma omp parallel for private(i, j)
	for (j = 0; j < h; j++){
		for (i = 0; i < w; i++){
			double x = i, y = j;
			WarpCylindrical_CoorTranfer_Inve(R, w, h, x, y);
			if (x >= 0 && y >= 0 && x < w - 1 && y < h - 1) {
				unsigned char* p = &dst.raw_img[((j+moveH)*(w+moveW) + (i+moveW)) *3];
				fast_Bilinear_rgb(p, src, y, x);
			}
		}
	}
}
// 圓柱投影縫合範例
void test_WarpCyli_AlphaBlend()
{
	double ft = 672.673, Ax=752-533, Ay=500-496;

	Timer t1;
	// 讀取影像
	basic_ImgData img1, dst1;
	Raw2Img::read_bmp(img1.raw_img, "sc02.bmp", &img1.width, &img1.height, &img1.bits);
	basic_ImgData img2, dst2;
	Raw2Img::read_bmp(img2.raw_img, "sc03.bmp", &img2.width, &img2.height, &img2.bits);

	t1.start();
	WarpCylindrical(dst1, img1, ft, 0, 0, 0.25);
	//Raw2Img::raw2bmp("WarpCyli1.bmp", dst1.raw_img, dst1.width, dst1.height);
	t1.print("WarpCylindrical_1");
	t1.start();
	WarpCylindrical(dst2, img2, ft, Ax, Ay, 0.25);
	//Raw2Img::raw2bmp("WarpCyli2.bmp", dst2.raw_img, dst2.width, dst2.height);
	t1.print("WarpCylindrical_2");

	basic_ImgData matchImg;

	t1.start();
	AlphaBlend(matchImg, dst1, dst2);
	t1.print("WarpCyli_AlphaBlend");
	
	Raw2Img::raw2bmp("WarpCyli_AlphaBlend.bmp", matchImg.raw_img, matchImg.width, matchImg.height);
}