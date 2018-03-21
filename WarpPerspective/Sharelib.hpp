/*****************************************************************
Name :
Date : 2018/03/15
By   : CharlotteHonG
Final: 2018/03/21
*****************************************************************/
#pragma once
#include "Raw2Img.hpp"

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
// 快速 線性插值 (不做任何檢查可能會超出邊界)
static void fast_Bilinear_rgb(unsigned char* p, 
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

/*
struct basic_ImgData {
	std::vector<unsigned char> raw_img;
	uint32_t width;
	uint32_t height;
	uint16_t bits;
};*/

/*
// 線性取值
inline static double atBilinear_rgb(const vector<unsigned char>& img, 
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
static void fast_Bilinear_rgb(unsigned char* p, 
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
}*/
