/*****************************************************************
Name :
Date : 2018/03/15
By   : CharlotteHonG
Final: 2018/03/21
*****************************************************************/
#pragma once
#include "Raw2Img.hpp"

#define M_PI 3.14159265358979323846



// 快速 線性插值
inline static void fast_Bilinear_rgb(unsigned char* p, 
	const basic_ImgData& src, double y, double x)
{
	// 起點
	int _x = (int)x;
	int _y = (int)y;
	// 左邊比值
	double l_x = x - (double)_x;
	double r_x = 1.f - l_x;
	double t_y = y - (double)_y;
	double b_y = 1.f - t_y;
	int srcW = src.width;
	int srcH = src.height;

	// 計算RGB
	double R , G, B;
	int x2 = (_x+1) > src.width -1? src.width -1: _x+1;
	int y2 = (_y+1) > src.height-1? src.height-1: _y+1;
	R  = (double)src.raw_img[(_y * srcW + _x) *3 + 0] * (r_x * b_y);
	G  = (double)src.raw_img[(_y * srcW + _x) *3 + 1] * (r_x * b_y);
	B  = (double)src.raw_img[(_y * srcW + _x) *3 + 2] * (r_x * b_y);
	R += (double)src.raw_img[(_y * srcW + x2) *3 + 0] * (l_x * b_y);
	G += (double)src.raw_img[(_y * srcW + x2) *3 + 1] * (l_x * b_y);
	B += (double)src.raw_img[(_y * srcW + x2) *3 + 2] * (l_x * b_y);
	R += (double)src.raw_img[(y2 * srcW + _x) *3 + 0] * (r_x * t_y);
	G += (double)src.raw_img[(y2 * srcW + _x) *3 + 1] * (r_x * t_y);
	B += (double)src.raw_img[(y2 * srcW + _x) *3 + 2] * (r_x * t_y);
	R += (double)src.raw_img[(y2 * srcW + x2) *3 + 0] * (l_x * t_y);
	G += (double)src.raw_img[(y2 * srcW + x2) *3 + 1] * (l_x * t_y);
	B += (double)src.raw_img[(y2 * srcW + x2) *3 + 2] * (l_x * t_y);

	*(p+0) = (unsigned char) R;
	*(p+1) = (unsigned char) G;
	*(p+2) = (unsigned char) B;
}

// 快速補值
inline static void fast_NearestNeighbor_rgb(unsigned char* p,
	const basic_ImgData& src, double y, double x) 
{
	// 位置(四捨五入)
	int _x = (int)(x+0.5);
	int _y = (int)(y+0.5);
	int srcW = src.width;
	int srcH = src.height;

	// 計算RGB
	double R , G, B;
	int x2 = (_x+1) > src.width -1? src.width -1: _x+1;
	int y2 = (_y+1) > src.height-1? src.height-1: _y+1;
	R  = (double)src.raw_img[(y2 * srcW + x2) *3 + 0];
	G  = (double)src.raw_img[(y2 * srcW + x2) *3 + 1];
	B  = (double)src.raw_img[(y2 * srcW + x2) *3 + 2];

	*(p+0) = (unsigned char) R;
	*(p+1) = (unsigned char) G;
	*(p+2) = (unsigned char) B;
}

// 比例混合
inline static void AlphaBlend(basic_ImgData& matchImg, 
	const basic_ImgData& imgL, const basic_ImgData& imgR) {
	// R 圖先補上去
	matchImg=imgR;
	// 比例混合
	int i, j, start, end;
#pragma omp parallel for private(i, j, start, end)
	for(j = 0; j < (int)imgL.height; j++) {
		start = imgL.width;
		end = imgL.width;
		for(i = 0; i <= (int)(imgL.width-1); i++) {
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
						matchImg.raw_img[j*matchImg.width*3 + i*3+0] = (unsigned char)//100;
							(imgL.raw_img[j*imgL.width*3 + i*3+0]*ratioL + 
							imgR.raw_img[j*imgR.width*3 + i*3+0]*ratioR);
						matchImg.raw_img[j*matchImg.width*3 + i*3+1] = (unsigned char)//0;
							(imgL.raw_img[j*imgL.width*3 + i*3+1]*ratioL + 
							imgR.raw_img[j*imgR.width*3 + i*3+1]*ratioR);
						matchImg.raw_img[j*matchImg.width*3 + i*3+2] = (unsigned char)//0;
							(imgL.raw_img[j*imgL.width*3 + i*3+2]*ratioL + 
							imgR.raw_img[j*imgR.width*3 + i*3+2]*ratioR);
					}
				}
			}
		}
	}
}