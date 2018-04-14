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

#include <opencv2\opencv.hpp>
using namespace cv;

#include "Raw2Img.hpp"
#include "Sharelib.hpp"
#include "WarpCyli.hpp"
#include "Pyramid.hpp"
#include "GauBlur.hpp"


//==================================================================================
// 私有函式
//==================================================================================
// 重設 ImgData 大小
static void ImgData_resize(basic_ImgData &dst, int newW, int newH, int bits) {
	dst.raw_img.resize(newW*newH*3);
	dst.width = newW;
	dst.height = newH;
	dst.bits = bits;
};
static void ImgData_resize(const basic_ImgData& src, basic_ImgData &dst) {
	dst.raw_img.resize(src.width*src.height*3);
	dst.width = src.width;
	dst.height = src.height;
	dst.bits = src.bits;
};
// 輸出 bmp
static void ImgData_write(basic_ImgData &dst, string name) {
	Raw2Img::raw2bmp(name, dst.raw_img, dst.width, dst.height);
};
// 讀取bmp
static void ImgData_read(basic_ImgData &src, std::string name) {
	Raw2Img::read_bmp(src.raw_img, name, &src.width, &src.height, &src.bits);
}
//==================================================================================


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

	// 圓柱投影
#pragma omp parallel for
	for (int j = 0; j < h; j++){
		for (int i = 0; i < w; i++){
			double x = i, y = j;
			WarpCylindrical_CoorTranfer_Inve(R, w, h, x, y);
			if (x >= 0 && y >= 0 && x < w - 1 && y < h - 1) {
				unsigned char* p = &dst.raw_img[((j+moveH)*(w+moveW) + (i+moveW)) *3];
				fast_Bilinear_rgb(p, src, y, x);
			}
		}
	}
}
// 找到圓柱投影角點
void WarpCyliCorner(const basic_ImgData &src, vector<int>& corner) {
	corner.resize(4);
	// 左上角角點
	for (int i = 0; i < src.width; i++) {
		int pix = (int)src.raw_img[(src.height/2*src.width +i)*3 +0];
		if (i<src.width/2 and pix != 0) {
			corner[0]=i;
			//cout << "corner=" << corner[0] << endl;
			i=src.width/2;
		} else if (i>src.width/2 and pix == 0) {
			corner[2] = i-1;
			//cout << "corner=" << corner[2] << endl;
			break;
		}
	}
	// 右上角角點
	for (int i = 0; i < src.height; i++) {
		int pix = (int)src.raw_img[(i*src.width +corner[0])*3 +0];
		if (i<src.height/2 and pix != 0) {
			corner[1] = i;
			//cout << "corner=" << corner[2] << endl;
			i=src.height/2;
		} else if (i>src.height/2 and pix == 0) {
			corner[3] = i-1;
			//cout << "corner=" << corner[3] << endl;
			break;
		}
	}
}


// 取出重疊區
void getOverlap(const basic_ImgData &src1, const basic_ImgData &src2,
	basic_ImgData& cut1, basic_ImgData& cut2, vector<int> corner)
{
	// 偏移量
	int mx=corner[4];
	int my=corner[5];
	// 新圖大小
	int newH=corner[3]-corner[1]-my;
	int newW=corner[2]-corner[0]+mx;
	// 重疊區大小
	int lapH=newH;
	int lapW=corner[2]-corner[0]-mx;
	// 兩張圖的高度偏差值
	int myA = my>0? 0:my;
	int myB = my<0? 0:my;

	// 重疊區
	ImgData_resize(cut1, lapW, lapH, 24);
	ImgData_resize(cut2, lapW, lapH, 24);
#pragma omp parallel for
	for (int j = 0; j < newH; j++) {
		for (int i = 0; i < newW-mx; i++) {
			// 圖1
			if (i < corner[2]-corner[0]-mx) {
				for (int  rgb = 0; rgb < 3; rgb++) {
					cut1.raw_img[(j*cut1.width +i) *3+rgb] = src1.raw_img[(((j+myA)+corner[1])*src1.width +(i+corner[0]+mx)) *3+rgb];
				}
			}
			// 圖2
			if (i >= mx) {
				for (int  rgb = 0; rgb < 3; rgb++) {
					cut2.raw_img[(j*cut2.width +(i-mx)) *3+rgb] = src2.raw_img[(((j+myB)+corner[1])*src1.width +((i-mx)+corner[0])) *3+rgb];
				}
			}
		}
	}
	//write_img(cut1, "__cut1.bmp");
	//write_img(cut2, "__cut2.bmp");
}
// 重疊區與兩張原圖合併
void mergeOverlap(const basic_ImgData &src1, const basic_ImgData &src2,
	const basic_ImgData &blend, basic_ImgData &dst, vector<int> corner)
{
	// 偏移量
	int mx=corner[4];
	int my=corner[5];
	// 新圖大小
	int newH=corner[3]-corner[1]-my;
	int newW=corner[2]-corner[0]+mx;
	// 兩張圖的高度偏差值
	int myA = my>0? 0:my;
	int myB = my<0? 0:my;

	// 合併圖片
#pragma omp parallel for
	for (int j = 0; j < newH; j++) {
		for (int i = 0; i < newW; i++) {
			// 圖1
			if (i < mx) {
				for (int  rgb = 0; rgb < 3; rgb++) {
					dst.raw_img[(j*dst.width +i) *3+rgb] = src1.raw_img[(((j+myA)+corner[1])*src1.width +(i+corner[0])) *3+rgb];
				}
			}
			// 重疊區
			else if (i >= mx and i < corner[2]-corner[0]) {
				for (int  rgb = 0; rgb < 3; rgb++) {
					dst.raw_img[(j*dst.width +i) *3+rgb] = blend.raw_img[(j*blend.width+(i-mx)) *3+rgb];
				}
			}
			// 圖2
			else if (i >= corner[2]-corner[0]) {
				for (int  rgb = 0; rgb < 3; rgb++) {
					dst.raw_img[(j*dst.width +i) *3+rgb] = src2.raw_img[(((j+myB)+corner[1])*src1.width +((i-mx)+corner[0])) *3+rgb];
				}
			}
		}
	}
}

// 混合兩張圓柱
void WarpCyliMuitBlend(basic_ImgData &dst, const 
	basic_ImgData &src1, const basic_ImgData &src2,
	int mx, int my) 
{
	// 檢測圓柱圖角點(minX, minY, maxX, maxY, mx, my)
	vector<int> corner;
	WarpCyliCorner(src1, corner);
	corner.push_back(mx);
	corner.push_back(my);

	// 新圖大小
	int newH=corner[3]-corner[1]-my;
	int newW=corner[2]-corner[0]+mx;
	// 兩張圖的高度偏差值
	int myA = my>0? 0:my;
	int myB = my<0? 0:my;

	// 取出重疊區
	basic_ImgData cut1, cut2;
	getOverlap(src1, src2, cut1, cut2, corner);
	// 混合重疊區
	basic_ImgData blend;
	blendImg(blend, cut1, cut2);
	ImgData_write(blend, "___lapblend.bmp");

	// 測試:加入柏松混合圖
	basic_ImgData blend2;
	//Raw2Img::read_bmp(blend2.raw_img, "posiblend.bmp", &blend2.width, &blend2.height, &blend2.bits);
	//blend = blend2;

	// 合併三張圖片
	mergeOverlap(src1, src2, blend, dst, corner);
}



// 切割圓柱投影
void cutWarpCyliImg(
	const basic_ImgData &src1, const basic_ImgData &src2, 
	basic_ImgData &dst,
	const vector<int>& corner)
{
	// 偏移量
	int mx=corner[4];
	int my=corner[5];
	// 新圖大小
	int newH=corner[3]-corner[1]-my;
	int newW=corner[2]-corner[0]+mx;
	// 兩張圖的高度偏差值
	int myA = my>0? 0:my;
	int myB = my<0? 0:my;




}
// 柏松混合圓柱(測試中)
void WarpCyliMuitBlend_pos(basic_ImgData &dst, const 
	basic_ImgData &src1, const basic_ImgData &src2,
	int mx, int my) 
{
	// 檢測圓柱圖角點(minX, minY, maxX, maxY, mx, my)
	vector<int> corner;
	WarpCyliCorner(src1, corner);
	corner.push_back(mx);
	corner.push_back(my);
	// 新圖大小
	int newH=corner[3]-corner[1]-my;
	int newW=corner[2]-corner[0]+mx;
	// 兩張圖的高度偏差值
	int myA = my>0? 0:my;
	int myB = my<0? 0:my;


	//--------------------------------------------------
	// 整張圖(test用)
	ImgData_resize(dst, newW, newH, 24);
	basic_ImgData all=dst;

	// 像右延伸像素
	basic_ImgData right, right2;
	ImgData_resize(right, newH, 1, 24);

//#pragma omp parallel for
	for (int j = 0; j < newH; j++) {
		for (int i = 0; i < newW; i++) {
			int idx=(j*dst.width +i) *3;
			int src1idx;

			// 圖1
			if (i < /*corner[2]-corner[0]*/ newW/2.0) {
				src1idx=(((j+myA)+corner[1])*src1.width +(i+corner[0])) *3;
				for (int  rgb = 0; rgb < 3; rgb++) {
					all.raw_img[idx+rgb] = src1.raw_img[src1idx+rgb];
				}
			}
			// 向右拉平
			if (i >= /*mx*/ newW/2.0 ) {
				for (int  rgb = 0; rgb < 3; rgb++) {
					//all.raw_img[idx+rgb] = src1.raw_img[src1idx+rgb];
					right.raw_img[j*3+rgb] = src1.raw_img[src1idx+rgb];
				}
			}
			// 圖2
			if (i >= /*mx*/ newW/2.0) {
				for (int  rgb = 0; rgb < 3; rgb++) {
					//all.raw_img[idx+rgb] = src2.raw_img[(((j+myB)+corner[1])*src1.width +((i-mx)+corner[0])) *3+rgb];
				}
			}
		}
	}
	GauBlur(right, right2, 1.6, 9);
	for (int j = 0; j < newH; j++) {
		for (int i = 0; i < newW; i++) {
			int idx=(j*dst.width +i) *3;
			int src1idx;

			// 圖1
			if (i < /*corner[2]-corner[0]*/ newW/2.0) {
				src1idx=(((j+myA)+corner[1])*src1.width +(i+corner[0])) *3;
				for (int  rgb = 0; rgb < 3; rgb++) {
					//all.raw_img[idx+rgb] = src1.raw_img[src1idx+rgb];
				}
			}
			// 向右拉平
			if (i >= /*mx*/ newW/2.0 ) {
				for (int  rgb = 0; rgb < 3; rgb++) {
					all.raw_img[idx+rgb] = right2.raw_img[j*3+rgb];
				}
			}
			// 圖2
			if (i >= /*mx*/ newW/2.0) {
				for (int  rgb = 0; rgb < 3; rgb++) {
					//all.raw_img[idx+rgb] = src2.raw_img[(((j+myB)+corner[1])*src1.width +((i-mx)+corner[0])) *3+rgb];
				}
			}
		}
	}
	ImgData_write(all, "___all.bmp");
	//--------------------------------------------------

	// 取出重疊區
	basic_ImgData cut1, cut2;
	getOverlap(src1, src2, cut1, cut2, corner);
	// 混合重疊區
	basic_ImgData blend;
	blendImg(blend, cut1, cut2);
	ImgData_write(blend, "___lapblend.bmp");

	// 測試:加入柏松混合圖
	basic_ImgData blend2;
	//Raw2Img::read_bmp(blend2.raw_img, "posiblend.bmp", &blend2.width, &blend2.height, &blend2.bits);
	//blend = blend2;

	// 合併三張圖片
	mergeOverlap(src1, src2, blend, dst, corner);
}

// 圓柱投影縫合範例
void test_WarpCyli_AlphaBlend()
{
	double ft = 672.673, Ax=219, Ay=3;

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
void test_WarpCyli_MuitBlend()
{
	
	// sc
	//double ft = 672.673, Ax=219, Ay=3;
	//string name1 = "sc02.bmp", name2 = "sc03.bmp";
	// ball
	string name1 = "ball_01.bmp", name2 = "ball_02.bmp";
	double ft = 2252.97, Ax = 539, Ay = 37;

	Timer t1;
	// 讀取影像
	basic_ImgData img1, dst1;
	Raw2Img::read_bmp(img1.raw_img, name1, &img1.width, &img1.height, &img1.bits);
	basic_ImgData img2, dst2;  
	Raw2Img::read_bmp(img2.raw_img, name2, &img2.width, &img2.height, &img2.bits);

	// 投影圖1
	t1.start();
	WarpCylindrical(dst1, img1, ft);
	t1.print(" WarpCylindrical");
	//Raw2Img::raw2bmp("WarpCyliA.bmp", dst1.raw_img, dst1.width, dst1.height);

	// 投影圖2
	t1.start();
	WarpCylindrical(dst2, img2, ft);
	t1.print(" WarpCylindrical");
	//Raw2Img::raw2bmp("WarpCyliB.bmp", dst2.raw_img, dst2.width, dst2.height);

	// 縫合圖片.
	basic_ImgData matchImg;
	t1.start();
	WarpCyliMuitBlend_pos(matchImg, dst1, dst2, Ax, Ay);
	t1.print(" WarpCyliMuitBlend");
	Raw2Img::raw2bmp("WarpCyliMuitBlend.bmp", matchImg.raw_img, matchImg.width, matchImg.height);
}