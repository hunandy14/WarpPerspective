// http://blog.csdn.net/weixinhum/article/details/50611750

#include <iostream>
#include <vector>
#include <algorithm>
#include <timer.hpp>
using namespace std;

#include "imagedata.hpp"
#include "WarpCyli.hpp"

#define M_PI 3.14159265358979323846

// 線性取值
using at_T=float;
//template<class at_T>
static at_T atBilinear_rgb(const vector<unsigned char>& img, 
	size_t width, at_T y, at_T x, size_t rgb)
{
	// 獲取鄰點(不能用 1+)
	at_T x0 = floor(x);
	at_T x1 = ceil(x);
	at_T y0 = floor(y);
	at_T y1 = ceil(y);
	// 獲取比例(只能用 1-)
	at_T dx1 = x - x0;
	at_T dx2 = 1 - dx1;
	at_T dy1 = y - y0;
	at_T dy2 = 1 - dy1;
	// 獲取點
	const at_T A = img[(y0*width + x0)*3 +rgb];
	const at_T B = img[(y0*width + x1)*3 +rgb];
	const at_T C = img[(y1*width + x0)*3 +rgb];
	const at_T D = img[(y1*width + x1)*3 +rgb];
	// 乘出比例(要交叉)
	at_T X = 0;
	X += A*dx2*dy2;
	X += B*dx1*dy2;
	X += C*dx2*dy1;
	X += D*dx1*dy1;
	return X;
}

using lineType=double;
static void Bilinear_rgb(unsigned char* p, const vector<unsigned char>& img, 
	size_t width, lineType y, lineType x)
{
	// 獲取鄰點(不能用 1+)
	lineType x0 = floor(x);
	lineType x1 = ceil(x);
	lineType y0 = floor(y);
	lineType y1 = ceil(y);
	// 獲取比例(只能用 1-)
	lineType dx1 = x - x0;
	lineType dx2 = 1 - dx1;
	lineType dy1 = y - y0;
	lineType dy2 = 1 - dy1;
	// 獲取點
	lineType A, B, C, D;
	int rgb;
	
	rgb=0;
	A = img[(y0*width + x0)*3 +rgb];
	B = img[(y0*width + x1)*3 +rgb];
	C = img[(y1*width + x0)*3 +rgb];
	D = img[(y1*width + x1)*3 +rgb];
	*(p+rgb) = (unsigned char)(A*dx2*dy2 + B*dx1*dy2 + C*dx2*dy1 + D*dx1*dy1);

	rgb=1;
	A = img[(y0*width + x0)*3 +rgb];
	B = img[(y0*width + x1)*3 +rgb];
	C = img[(y1*width + x0)*3 +rgb];
	D = img[(y1*width + x1)*3 +rgb];
	*(p+rgb) = (unsigned char)(A*dx2*dy2 + B*dx1*dy2 + C*dx2*dy1 + D*dx1*dy1);

	rgb=2;
	A = img[(y0*width + x0)*3 +rgb];
	B = img[(y0*width + x1)*3 +rgb];
	C = img[(y1*width + x0)*3 +rgb];
	D = img[(y1*width + x1)*3 +rgb];
	*(p+rgb) = (unsigned char)(A*dx2*dy2 + B*dx1*dy2 + C*dx2*dy1 + D*dx1*dy1);
}


// bilinear補點
struct Color {
	unsigned char R;
	unsigned char G;
	unsigned char B;
};

static Color bilinear(const Raw &src, float _x, float _y)
{
	Color color;
	int x, y;
	x = (int)_x;
	y = (int)_y;

	float l_x = 0.f, r_x = 0.f;
	float t_y = 0.f, b_y = 0.f;

	l_x = _x - (float)x;
	r_x = 1.f - l_x;

	t_y = _y - (float)y;
	b_y = 1.f - t_y;

	float R = 0.f, G = 0.f, B = 0.f;
	int x1 = (x + 1) > (src.getCol() - 1) ? (x + 1) : (src.getCol() - 1);
	int y1 = (y + 1) > (src.getRow() - 1) ? (y + 1) : (src.getRow() - 1);

	if (x >= 0 && x < src.getCol() - 1 && y >= 0 && y < src.getRow() - 1)
	{
		R = (float)src.RGB[((y)* src.getCol() + (x)) * 3 + 0] * (r_x * b_y);
		G = (float)src.RGB[((y)* src.getCol() + (x)) * 3 + 1] * (r_x * b_y);
		B = (float)src.RGB[((y)* src.getCol() + (x)) * 3 + 2] * (r_x * b_y);

		R += (float)src.RGB[((y)* src.getCol() + (x + 1)) * 3 + 0] * (l_x * b_y);
		G += (float)src.RGB[((y)* src.getCol() + (x + 1)) * 3 + 1] * (l_x * b_y);
		B += (float)src.RGB[((y)* src.getCol() + (x + 1)) * 3 + 2] * (l_x * b_y);

		R += (float)src.RGB[((y + 1) * src.getCol() + (x)) * 3 + 0] * (r_x * t_y);
		G += (float)src.RGB[((y + 1) * src.getCol() + (x)) * 3 + 1] * (r_x * t_y);
		B += (float)src.RGB[((y + 1) * src.getCol() + (x)) * 3 + 2] * (r_x * t_y);

		R += (float)src.RGB[((y + 1) * src.getCol() + (x + 1)) * 3 + 0] * (l_x * t_y);
		G += (float)src.RGB[((y + 1) * src.getCol() + (x + 1)) * 3 + 1] * (l_x * t_y);
		B += (float)src.RGB[((y + 1) * src.getCol() + (x + 1)) * 3 + 2] * (l_x * t_y);
	}

	color.R = (unsigned char)(R > 255.0 ? 255 : R < 0.0 ? 0 : R);
	color.G = (unsigned char)(G > 255.0 ? 255 : G < 0.0 ? 0 : G);
	color.B = (unsigned char)(B > 255.0 ? 255 : B < 0.0 ? 0 : B);
	return color;
}
/************************ Functions prototyped here **************************/

// 圓柱投影座標反轉換
inline static void WarpPerspective_CoorTranfer_Inve(
	double R, size_t width, size_t height, 
	double& x, double& y)
{
	double r2 = (x - width*.5);
	double k = R / sqrt(R*R + r2*r2);
	x = (x - width *.5)/k + width *.5;
	y = (y - height*.5)/k + height*.5;
}

// 圓柱投影
Raw DealWithImgData(Raw &srcdata, int width, int height, double R)
{
	Raw drcdata(width, height * 2);
	
	int j, i;
	double* k_num = new double[width];

	// 圓柱投影
#pragma omp parallel for private(i, j)
	for (j = -(height / 2); j < (height * 3 / 2); j++){
		unsigned char* drcdata_RGB = &drcdata.RGB[(j + (height / 2)) * width * 3];
		for (i = 0; i < width; i++){
			double x = i, y = j;
			WarpPerspective_CoorTranfer_Inve(R, width, height, x, y);

			if (x >= 0 && y >= 0 && x < width - 1 && y < height - 1)
			{
				Color color = bilinear(srcdata, x, y);
				drcdata_RGB[i * 3 + 0] = color.R;
				drcdata_RGB[i * 3 + 1] = color.G;
				drcdata_RGB[i * 3 + 2] = color.B;
				// 不懂為什麼比較慢
				/*unsigned char* p = &drcdata.RGB[((j+(height / 2))*width + (i))*3 + 0];
				Bilinear_rgb(p, srcdata.RGB, width, y, x);*/
			}
		}
	}
	return drcdata;
}

void warping(const vector<Raw> &inputArrays, float FL2, 
	vector<Raw> &Output, vector<fpoint> &upedge, vector<fpoint> &downedge)
{
	float FL = FL2;
	Output.resize(inputArrays.size());

	Timer t;
	for (int idx = 0; idx < inputArrays.size(); idx++)
	{
		const Raw& image = inputArrays[idx];
		Output[idx].resize(image.getCol(), image.getRow());
		float mid_x = (float)image.getCol() / 2.f;
		float mid_y = (float)image.getRow() / 2.f;
		size_t curr_h = image.getRow();
		size_t curr_w = image.getCol();

		vector<float> k_core(curr_w);
		Timer t;
		for(int j = 0; j < curr_h; j++){
			for(int i = 0; i < curr_w; i++){
				if(j == 0)
					k_core[i] = sqrtf(FL * FL + (i - mid_x) * (i - mid_x)) / FL;
				float k = k_core[i];
				float x = (i-mid_x)*k + mid_x;
				float y = (j-mid_y)*k + mid_y;

				if(x >= 0 && y >= 0 &&
					x < curr_w && y < curr_h)
				{
					Color color = bilinear(image, x, y);
					Output[idx].RGB[(j*curr_w + i) * 3 + 0] = color.R;
					Output[idx].RGB[(j*curr_w + i) * 3 + 1] = color.G;
					Output[idx].RGB[(j*curr_w + i) * 3 + 2] = color.B;
				}
				if(j == 0) { // 圖上邊邊界.
					upedge.emplace_back(fpoint((float)x, (float)y));
				} else if(j == curr_h - 1) { // 圖下邊邊界.
					downedge.emplace_back(fpoint((float)x, (float)y));
				}
			}
		}
	}
	//t.print("warping time = ");
}




