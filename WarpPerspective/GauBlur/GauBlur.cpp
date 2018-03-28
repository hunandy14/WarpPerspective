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
#include "Sharelib.hpp"
#include "WarpScale.hpp"
#include "GauBlur.hpp"

// 來源相同的例外
class file_same : public std::runtime_error {
public:
	file_same(const string& str): std::runtime_error(str) {}
};
// 高斯公式
float gau_meth(size_t r, float p) {
	float two = 2.0;
	float num = exp(-pow(r, two) / (two*pow(p, two)));
	num /= sqrt(two*M_PI)*p;
	return num;
}
// 高斯矩陣 (mat_len defa=3)
vector<float> gau_matrix(float p, size_t mat_len) {
	vector<float> gau_mat;
	// 計算矩陣長度
	if (mat_len == 0) {
		//mat_len = (int)(((p - 0.8) / 0.3 + 1.0) * 2.0);// (顏瑞穎給的公式)
		mat_len = (int)(round((p*6 + 1))) | 1; // (opencv的公式)
	}
	// 奇數修正
	if (mat_len % 2 == 0) { ++mat_len; }
	// 一維高斯矩陣
	gau_mat.resize(mat_len);
	float sum = 0;
	for (int i = 0, j = mat_len / 2; j < mat_len; ++i, ++j) {
		float temp;
		if (i) {
			temp = gau_meth(i, p);
			gau_mat[j] = temp;
			gau_mat[mat_len - j - 1] = temp;
			sum += temp += temp;
		}
		else {
			gau_mat[j] = gau_meth(i, p);
			sum += gau_mat[j];
		}
	}
	// 歸一化
	for (auto&& i : gau_mat) { i /= sum; }
	return gau_mat;
}

// 高斯模糊
using types=float;
void GauBlur(const basic_ImgData& src, basic_ImgData& dst, float p, size_t mat_len)
{
	/*vector<unsigned char> raw_img = src.raw_img;
	Raw2Img::raw2gray(raw_img);*/

	size_t width  = src.width;
	size_t height = src.height;


	// 來源不可相同
	if (&dst == &src) {
		throw file_same("## Erroe! in and out is same.");
	}
	vector<types> gau_mat = gau_matrix(p, mat_len);

	// 初始化 dst
	dst.raw_img.resize(width*height * src.bits/8.0);
	dst.width  = width;
	dst.height = height;
	dst.bits   = src.bits;


	// 緩存
	vector<types> img_gauX(width*height*3);

	// 高斯模糊 X 軸
	const size_t r = gau_mat.size() / 2;

	int i, j, k;
#pragma omp parallel for private(i, j, k)
	for (j = 0; j < height; ++j) {
		for (i = 0; i < width; ++i) {
			double sumR = 0;
			double sumG = 0;
			double sumB = 0;
			for (k = 0; k < gau_mat.size(); ++k) {
				int idx = i-r + k;
				// idx超出邊緣處理
				if (idx < 0) {
					idx = 0;
				} else if (idx >(int)(width-1)) {
					idx = (width-1);
				}
				sumR += (double)src.raw_img[(j*width + idx)*3 + 0] * gau_mat[k];
				sumG += (double)src.raw_img[(j*width + idx)*3 + 1] * gau_mat[k];
				sumB += (double)src.raw_img[(j*width + idx)*3 + 2] * gau_mat[k];
			}
			img_gauX[(j*width + i)*3 + 0] = sumR;
			img_gauX[(j*width + i)*3 + 1] = sumG;
			img_gauX[(j*width + i)*3 + 2] = sumB;
		}
	}
	// 高斯模糊 Y 軸
#pragma omp parallel for private(i, j, k)
	for (j = 0; j < height; ++j) {
		for (i = 0; i < width; ++i) {
			double sumR = 0;
			double sumG = 0;
			double sumB = 0;
			for (k = 0; k < gau_mat.size(); ++k) {
				int idx = j-r + k;
				// idx超出邊緣處理
				if (idx < 0) {
					idx = 0;
				} else if (idx > (int)(height-1)) {
					idx = (height-1);
				}
				sumR += img_gauX[(idx*width + i)*3 + 0] * gau_mat[k];
				sumG += img_gauX[(idx*width + i)*3 + 1] * gau_mat[k];
				sumB += img_gauX[(idx*width + i)*3 + 2] * gau_mat[k];

			}
			dst.raw_img[(j*width + i)*3 + 0] = sumR;
			dst.raw_img[(j*width + i)*3 + 1] = sumG;
			dst.raw_img[(j*width + i)*3 + 2] = sumB;
		}
	}
}

void test_GauBlur() {
	Timer t1;
	// 讀取影像
	basic_ImgData img1, img2;
	Raw2Img::read_bmp(img1.raw_img, "sc02.bmp", &img1.width, &img1.height, &img1.bits);
	
	t1.start();
	GauBlur(img1, img2, 1.6, 0);
	t1.print(" WarpScale");

	Raw2Img::raw2bmp("test_GauBlur.bmp", img2.raw_img, img2.width, img2.height, 24);
}