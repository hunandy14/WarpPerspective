/*****************************************************************
Name :
Date : 2017/06/11
By   : CharlotteHonG
Final: 2017/06/14
*****************************************************************/
#include <iostream>
#include <vector>
#include "Raw2Img.hpp"
#include <timer.hpp>
using namespace std;

class Raw {
public:
	Raw() = default;
	Raw(size_t w, size_t h) :
		RGB(w*h * 3), col(w), row(h) {}
	void resize(size_t w, size_t h){
		RGB.resize(w*h*3);
		col=w, row=h;
	}
	int getCol() const { return (int)col; }
	int getRow() const { return (int)row; }
public:
	std::vector<unsigned char> RGB;
protected:
	size_t col;
	size_t row;
};

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

// 輸入 dst 座標, 反轉 scr 輸出.
void WarpPerspective_CoorTranfer(const vector<double>& HomogMat, double& x, double& y) {
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
void _WarpPerspective(const Raw &src, Raw &dst, const vector<double> &H)
{
	dst.resize(src.getCol(), src.getRow());
	//-------------------------------------
	Color color;
	int d_Row = dst.getRow();
	int d_Col = dst.getCol();
	int s_Row = src.getRow();
	int s_Col = src.getCol();
	unsigned char *dst_RGB;
	//-------------------------------------
	int j, i;
	double x, y;

#pragma omp parallel for private(i, j, x, y)
	for (j = 0; j < d_Row; ++j)
	{
		dst_RGB = &dst.RGB[j * d_Col * 3];
		for (i = 0; i < d_Col; ++i)
		{
			x = i, y = j;
			WarpPerspective_CoorTranfer(H, x, y);

			if (x < (double)s_Col && x >= 0.0 && y < (double)s_Row && y >= 0.0)
			{
				color = bilinear(src, x, y);

				dst_RGB[i * 3 + 0] = color.R;
				dst_RGB[i * 3 + 1] = color.G;
				dst_RGB[i * 3 + 2] = color.B;
			}
		}
	}
}
//================================================================
int main(int argc, char const *argv[]) {
	vector<unsigned char> raw_img;
	uint32_t weidth, heigh;
	uint16_t bits;
	size_t imgSize, imgPixSize;
	/* 讀寫 Bmp */
	/*Raw2Img::read_bmp(raw_img, "ImgInput/Einstein.bmp", &weidth, &heigh, &bits);
	Raw2Img::write_raw("ImgOutput/out_Einstein.raw", raw_img);
	Raw2Img::raw2bmp("ImgOutput/out_Einstein.bmp", raw_img, weidth, heigh, bits);*/
	/* 轉灰階 */
	/*Raw2Img::raw2gray(raw_img);
	Raw2Img::raw2bmp("ImgOutput/out_Einstein_gray.bmp", raw_img, weidth, heigh, 8);*/

	Raw2Img::read_bmp(raw_img, "kanna.bmp", &weidth, &heigh, &bits);
	imgSize = weidth * heigh * bits/8.0;
	imgPixSize = weidth * heigh/8.0;

	vector<double> HomogMat{
		 0.708484   ,  0.00428145 , 245.901,
		-0.103356   ,  0.888676   , 31.6815,
		-0.000390072, -1.61619e-05, 1
	};

	Raw src(weidth, heigh), dst;
	src.RGB=raw_img;

	Timer t1;
	_WarpPerspective(src, dst, HomogMat);
	t1.print("_WarpPerspective");

	Raw2Img::raw2bmp("kanna2.bmp", dst.RGB, dst.getCol(), dst.getRow());

	return 0;
}
//================================================================
