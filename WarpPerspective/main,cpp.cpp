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

// 線性取值

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

void WarpPerspective(const Raw &src, Raw &dst, const vector<double> &H)
{
	dst.resize(src.getCol(), src.getRow());
	//-------------------------------------
	int d_Row = dst.getRow();
	int d_Col = dst.getCol();
	int s_Row = src.getRow();
	int s_Col = src.getCol();
	//-------------------------------------
	int j, i;
	double x, y;
#pragma omp parallel for private(i, j, x, y)
	for (j = 0; j < d_Row; ++j) {
		for (i = 0; i < d_Col; ++i){
			x = i, y = j;
			WarpPerspective_CoorTranfer(H, x, y);
			if ((x <= (double)s_Col-1.0 and x >= 0.0) and
				(y <= (double)s_Row-1.0 and y >= 0.0))
			{
				dst.RGB[(j*d_Col + i)*3 + 0] = atBilinear_rgb(src.RGB, src.getCol(), y, x, 0);
				dst.RGB[(j*d_Col + i)*3 + 1] = atBilinear_rgb(src.RGB, src.getCol(), y, x, 1);
				dst.RGB[(j*d_Col + i)*3 + 2] = atBilinear_rgb(src.RGB, src.getCol(), y, x, 2);

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
	WarpPerspective(src, dst, HomogMat);
	t1.print("_WarpPerspective");

	Raw2Img::raw2bmp("kanna2.bmp", dst.RGB, dst.getCol(), dst.getRow());

	return 0;
}
//================================================================
