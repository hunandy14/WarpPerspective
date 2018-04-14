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




class LaplacianBlending {  
private:  
	Mat_<Vec3f> left;  
	Mat_<Vec3f> right;  
	Mat_<float> blendMask;  

	vector<Mat_<Vec3f> > leftLapPyr,rightLapPyr,resultLapPyr;//Laplacian Pyramids  
	Mat leftHighestLevel, rightHighestLevel, resultHighestLevel;  
	vector<Mat_<Vec3f> > maskGaussianPyramid; //masks are 3-channels for easier multiplication with RGB  
	int levels;  
	void buildPyramids() {  
		buildLaplacianPyramid(left,leftLapPyr,leftHighestLevel);  
		buildLaplacianPyramid(right,rightLapPyr,rightHighestLevel);  
		buildGaussianPyramid();  
	}  
	void buildGaussianPyramid() {//金字塔內容爲每一層的掩模  
		assert(leftLapPyr.size()>0);  

		maskGaussianPyramid.clear();  
		Mat currentImg;  
		cvtColor(blendMask, currentImg, CV_GRAY2BGR);//store color img of blend mask into maskGaussianPyramid  
		maskGaussianPyramid.push_back(currentImg); //0-level  

		currentImg = blendMask;  
		for (int l=1; l<levels+1; l++) {  
			Mat _down;  
			if (leftLapPyr.size() > l)  
				pyrDown(currentImg, _down, leftLapPyr[l].size());  
			else  
				pyrDown(currentImg, _down, leftHighestLevel.size()); //lowest level  

			Mat down;  
			cvtColor(_down, down, CV_GRAY2BGR);  
			maskGaussianPyramid.push_back(down);//add color blend mask into mask Pyramid  
			currentImg = _down;  
		}  
	}  
	void buildLaplacianPyramid(const Mat& img, vector<Mat_<Vec3f> >& lapPyr, Mat& HighestLevel) {  
		lapPyr.clear();  
		Mat currentImg = img;  
		for (int l=0; l<levels; l++) {  
			Mat down,up;  
			pyrDown(currentImg, down);  
			pyrUp(down, up,currentImg.size());  
			Mat lap = currentImg - up;  
			lapPyr.push_back(lap);  
			currentImg = down;  
		}  
		currentImg.copyTo(HighestLevel);  
	}  
	Mat_<Vec3f> reconstructImgFromLapPyramid() {  
		//將左右laplacian圖像拼成的resultLapPyr金字塔中每一層  
		//從上到下插值放大並相加，即得blend圖像結果  
		Mat currentImg = resultHighestLevel;  
		for (int l=levels-1; l>=0; l--) {  
			Mat up;  

			pyrUp(currentImg, up, resultLapPyr[l].size());  
			currentImg = up + resultLapPyr[l];  
		}  
		return currentImg;  
	}  
	void blendLapPyrs() {  
		//獲得每層金字塔中直接用左右兩圖Laplacian變換拼成的圖像resultLapPyr  
		resultHighestLevel = 
			leftHighestLevel.mul(maskGaussianPyramid.back()) +  
			rightHighestLevel.mul(Scalar(1.0,1.0,1.0) - maskGaussianPyramid.back());  
		for (int l=0; l<levels; l++) {  
			Mat A = leftLapPyr[l].mul(maskGaussianPyramid[l]);  
			Mat antiMask = Scalar(1.0,1.0,1.0) - maskGaussianPyramid[l];  
			Mat B = rightLapPyr[l].mul(antiMask);  
			Mat_<Vec3f> blendedLevel = A + B;  

			resultLapPyr.push_back(blendedLevel);  
		}  
	}  

public:  
	LaplacianBlending(const Mat_<Vec3f>& _left, const Mat_<Vec3f>& _right, 
		const Mat_<float>& _blendMask, int _levels)://construct function, used in LaplacianBlending lb(l,r,m,4);  
		left(_left),right(_right),blendMask(_blendMask),levels(_levels)  
	{  
		assert(_left.size() == _right.size());  
		assert(_left.size() == _blendMask.size());  
		buildPyramids();  //construct Laplacian Pyramid and Gaussian Pyramid  
		blendLapPyrs();   //blend left & right Pyramids into one Pyramid  
	};  

	Mat_<Vec3f> blend() {  
		return reconstructImgFromLapPyramid();//reconstruct Image from Laplacian Pyramid  
	}  
};  
Mat_<Vec3f> LaplacianBlend(const Mat_<Vec3f>& l, const Mat_<Vec3f>& r, const Mat_<float>& m) {  
	LaplacianBlending lb(l,r,m,4);
	Mat_<Vec3f> blend=lb.blend();
	Mat re; 
	blend.convertTo(re,CV_8UC3,255);  
	return re;
}  
void mutBlender() {
	Mat l8u = imread("srcIMG\\LA.bmp");
	Mat r8u = imread("srcIMG\\LB.bmp");
	//imshow("left",l8u);   
	//imshow("right",r8u);

	Mat_<Vec3f> l;
	l8u.convertTo(l,CV_32F,1.0/255.0);
	Mat_<Vec3f> r;
	r8u.convertTo(r,CV_32F,1.0/255.0);

	Mat_<float> m(l.rows,l.cols,0.0);
	m(Range::all(),Range(0,m.cols/2)) = 1.0;

	Mat blend = LaplacianBlend(l, r, m); 
	imwrite("LS.bmp", blend);

	waitKey(0); 
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
	ImgData_resize(dst, newW, newH, 24);
	// 取出重疊區
	basic_ImgData cut1, cut2;
	getOverlap(src1, src2, cut1, cut2, corner);
	// 混合重疊區
	basic_ImgData blend;
	blendImg(blend, cut1, cut2);
	ImgData_write(blend, "___lapblend.bmp");
	// 合併三張圖片
	mergeOverlap(src1, src2, blend, dst, corner);
}




// 輸出圓柱投影AB
void cutWarpCyliImgAB(
	const basic_ImgData &src1, const basic_ImgData &src2, 
	basic_ImgData &dst, basic_ImgData &dst1, basic_ImgData &dst2, 
	const vector<int>& corner, int mode=0)
{
	// 偏移量
	int mx=corner[4];
	int my=corner[5];
	// 新圖大小
	int newH=corner[3]-corner[1]-my;
	int newW=corner[2]-corner[0]+mx;
	int img2mx = newW/2.0;
	// 重建大小
	ImgData_resize(dst, newW, newH, 24);
	ImgData_resize(dst1, newW-mx, newH, 24);
	ImgData_resize(dst2, newW-img2mx, newH, 24);
	// 兩張圖的高度偏差值
	int myA = my>0? 0:my;
	int myB = my<0? 0:my;


#pragma omp parallel for
	for (int j = 0; j < newH; j++) {
		for (int i = 0; i < newW; i++) {
			int src1idx;
			// 圖1
			if (i < corner[2]-corner[0]) {
				src1idx=(((j+myA)+corner[1])*src1.width +(i+corner[0])) *3;
				for (int  rgb = 0; rgb < 3; rgb++) {
					dst1.raw_img[(j*dst1.width +i)*3 +rgb] = src1.raw_img[src1idx+rgb];
				}
			}
			// 圖2
			if (i >= img2mx) {
				for (int  rgb = 0; rgb < 3; rgb++) {
					dst2.raw_img[(j*dst2.width +(i-img2mx)) *3+rgb] = 
						src2.raw_img[(((j+myB)+corner[1])*src1.width +((i-mx)+corner[0])) *3+rgb];
				}
			}
		}
	}
}

// 延伸圓柱投影A
void cutWarpCyliImgA(
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
	ImgData_resize(dst, newW, newH, 24);
	// 兩張圖的高度偏差值
	int myA = my>0? 0:my;
	int myB = my<0? 0:my;

	// 整張圖(test用)
	basic_ImgData& all=dst;

	// 像右延伸像素
	basic_ImgData right, right_gau;
	ImgData_resize(right, newH, 1, 24);
	//#pragma omp parallel for
	for (int j = 0; j < newH; j++) {
		int src1idx;
		for (int i = 0; i < newW; i++) {
			int idx=(j*dst.width +i) *3;
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

	GauBlur(right, right_gau, 1.6*3, 9);
	
	for (int j = 0; j < newH; j++) {
		for (int i = 0; i < newW; i++) {
			int idx=(j*dst.width +i) *3;
			int src1idx;
			// 向右拉平
			if (i >= /*mx*/ newW/2.0 ) {
				for (int  rgb = 0; rgb < 3; rgb++) {
					all.raw_img[idx+rgb] = right_gau.raw_img[j*3+rgb];
				}
			}
		}
	}
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
	ImgData_resize(dst, newW, newH, 24);

	//--------------------------------------------------
	// 整張圖(test用)
	// 擷取AB
	basic_ImgData dst1, dst2;
	cutWarpCyliImgAB(src1, src2, dst, dst1, dst2, corner);
	ImgData_write(dst1, "___dst1.bmp");
	ImgData_write(dst2, "___dst2.bmp");

	// 延長dstp
	basic_ImgData dst1p;
	cutWarpCyliImgA(src1, src2, dst1p, corner);
	ImgData_write(dst1p, "___dst1p.bmp");


	mutBlender();
	// 柏松混合
	Mat matsrc = imread("___dst2.bmp");
	Mat matsrcbg = imread("___dst1p.bmp");
	Mat src_mask = 255 * Mat::ones(matsrc.rows, matsrc.cols, matsrc.depth());
	Point mvPosi((matsrcbg.cols*2-matsrc.cols)/2, matsrcbg.rows / 2);
	// 一般柏松融合
	Timer t1;
	Mat normal_clone;
	seamlessClone(matsrc, matsrcbg, src_mask, mvPosi, normal_clone, NORMAL_CLONE);
	//imshow("normal_clone", normal_clone);
	imwrite("WarpCyliMuitBlend_pos.bmp", normal_clone);
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
	string name1 = "srcIMG\\ball_01.bmp", name2 = "srcIMG\\ball_02.bmp";
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
	//WarpCyliMuitBlend(matchImg, dst1, dst2, Ax, Ay);
	t1.print(" WarpCyliMuitBlend");
	//Raw2Img::raw2bmp("WarpCyliMuitBlend.bmp", matchImg.raw_img, matchImg.width, matchImg.height);
}