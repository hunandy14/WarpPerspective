/*****************************************************************
Name :
Date : 2018/03/15
By   : CharlotteHonG
Final: 2018/03/16
*****************************************************************/
#pragma once

using std::vector;
using std::string;

/*class Raw {
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
};*/

/*class _ImgRaw {
public:
	using uch = unsigned char;
public:
	_ImgRaw() = default;
	_ImgRaw(std::vector<uch> raw_img,	
		uint32_t width, uint32_t heigh,	uint16_t bits):
		raw_img(raw_img), width(width), height(heigh), bits() {}
public:
	void resize(uint32_t width, uint32_t heigh,	uint16_t bits) {
		raw_img.resize(width * heigh * bits/8.0);
		this->width=width;
		this->height=heigh;
		this->bits=bits;
	}
	size_t sizePix() {
		return width*height;
	}
	size_t size() {
		return sizePix()*bits;
	}
public:
	operator std::vector<uch>&() {
		return raw_img;
	}
	operator const std::vector<uch>& () const {
		return raw_img;
	}
	uch& operator[](size_t idx) {
		return raw_img[idx];
	}
	const uch& operator[](size_t idx) const {
		return raw_img[idx];
	}
public:
	std::vector<uch> raw_img;
	uint32_t width;
	uint32_t height;
	uint16_t bits;
};*/

void WarpPerspective(const basic_ImgData & src, basic_ImgData & dst, const vector<double>& H, bool clip);
void test1(string name, const vector<double>& HomogMat);

/*void WarpPerspective(const Raw & src, Raw & dst, const vector<double>& H, bool clip);
void test2(string name, const vector<double>& HomogMat);*/

void WarpPerspective(const vector<unsigned char>& src, const uint32_t srcW, const uint32_t srcH, vector<unsigned char>& dst, uint32_t & dstW, uint32_t & dstH, const vector<double>& H, bool clip);
void test3(string name, const vector<double>& HomogMat);

void WarpPers_Stitch(basic_ImgData & matchImg, const basic_ImgData & imgL, const basic_ImgData & imgR, const vector<double>& HomogMat);
void test_WarpPers_Stitch();
