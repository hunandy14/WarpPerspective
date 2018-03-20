/*****************************************************************
Name : 
Date : 2018/01/04
By   : CharlotteHonG
Final: 2018/01/04
*****************************************************************/
#pragma once

struct fpoint{
public:
	fpoint() :x(0.0), y(0.0) {}
	fpoint(float p_x, float p_y) : x(p_x), y(p_y) {}
	float x;
	float y;
};

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

class Image_Data {
public:
	Image_Data() {}
	Image_Data(size_t w, size_t h) : col(w), row(h), array(w * h, 0.0) {}
	~Image_Data(){
		array.clear();
	}
	float& operator[](size_t idx){
		return array[idx];
	}
	const float& operator[](size_t idx) const{
		return array[idx];
	}
	size_t getArraysize(){
		return array.size();
	}
	int getCol(){
		return (int)col;
	}
	int getRow(){
		return (int)row;
	}
private:
	std::vector<float> array;
	size_t col;
	size_t row;
};

struct Feature {
private:
	using Desc = std::vector<std::vector<std::vector<float>>>;
public:
	float size;                     // 階
	int kai;                        // 層
	float sigmaOCT;                 // 高斯模糊係數
	int x, y;                       // 各所在階層的座標
	float mm;                       // 強度
	int sita;                       // 包含主方向與負方向的角度
	Feature* nextptr = nullptr;     // 下一個鏈結點
									// 描述子相關
	float descr[128] = {};          // 統計完成後的描述子
	int d = 0;                      // 描述子長度                    
	Feature* fwd_match = nullptr;   // 匹配點
	fpoint img_pt;                  // 縮放回原圖的點
	fpoint mdl_pt;                  // 
	void* feature_data = nullptr;   // rob函式運算時的暫存
									// 顏穎
	float d_l;                      // 這個不知道幹嘛的，可能是測試用的
public:
	float rX() const {return x/size;}
	float rY() const {return y/size;}
};

// 保存檢測的相關數據
struct detection_data
{
	int r;
	int c;
	int octv;
	int intvl;
	float subintvl;
	float scl_octv;
};
