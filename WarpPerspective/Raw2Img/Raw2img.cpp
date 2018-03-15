/*****************************************************************
Name :
Date : 2017/06/14
By   : CharlotteHonG
Final: 2017/06/14
*****************************************************************/
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "Raw2img.hpp"

using namespace std;
using uch = unsigned char;
//----------------------------------------------------------------
// 寫 Bmp 檔
void Raw2Img::raw2bmp(
    string name, vector<uch>& raw,
    uint32_t width, uint32_t height, uint16_t bits)
{
    // 檔案資訊
    BmpFileHeader file_h = makeFH(width, height, bits);
    // 圖片資訊
    BmpInfoHeader info_h = makeIH(width, height, bits);
    // 寫檔
    ofstream bmp(name, ios::binary);
    bmp.exceptions(ifstream::failbit|ifstream::badbit);
    bmp << file_h << info_h;
    // 寫調色盤
    if(bits == 8) {
        for(unsigned i = 0; i < 256; ++i) {
            bmp << uch(i) << uch(i) << uch(i) << uch(0);
        }
    }
    // 寫入圖片資訊
	size_t realW = info_h.biWidth * info_h.biBitCount/8.0;
	size_t alig = (realW*3) % 4;
    for(int j = height-1; j >= 0; --j) {
        for(unsigned i = 0; i < width; ++i) {
            if(bits==24) {
                bmp << raw[(j*width+i)*3 + B];
                bmp << raw[(j*width+i)*3 + G];
                bmp << raw[(j*width+i)*3 + R];
            } else if(bits==8) {
                bmp << raw[(j*width+i)];
            }
        }
        // 對齊4byte
        for(unsigned i = 0; i < alig; ++i) {
            bmp << uch(0);
        }
    }
}
// 寫 Raw 檔
void Raw2Img::write_raw(std::string name, std::vector<uch>& raw) {
    std::ofstream raw_file(name.c_str(), std::ios::binary);
    raw_file.exceptions(std::ifstream::failbit|std::ifstream::badbit);
    raw_file.write(reinterpret_cast<char*>(raw.data()), raw.size());
}
//----------------------------------------------------------------
// 讀 Bmp 檔案
void Raw2Img::read_bmp(vector<uch>& raw, string name,
    uint32_t* width, uint32_t* height, uint16_t* bits) {
    ifstream bmp(name.c_str(), ios::binary);
    bmp.exceptions(ifstream::failbit|ifstream::badbit);
    bmp.seekg(0, ios::beg);
    // 讀檔頭
    BmpFileHeader file_h;
    bmp >> file_h;
    BmpInfoHeader info_h;
    bmp >> info_h;
    // 回傳資訊
    if (width  != nullptr && 
        height != nullptr && 
        bits   != nullptr)
    {
        *width  = info_h.biWidth;
        *height = info_h.biHeight;
        *bits   = info_h.biBitCount;
    }
    // 讀 Raw
    bmp.seekg(file_h.bfOffBits, ios::beg);
    raw.resize(info_h.biWidth * info_h.biHeight * (info_h.biBitCount/8));
	size_t realW = info_h.biWidth * info_h.biBitCount/8.0;
    size_t alig = (realW*3) % 4;
    char* p = reinterpret_cast<char*>(raw.data());
    for(int j = info_h.biHeight-1; j >= 0; --j) {
        for(unsigned i = 0; i < info_h.biWidth; ++i) {
            // 來源是 rgb
            if(info_h.biBitCount == 24) {
                bmp.read(p + j*info_h.biWidth*3+i*3 + B, 1);
                bmp.read(p + j*info_h.biWidth*3+i*3 + G, 1);
                bmp.read(p + j*info_h.biWidth*3+i*3 + R, 1);
            }
            // 來源是 gray
            else if(info_h.biBitCount == 8) {
                bmp.read(p + j*info_h.biWidth+i, 1);
            }
        }
        bmp.seekg(alig, ios::cur); // 跳開 4bite 對齊的空格
    }
}
// 讀 Raw 檔
void Raw2Img::read_raw(std::vector<uch>& raw, std::string name) {
    std::ifstream raw_file(name.c_str(), 
        std::ios::binary | std::ios::ate);
    raw_file.exceptions(std::ifstream::failbit|std::ifstream::badbit);
    raw.resize(static_cast<size_t>(raw_file.tellg()));
    raw_file.seekg(0, std::ios::beg);
    raw_file.read(reinterpret_cast<char*>(raw.data()), raw.size());
    raw_file.close();
}
