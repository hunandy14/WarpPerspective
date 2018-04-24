/*****************************************************************
Name :
Date : 2018/03/15
By   : CharlotteHonG
Final: 2018/03/16
*****************************************************************/
#pragma once

using std::vector;
using std::string;

// 其他細節
void WarpPerspective(const basic_ImgData & src, basic_ImgData & dst, const vector<double>& H, bool clip);
void WarpPerspective_cut(const basic_ImgData & src, basic_ImgData & dst, const vector<double>& H, bool clip);
void AlphaBlend(basic_ImgData & matchImg, const basic_ImgData & imgL, const basic_ImgData & imgR);
void PasteBlend(basic_ImgData & matchImg, const basic_ImgData & imgL, const basic_ImgData & imgR);

// 主要函式
void WarpPers_Stitch(basic_ImgData & matchImg, const basic_ImgData & imgL, const basic_ImgData & imgR, const vector<double>& HomogMat);

// 測試函式
void test1(string name, const vector<double>& HomogMat);
void test_WarpPers_Stitch();
void test_WarpPers_Stitch2(string name1, string name2);
