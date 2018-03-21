/*****************************************************************
Name :
Date : 2018/03/15
By   : CharlotteHonG
Final: 2018/03/16
*****************************************************************/
#pragma once

using std::vector;
using std::string;

void WarpPerspective(const basic_ImgData & src, basic_ImgData & dst, const vector<double>& H, bool clip);
void test1(string name, const vector<double>& HomogMat);
void test_WarpPers_Stitch();
