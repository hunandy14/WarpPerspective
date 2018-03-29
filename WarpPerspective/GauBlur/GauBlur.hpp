/*****************************************************************
Name :
Date : 2018/03/28
By   : CharlotteHonG
Final: 2018/03/28
*****************************************************************/

#pragma once


vector<float> gau_matrix(float p, size_t mat_len);

vector<float> gau_matrix2d(float p, size_t mat_len);

void GauBlur(const basic_ImgData & src, basic_ImgData & dst, float p, size_t mat_len);

void test_GauBlur();
