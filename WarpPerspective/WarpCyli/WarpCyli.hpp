/*****************************************************************
Name :
Date : 2018/03/15
By   : CharlotteHonG
Final: 2018/03/21
*****************************************************************/
#pragma once

void WarpCylindrical(basic_ImgData &dst, const basic_ImgData &src, 
	double R ,int mx=0, int my=0, double edge=0);
void test_WarpCyli_AlphaBlend();
void test_WarpCyli_MuitBlend();
