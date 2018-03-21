#pragma once

/*************************** Function Prototypes *****************************/
void DealWithImgData(Raw& dst, const Raw& src, double R);
void WarpCylindrical(basic_ImgData &dst, const basic_ImgData &src, 
	double R ,int mx=0, int my=0, double edge=0.5);

void warping(const std::vector<Raw> &inputArrays, float FL2, 
	std::vector<Raw> &Output, std::vector<fpoint> &upedge, std::vector<fpoint> &downedge);
