#pragma once

/*************************** Function Prototypes *****************************/
void DealWithImgData(Raw &dst, const Raw &src, double R);

void warping(const std::vector<Raw> &inputArrays, float FL2, 
	std::vector<Raw> &Output, std::vector<fpoint> &upedge, std::vector<fpoint> &downedge);
