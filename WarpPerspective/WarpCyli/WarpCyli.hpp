#pragma once

/*************************** Function Prototypes *****************************/
Raw DealWithImgData(Raw &srcdata, int width, int height, double R);

void warping(const std::vector<Raw> &inputArrays, float FL2, 
	std::vector<Raw> &Output, std::vector<fpoint> &upedge, std::vector<fpoint> &downedge);
