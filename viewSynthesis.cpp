#define _CRT_SECURE_NO_WARNINGS
#include "bmp.h"		//	Simple .bmp library
#include<iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

#define Baseline 30.0
#define Focal_Length 100
#define Image_Width 35.0
#define Image_Height 35.0
#define Resolution_Row 512
#define Resolution_Col 512
#define View_Grid_Row 9
#define View_Grid_Col 9

struct Point3d
{
	double x;
	double y;
	double z;
	Point3d(double x_, double y_, double z_) :x(x_), y(y_), z(z_) {}
};

struct Point2d
{
	double x;
	double y;
	Point2d(double x_, double y_) :x(x_), y(y_) {}
};

struct RGB
{
	unsigned char R;
	unsigned char G;
	unsigned char B;
	RGB(unsigned char _R, unsigned char _G, unsigned char _B) : R(_R), G(_G), B(_B) {}
};

void getUpperLeft(int* x, int* y);
void Basic(Point3d& rayRGB, double Vx, double Vy, int r, int c, vector<Bitmap>& viewImageList);
void Enhance(Point3d& rayRGB, const Point2d& intersectionP, int r, int c, vector<Bitmap>& viewImageList, int targetFocalLen);
void neighborInterpolation(Bitmap& image, RGB& tmp, int r, int c, int targetFocalLen);
void BilinearInterpolation(Point3d& rayRGB, vector<RGB>& rgb, double alpha, double beta);

int main(int argc, char** argv)
{
	/*
	if (argc < 5 || argc > 6)
	{
		cout << "Arguments prompt: viewSynthesis.exe <LF_dir> <X Y Z> OR: viewSynthesis.exe <LF_dir> <X Y Z> <focal_length>" << endl;
		return 0;
	}
	*/
	string LFDir = R"(C:\Users\danci\Downloads\ViewSynthesis-master\ViewSynthesis\LF_views)";
	double Vx = 0, Vy = 0, Vz = 100;
	double targetFocalLen = 100;
	if (argc == 6)
	{
		targetFocalLen = stod(argv[5]);
	}


	vector<Bitmap> viewImageList;
	//! loading light field views
	for (int i = 0; i < View_Grid_Col * View_Grid_Row; i++)
	{
		char name[128];
		sprintf(name, "/cam%03d.bmp", i + 1);
		string filePath = LFDir + name;
		Bitmap view_i(filePath.c_str());
		viewImageList.push_back(view_i);
	}

	Bitmap targetView(Resolution_Col, Resolution_Row);
	cout << "Synthesizing image from viewpoint: (" << Vx << "," << Vy << "," << Vz << ") with focal length: " << targetFocalLen << endl;
	cout << Resolution_Col << ' ' << Resolution_Row << endl;
	
	for (int r = 0; r < Resolution_Row; r++)
	{
		for (int c = 0; c < Resolution_Col; c++)
		{
			Point3d rayRGB(0, 0, 0);
			Point2d topLeft(Vx + Image_Width / 2, Vy - Image_Height / 2);
			Point2d line(Vx - topLeft.x, Vy - topLeft.y);
			Point2d intersection2();
			
			Point3d directionV(c - ((Resolution_Col - 1.) / 2.), r - ((Resolution_Row - 1.) / 2.), -targetFocalLen * (Resolution_Col / Image_Width));
			cout << directionV.x << ' ' << directionV.y << ' ' << directionV.z << endl;
			double k = (-Vz) / directionV.z;
			Point2d intersectionP(Vx + k * directionV.x, Vy + k * directionV.y);
			
			Enhance(rayRGB, intersectionP, r, c, viewImageList, targetFocalLen);

			//! record the resampled pixel value
			targetView.setColor(c, r, (unsigned char)rayRGB.x, (unsigned char)rayRGB.y, (unsigned char)rayRGB.z);
		}
	}
	string savePath = "newView.bmp";
	targetView.save(savePath.c_str());
	cout << "Result saved!" << endl;
	return 0;
}


void getUpperLeft(int* x, int* y)
{
	*x /= Baseline; *y /= Baseline;
	if (*x < 0) *x -= 1;
	if (*y > 0) *y += 1;
	*x += 4; *y += 4;
}


void Enhance(Point3d& rayRGB, const Point2d& intersectionP, int r, int c, vector<Bitmap>& viewImageList, int targetFocalLen)
{
	int x = (int)intersectionP.x; int y = (int)intersectionP.y;
	vector<pair<int, int>> neighbors;
	getUpperLeft(&x, &y);
	if ((int)intersectionP.x == 120)	x = 7;
	else if ((int)intersectionP.x == -120) x = 0;
	if ((int)intersectionP.y == 120) y = 7;
	else if ((int)intersectionP.y == -120) y = 0;
	if (x < 0 || 8 <= x || y <= 0 || 8 < y)
		return;
	neighbors.push_back(make_pair(x, y));
	neighbors.push_back(make_pair(x + 1, y));
	neighbors.push_back(make_pair(x, y - 1));
	neighbors.push_back(make_pair(x + 1, y - 1));
	vector<RGB> neighborPixel;
	for (int i = 0; i < 4; i++)
	{
		RGB t(0, 0, 0);
		neighborPixel.push_back(t);
	}
	double alpha = abs(intersectionP.x / 30 - (x - 4.)); double beta = abs((y - 4.) - intersectionP.y / 30);
	for (int i = 0; i < 4; i++)
		neighborInterpolation(viewImageList[(8 - neighbors[i].second) * View_Grid_Row + neighbors[i].first], neighborPixel[i], r, c, targetFocalLen);
	BilinearInterpolation(rayRGB, neighborPixel, alpha, beta);
}

void neighborInterpolation(Bitmap& image, RGB& res, int r, int c, int targetFocalLen)
{
	Point3d directionV(c - ((Resolution_Col - 1.) / 2.), r - ((Resolution_Row - 1.) / 2.), -targetFocalLen);
	double k = (-Focal_Length) / directionV.z;
	Point2d intersectionP(k * directionV.x + (Resolution_Col - 1.) / 2., k * directionV.y + (Resolution_Col - 1.) / 2.);
	int x = (int)intersectionP.x; int y = (int)intersectionP.y + 1;
	if (x < 0 || Resolution_Col - 1 <= x || y <= 0 || Resolution_Row - 1 < y)
	{
		res.R = 0; res.G = 0; res.B = 0;
		return;
	}
	double alpha = abs(intersectionP.x - x); double beta = abs(y - intersectionP.y);

	vector<RGB> neighborPoint;
	for (int i = 0; i < 4; i++)
	{
		RGB t(0, 0, 0);
		neighborPoint.push_back(t);
	}
	image.getColor(x, y, neighborPoint[0].R, neighborPoint[0].G, neighborPoint[0].B);
	image.getColor(x + 1, y, neighborPoint[1].R, neighborPoint[1].G, neighborPoint[1].B);
	image.getColor(x, y - 1, neighborPoint[2].R, neighborPoint[2].G, neighborPoint[2].B);
	image.getColor(x + 1, y - 1, neighborPoint[3].R, neighborPoint[3].G, neighborPoint[3].B);
	Point3d tmp(0, 0, 0);
	BilinearInterpolation(tmp, neighborPoint, alpha, beta);
	res.R = (unsigned char)tmp.x; res.G = (unsigned char)tmp.y; res.B = (unsigned char)tmp.z;
}

void BilinearInterpolation(Point3d& rayRGB, vector<RGB>& rgb, double alpha, double beta)
{
	Point3d p((1. - alpha) * rgb[0].R + alpha * rgb[1].R, (1. - alpha) * rgb[0].G + alpha * rgb[1].G, (1. - alpha) * rgb[0].B + alpha * rgb[1].B);
	Point3d q((1. - alpha) * rgb[2].R + alpha * rgb[3].R, (1. - alpha) * rgb[2].G + alpha * rgb[3].G, (1. - alpha) * rgb[2].B + alpha * rgb[3].B);
	Point3d t((1. - beta) * p.x + beta * q.x, (1. - beta) * p.y + beta * q.y, (1. - beta) * p.z + beta * q.z);
	rayRGB.x = t.x; rayRGB.y = t.y; rayRGB.z = t.z;
}