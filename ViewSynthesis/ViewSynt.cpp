#define _CRT_SECURE_NO_WARNINGS

#include "bmp.h"        //	Simple .bmp library
#include<iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

#define Baseline 30.0
#define Focal_Length 100
#define Image_Width 35.0
#define Image_Height 35.0
#define Resolution_Row 512
#define Resolution_Col 512
#define View_Grid_Row 9
#define View_Grid_Col 9

struct Point3d {
    double x;
    double y;
    double z;

    Point3d(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
};

struct Point2d {
    double x;
    double y;

    Point2d(double x_, double y_) : x(x_), y(y_) {}
};

struct RGB {
    unsigned char R;
    unsigned char G;
    unsigned char B;

    RGB(unsigned char _R, unsigned char _G, unsigned char _B) : R(_R), G(_G), B(_B) {}

    RGB() : R(0), G(0), B(0) {}
};

class ViewGridPoint {
public:
    Point2d coordinate;
    Bitmap image;
    double distance;

    ViewGridPoint(Point2d _coordinate, const Bitmap& _image, double _distance) : coordinate(_coordinate), image(_image),
        distance(_distance) {}
    bool operator<(const ViewGridPoint& other) const {
        return distance < other.distance;
    }
};


void getPixel(Bitmap& picture, RGB& rayRGB, int r, int c, int focal);

void bilinearInterpolation(RGB& rayRGB, RGB rgb[], double alpha, double beta);

void getNeigboringCameras(RGB& rayRGB, Point3d topLeftCorner, vector<Bitmap> viewImageList, Point2d intersection, int r,
    int c, int focal);

void getNeigboringCameras(RGB& rayRGB, Point3d topLeftCorner, vector<Bitmap> viewImageList, Point2d intersection, int r,
    int c, int focal) {
    vector<ViewGridPoint> gridPoints;
    for (int i = 0; i < View_Grid_Row; i++) {
        for (int j = 0; j < View_Grid_Col; j++) {
            Point2d coordinate(topLeftCorner.x + Baseline * i, topLeftCorner.y + Baseline * j);
            ViewGridPoint gridPoint(coordinate, viewImageList.at(i * j), sqrt(pow(coordinate.x - intersection.x, 2) +
                pow(coordinate.y - intersection.y, 2)));

            gridPoints.push_back(gridPoint);
        }
    }
    std::sort(gridPoints.begin(), gridPoints.end());
    RGB pixels[4];
    double alpha = (abs(abs(round(intersection.x)) - abs(gridPoints.at(0).coordinate.x))) / Image_Width;
    double beta = (abs(abs(gridPoints.at(0).coordinate.y) - abs(round(intersection.y)))) / Image_Height;
    for (int i = 0; i < 4; i++) {
        getPixel(gridPoints.at(i).image, pixels[i], r, c, focal);
    }
    bilinearInterpolation(rayRGB, pixels, alpha, beta);

}

int main(int argc, char** argv) {

    if (argc < 5 || argc > 6)
    {
        cout << "Arguments prompt: viewSynthesis.exe <LF_dir> <X Y Z> OR: viewSynthesis.exe <LF_dir> <X Y Z> <focal_length>" << endl;
        return 0;
    }
    string LFDir = argv[1];
    double Vx = stod(argv[2]), Vy = stod(argv[3]), Vz = stod(argv[4]);
    double targetFocalLen = 100; // default focal length for "basic requirement" part
    if (argc == 6)
    {
        targetFocalLen = stod(argv[5]);
    }

    vector<Bitmap> viewImageList;
    //! loading light field views
    for (int i = 0; i < View_Grid_Col * View_Grid_Row; i++) {
        char name[128];
        sprintf(name, "/cam%03d.bmp", i + 1);
        string filePath = LFDir + name;
        Bitmap view_i(filePath.c_str());
        viewImageList.push_back(view_i);
    }

    Point3d topLeft(Vx - (Baseline * (View_Grid_Row / 2.0)), Vy - (Baseline * (View_Grid_Col / 2.0)),
        Vz - Focal_Length);
    Bitmap targetView(Resolution_Col, Resolution_Row);
    cout << "Synthesizing image from viewpoint: (" << Vx << "," << Vy << "," << Vz << ") with focal length: "
        << targetFocalLen << endl;
    for (int r = 0; r < Resolution_Row; r++) {
        for (int c = 0; c < Resolution_Col; c++) {
            /*Calculating the intersection point where the light ray intersects the ViewPlane*/
            Point3d line(((r - 0.5) * (Baseline / Resolution_Row)) + topLeft.x,
            topLeft.y - ((c - 0.5) * (Baseline / Resolution_Col)), topLeft.z);
            double t = Vz / line.z * -1;
            Point2d intersection(Vx + line.x * t, Vy + line.y * t);
            RGB rayRGB(0, 0, 0);
            getNeigboringCameras(rayRGB, topLeft, viewImageList, intersection, r, c, targetFocalLen);
            //! record the resampled pixel value
            targetView.setColor(r, c, (unsigned char)rayRGB.R, (unsigned char)rayRGB.G, (unsigned char)rayRGB.B);
        }
    }

    string savePath = "newView.bmp";
    targetView.save(savePath.c_str());
    cout << "Result saved!" << endl;
    return 0;
}

void getPixel(Bitmap& picture, RGB& rayRGB, int row, int col, int focalLenght) {
    Point3d directionVector(row - ((Resolution_Row) / (double)2), col - ((Resolution_Col) / (double)2), -focalLenght);
    double t = Focal_Length / directionVector.z * -1;
    Point2d intersection(t * directionVector.x + Resolution_Col / (double)2, t * directionVector.y + Resolution_Col / (double)2);
    int x = round(intersection.x);
    int y = round(intersection.y);
    if (x < 0 || Resolution_Col - 1 <= x || y <= 0 || Resolution_Row - 1 < y) {
        rayRGB = RGB(0, 0, 0);
    }
    else {
        double alpha = abs(intersection.x - x);
        double beta = abs(y - intersection.y);

        RGB pixel[4];
       
        picture.getColor(x, y - 1, pixel[1].R, pixel[1].G, pixel[1].B);
        picture.getColor(x, y, pixel[1].R, pixel[1].G, pixel[1].B);
        picture.getColor(x + 1, y - 1, pixel[2].R, pixel[2].G, pixel[2].B);
        picture.getColor(x + 1, y, pixel[3].R, pixel[3].G, pixel[3].B);

        bilinearInterpolation(rayRGB, pixel, alpha, beta);
    }
}

void bilinearInterpolation(RGB& rayRGB, RGB rgb[], double alpha, double beta) {
    RGB pi((1 - alpha) * rgb[0].R + alpha * rgb[1].R, (1 - alpha) * rgb[0].G + alpha * rgb[1].G,
        (1 - alpha) * rgb[0].B + alpha * rgb[1].B);
    RGB pi1((1 - alpha) * rgb[2].R + alpha * rgb[3].R, (1 - alpha) * rgb[2].G + alpha * rgb[3].G,
        (1 - alpha) * rgb[2].B + alpha * rgb[3].B);
    rayRGB = RGB((1 - beta) * pi.R + beta * pi1.B, (1 - beta) * pi.B + beta * pi1.B, (1 - beta) * pi.G + beta * pi1.G);
}