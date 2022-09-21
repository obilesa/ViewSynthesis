#define _CRT_SECURE_NO_DEPRECATE
#include "bmp.h"
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <algorithm>

#define SafeFree(p)  { if(p) { free(p); p = NULL; } }


Bitmap::Bitmap()
{
	memset(this, 0, sizeof(Bitmap));
}



Bitmap::Bitmap(const char *path)
{
	memset(this, 0, sizeof(Bitmap));
	create(path);
}



Bitmap::Bitmap(int width, int height)
{
	memset(this, 0, sizeof(Bitmap));
	create(width, height);
}



Bitmap::~Bitmap()
{
//	SafeFree(_data);
}


bool Bitmap::create(int width, int height)
{

	SafeFree(_data);

	_width = width;
	_height = height;

	_data = (Color *)malloc(_width * _height * sizeof(Color));
	if (_data == NULL)
		return false;
	memset(_data, 0, _width * _height * sizeof(Color)); 

	return true;

}


bool Bitmap::create(const char *path)
{
	BmpHeaderInfo bhi;

	{
		FILE *fp = fopen(path, "rb");
		if (fp == NULL)
			return false;

		fread(&bhi, sizeof(BmpHeaderInfo), 1, fp);
		fclose(fp);
	}

	if (bhi.bfType != 'MB')
		return false;

	if (bhi.biCompression != 0)
		return false;

	if (bhi.biBitCount != 24)
		return false;

	if (!create(bhi.biWidth, bhi.biHeight))
		return false;

	{
		FILE *fp = fopen(path, "rb");
		fseek(fp, bhi.bfOffBits, SEEK_SET);

		int index = 0;
		for (int j = 0; j < _height; j++)
		{
			fread(&_data[index], sizeof(Color), _width, fp);
			fseek(fp, (4 - _width * 3 % 4) % 4, SEEK_CUR);
			index += _width;
		}
		fclose(fp);
	}

	return true;
}


bool Bitmap::save(const char *path)
{
	BmpHeaderInfo bhi;

	bhi.bfType = 'MB';
	bhi.bfSize = _width * _height * 3 * sizeof(unsigned char) + sizeof(bhi);
	bhi.bfReserved1 = 0;
	bhi.bfReserved2 = 0;
	bhi.bfOffBits = sizeof(bhi);

	bhi.biSize = 40;
	bhi.biWidth = _width;
	bhi.biHeight = _height;
	bhi.biPlanes = 1;
	bhi.biBitCount = 24;
	bhi.biCompression = 0;
	bhi.biSizeImage = 0;
	bhi.biXpelsPerMeter = 0;
	bhi.biYpelsPerMeter = 0;
	bhi.biClrUsed = 0;
	bhi.biClrImportant = 0;

	unsigned char pad[3] = {0, 0, 0};

	FILE * fp = fopen(path, "wb");
	if (fp == NULL) return false;

	fwrite(&bhi, sizeof(BmpHeaderInfo), 1, fp);
	int index = 0;
	for (int j = 0; j < _height; j++) {
		
		fwrite(&_data[index], sizeof(Color), _width, fp);
		fwrite(pad, sizeof(unsigned char), (4 - _width * 3 % 4) % 4, fp);
		index += _width;
	}
	fclose(fp);

	return(true);

}

bool Bitmap::resize(int maxlen)
{
	if (_width < maxlen && _height < maxlen)
		return true;
	int scale = (std::max(_width, _height) + maxlen - 1) / maxlen;

	int width = _width / scale;
	int height = _height / scale;

	double * _newdata = (double *)malloc(width * height * 3 * sizeof(double));
	if (_newdata == NULL)
		return false;
	memset(_newdata, 0, (width * height * 3 * sizeof(double)));

	for (int x = 0; x < width*scale; x++)
	{
		for (int y = 0; y < height*scale; y++)
		{
			int pos = x + y * _width;
			int newpos = x / scale + y / scale * width;
			_newdata[newpos*3] += (double)_data[pos].R / (scale*scale);
			_newdata[newpos*3+1] += (double)_data[pos].G / (scale*scale);
			_newdata[newpos*3+2] += (double)_data[pos].B / (scale*scale);
		}
	}

	SafeFree(_data);
	_data = (Color *)malloc(width * height * sizeof(Color));
	if (_data == NULL)
		return false;
	memset(_data, 0, (width * height * sizeof(Color)));

	for (int x = 0; x < width; x++)
	{
		for (int y = 0; y < height; y++)
		{
			int pos = x + y * width;
			_data[pos].R = (unsigned char)_newdata[pos*3];
			_data[pos].G = (unsigned char)_newdata[pos*3+1];
			_data[pos].B = (unsigned char)_newdata[pos*3+2];
		}
	}

	_width = width;
	_height = height;
	SafeFree(_newdata);
	return true;

}