/**
 * @file opencv_canny.cpp
 */

#include <stdio.h>

#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

#include "opencv_canny.h"

using namespace cv;
using namespace std;

void cv_canny(uint8_t *dest, uint8_t *src, int width, int height)
{
	Mat src1(height, width, CV_8UC1, src);

	//src1 = imread("/tmp/image.jpeg", CV_LOAD_IMAGE_UNCHANGED);
	//
	printf(" new image %d %d\n", width, height);

	namedWindow("Original image", CV_WINDOW_AUTOSIZE);
	imshow("Original image", src1);

//	Mat gray, edge, draw;
//	cvtColor(src1, gray, CV_BGR2GRAY);

//	Canny(gray, edge, 32, 144, 3);

//	edge.convertTo(draw, CV_8U);
//	namedWindow("image", CV_WINDOW_AUTOSIZE);
//	imshow("image", draw);

	waitKey(0);
	//IplImage* image = cvCreateImage(cvSize(width, height), 8, 1);

	/*
	IplImage* image = cvCreateImage(cvSize(width, height), 8, 1);
	memcpy(image->imageData, src, width * height);

	namedWindow("Gray image", CV_WINDOW_AUTOSIZE);
	imshow("example", image);
	waitKey(0);
	*/

	//imwrite( "/tmp/Gray_Image.jpg", image);
}

