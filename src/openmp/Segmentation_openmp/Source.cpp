#include <iostream>
#include <math.h>
#include <stdlib.h>
#include<string.h>
#include<msclr\marshal_cppstd.h>
#include <ctime>// include this header 
#include <omp.h>
#include<limits.h>
#include<stralign.h>
#define NUM_THREADS 4
#pragma once

#using <mscorlib.dll>
#using <System.dll>
#using <System.Drawing.dll>
#using <System.Windows.Forms.dll>
using namespace std;
using namespace msclr::interop;

typedef struct imageData {
	int* gray;
	int* red;
	int* blue;
	int* green;
}imageData;


imageData inputImage(int* w, int* h, System::String^ imagePath) //put the size of image in w & h
{
	int* gray;


	int OriginalImageWidth, OriginalImageHeight;

	//*********************************************************Read Image and save it to local arrays*************************	
	//Read Image and save it to local arrayss

	System::Drawing::Bitmap BM(imagePath);

	OriginalImageWidth = BM.Width;
	OriginalImageHeight = BM.Height;
	*w = BM.Width;
	*h = BM.Height;
	int* Red = new int[BM.Height * BM.Width];
	int* Green = new int[BM.Height * BM.Width];
	int* Blue = new int[BM.Height * BM.Width];
	gray = new int[BM.Height * BM.Width];

	for (int i = 0; i < BM.Height; i++)
	{
		for (int j = 0; j < BM.Width; j++)
		{
			System::Drawing::Color c = BM.GetPixel(j, i);

			Red[i * BM.Width + j] = c.R;
			Blue[i * BM.Width + j] = c.B;
			Green[i * BM.Width + j] = c.G;

			gray[i * BM.Width + j] = ((c.R + c.B + c.G) / 3); //gray scale value equals the average of RGB values

		}

	}

	imageData imgData;
	imgData.gray = gray;
	imgData.red = Red;
	imgData.blue = Blue;
	imgData.green = Green;

	return imgData;
}


void createImage(imageData imgData, int width, int height, int index)
{

	//Declaring the variables
	int start_s, stop_s, TotalTime = 0;
	start_s = clock();

	int cluster_num = 3;
	int channels = 3;
	int dims = width * height;
	int thresh = cluster_num * channels;
	int remainder = dims % NUM_THREADS;
	int* centers = new int[cluster_num * channels];
	int* means = new int[cluster_num * channels];
	int* cluster_size = new int[cluster_num];
	int* clusters = new int[dims];
	int flag = 1;

	//Initialize the centers
	int point;
	srand(time(NULL));
	for (int i = 0; i < cluster_num; i++)
	{
		point = rand() % (dims);
		centers[i * channels] = imgData.red[point];
		centers[i * channels + 1] = imgData.green[point];
		centers[i * channels + 2] = imgData.blue[point];
		means[i * channels] = 0;
		means[i * channels + 1] = 0;
		means[i * channels + 2] = 0;
		cluster_size[i] = 0;
	}

	omp_set_num_threads(NUM_THREADS);

#pragma omp parallel
	{
		int nthreads = omp_get_num_threads();
		int id = omp_get_thread_num();

		if (id == 0)
			cout << "Number of threads: " << nthreads << endl;
#pragma omp barrier
		int start = id * (dims / nthreads);
		int end = start + (dims / nthreads);

		if (id == nthreads - 1)
			end += remainder;

		while (flag) {

#pragma omp for schedule(static)
			for (int i = 0; i < dims; i++) {
				int min = INT_MAX;
				int min_index = 0;
				int dist = 0;
				for (int j = 0; j < cluster_num; j++) {

					dist += pow((imgData.red[i] - centers[j * channels]), 2);
					dist += pow((imgData.green[i] - centers[j * channels + 1]), 2);
					dist += pow((imgData.blue[i] - centers[j * channels + 2]), 2);
					if (dist < min) {
						min = dist;
						min_index = j;
					}
					dist = 0;
				}

				clusters[i] = min_index;
#pragma omp critical
				{
					cluster_size[min_index]++;
					means[min_index * channels] += imgData.red[i];
					means[min_index * channels + 1] += imgData.green[i];
					means[min_index * channels + 2] += imgData.blue[i];
				}

			}

#pragma omp barrier
#pragma omp single
			{
				int diff = 0;
				for (int i = 0; i < cluster_num; i++)
				{
					diff += abs(centers[i * channels] - means[i * channels] / cluster_size[i]);
					diff += abs(centers[i * channels + 1] - means[i * channels + 1] / cluster_size[i]);
					diff += abs(centers[i * channels + 2] - means[i * channels + 2] / cluster_size[i]);
					centers[i * channels] = means[i * channels] / cluster_size[i];
					centers[i * channels + 1] = means[i * channels + 1] / cluster_size[i];
					centers[i * channels + 2] = means[i * channels + 2] / cluster_size[i];
				}

				for (int i = 0; i < cluster_num; i++)
				{
					means[i * channels] = 0;
					means[i * channels + 1] = 0;
					means[i * channels + 2] = 0;
					cluster_size[i] = 0;
				}


				if (diff < thresh)
					flag = 0;
				else
					flag = 1;

			}

#pragma omp barrier
		}
	}

	System::Drawing::Bitmap MyNewImage(width, height);
	for (int i = 0; i < dims; i++) {
		System::Drawing::Color c = System::Drawing::Color::FromArgb(\
			centers[clusters[i] * channels], \
			centers[clusters[i] * channels + 1], \
			centers[clusters[i] * channels + 2]);
		int row = i / width;
		int col = i % width;
		MyNewImage.SetPixel(col, row, c);
	}

	MyNewImage.Save("E://GitHub//Parallel_Image_Segmentation//Data//Output//outputRes" + index + ".png");
	cout << "result Image Saved " << index << endl;

	stop_s = clock();
	TotalTime += (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000;

	cout << "time: " << TotalTime << endl;

	delete[] means;
	delete[] clusters;
	delete[] centers;
	delete[] cluster_size;
}


int main()
{

	int ImageWidth = 4, ImageHeight = 4;

	System::String^ imagePath;
	std::string img;
	//Grayscale_Segmentation.jpg
	img = "E://GitHub//Parallel_Image_Segmentation//Data//Input//Grayscale_Segmentation.jpg";

	imagePath = marshal_as<System::String^>(img);
	imageData imgData = inputImage(&ImageWidth, &ImageHeight, imagePath);

	createImage(imgData, ImageWidth, ImageHeight, 1);


	delete[] imgData.red;
	delete[] imgData.green;
	delete[] imgData.blue;
	delete[] imgData.gray;

	return 0;
}