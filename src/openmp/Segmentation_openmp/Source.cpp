#include <iostream>
#include <math.h>
#include <stdlib.h>
#include<string.h>
#include<msclr\marshal_cppstd.h>
#include <ctime>// include this header 
#include <omp.h>
#include<limits.h>
#include<stralign.h>
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


void createImage(imageData imgData, int width, int height, int thread_num, int cluster_num)
{

	//Declaring the variables
	int start_s, stop_s, TotalTime = 0;
	start_s = clock();

	int channels = 3;
	int dims = width * height;
	int thresh = cluster_num * channels;
	int remainder = dims % thread_num;
	int* centers = new int[cluster_num * channels];
	int* means = new int[cluster_num * channels];
	int* cluster_size = new int[cluster_num];
	int* clusters = new int[dims];
	int flag = 1;

	//Initialize the centers, the means and cluster sizes
	int point;
	srand(time(NULL));
	for (int i = 0; i < cluster_num; i++)
	{
		point = rand() % (dims);
		centers[i * channels] = imgData.red[point];
		centers[i * channels + 1] = imgData.green[point];
		centers[i * channels + 2] = imgData.blue[point];

		for (int c = 0; c < channels; c++) {
			means[i * channels + c] = 0;
		}

		cluster_size[i] = 0;
	}

	//Set the number of threads
	omp_set_num_threads(thread_num);

	//parallel region
	#pragma omp parallel
	{
		int nthreads = omp_get_num_threads();
		int id = omp_get_thread_num();

		if (id == 0)
			cout << "Number of threads: " << nthreads << endl;

		int start = id * (dims / nthreads);
		int end = start + (dims / nthreads);

		//the last thread will take care of the remainder 
		//in case the number of threads doesn't divide the number of pixels
		if (id == nthreads - 1)
			end += remainder;

		//the main loop iterating until the clusters converge
		while (flag) {

			//each thread will be assigned a portion of the image
			#pragma omp for schedule(static)
			for (int i = 0; i < dims; i++) {
				int min = INT_MAX;
				int min_index = 0;
				int dist = 0;
				for (int j = 0; j < cluster_num; j++) {
					//calculate the distance between the pixel and the cluster center
					dist += pow((imgData.red[i] - centers[j * channels]), 2);
					dist += pow((imgData.green[i] - centers[j * channels + 1]), 2);
					dist += pow((imgData.blue[i] - centers[j * channels + 2]), 2);

					if (dist < min) {
						min = dist;
						min_index = j;
					}
					dist = 0;
				}
				//assign the pixel to the cluster with the minimum distance
				clusters[i] = min_index;

				//update the means and the cluster sizes
				//this is a critical region since cluster_size and means are shared variables
				#pragma omp critical
				{
					cluster_size[min_index]++;
					means[min_index * channels] += imgData.red[i];
					means[min_index * channels + 1] += imgData.green[i];
					means[min_index * channels + 2] += imgData.blue[i];
				}

			}

			//Only one thread will enter this region to check if the clusters converged and update the centers
			#pragma omp single
			{
				int diff = 0;
				for (int i = 0; i < cluster_num; i++)
				{
					for (int c = 0; c < channels; c++) {

						if (cluster_size[i] == 0) { //In case a cluster has no elements
							continue;
						}

						//Calculate the difference between the new and the old cluster centers
						diff += abs(centers[i * channels + c] - means[i * channels + c] / cluster_size[i]);

						//Update the old center of each cluster by calculating the mean of its elements
						centers[i * channels + c] = means[i * channels + c] / cluster_size[i];
					}
				}

				//Reset the means and the cluster sizes
				for (int i = 0; i < cluster_num; i++)
				{
					for (int c = 0; c < channels; c++) {
						means[i * channels + c] = 0;
					}
					cluster_size[i] = 0;
				}

				//Check if the clusters converged
				if (diff < thresh)
					flag = 0;
				else
					flag = 1;

			}
			//Wait for all threads to finish to make sure that the centers and flag are updated
			#pragma omp barrier
		}
	}

	System::Drawing::Bitmap rgbImage(width, height);

	//Reconstruct the new image 
	//by replacing each pixel's RGP values with those of its cluster centroid
	for (int i = 0; i < dims; i++) {
		System::Drawing::Color c = System::Drawing::Color::FromArgb(\
			centers[clusters[i] * channels], \
			centers[clusters[i] * channels + 1], \
			centers[clusters[i] * channels + 2]);

		int row = i / width;
		int col = i % width;
		rgbImage.SetPixel(col, row, c);
	}

	//Save the reconstructed image
	rgbImage.Save("..//..//..//Data//Output//outputRes" + "rgb" + ".png");
	cout << "RGB Image Saved " << endl;

	stop_s = clock();
	TotalTime += (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000;

	cout << "RGB time: " << TotalTime << endl;

	delete[] means;
	delete[] clusters;
	delete[] centers;
	delete[] cluster_size;
}

void createImageGrayscale(imageData imgData, int width, int height, int thread_num, int cluster_num)
{

	//Declaring the variables
	int start_s, stop_s, TotalTime = 0;
	start_s = clock();

	int dims = width * height;
	int thresh = cluster_num;
	int remainder = dims % thread_num;
	int* centers = new int[cluster_num];
	int* means = new int[cluster_num];
	int* cluster_size = new int[cluster_num];
	int* clusters = new int[dims];
	int flag = 1;

	//Initialize the centers, the means and cluster sizes
	int point;
	srand(time(NULL));
	for (int i = 0; i < cluster_num; i++)
	{
		point = rand() % (dims);
		centers[i] = imgData.gray[point];

		means[i] = 0;

		cluster_size[i] = 0;
	}

	//Set the number of threads
	omp_set_num_threads(thread_num);

	//parallel region
	#pragma omp parallel
	{
		int nthreads = omp_get_num_threads();
		int id = omp_get_thread_num();

		if (id == 0)
			cout << "Number of threads: " << nthreads << endl;

		int start = id * (dims / nthreads);
		int end = start + (dims / nthreads);

		//the last thread will take care of the remainder 
		//in case the number of threads doesn't divide the number of pixels
		if (id == nthreads - 1)
			end += remainder;

		//the main loop iterating until the clusters converge
		while (flag) {

			//each thread will be assigned a portion of the image
			#pragma omp for schedule(static)
			for (int i = 0; i < dims; i++) {
				int min = INT_MAX;
				int min_index = 0;
				int dist = 0;
				for (int j = 0; j < cluster_num; j++) {
					//calculate the distance between the pixel and the cluster center
					dist += abs(imgData.gray[i] - centers[j]);

					if (dist < min) {
						min = dist;
						min_index = j;
					}
					dist = 0;
				}
				//assign the pixel to the cluster with the minimum distance
				clusters[i] = min_index;

				//update the means and the cluster sizes
				//this is a critical region since cluster_size and means are shared variables
				#pragma omp critical
				{
					cluster_size[min_index]++;
					means[min_index] += imgData.gray[i];
				}

			}

			//Only one thread will enter this region to check if the clusters converged and update the centers
			#pragma omp single
			{
				int diff = 0;
				for (int i = 0; i < cluster_num; i++)
				{
					if (cluster_size[i] == 0) { //In case a cluster has no elements
						continue;
					}

					//Calculate the difference between the new and the old cluster centers
					diff += abs(centers[i] - means[i] / cluster_size[i]);

					//Update the old center of each cluster by calculating the mean of its elements
					centers[i] = means[i] / cluster_size[i];
				}

				//Reset the means and the cluster sizes
				for (int i = 0; i < cluster_num; i++)
				{
					means[i] = 0;
					cluster_size[i] = 0;
				}

				//Check if the clusters converged
				if (diff < thresh)
					flag = 0;
				else
					flag = 1;

			}
			//Wait for all threads to finish to make sure that the centers and flag are updated
			#pragma omp barrier
		}
	}

	System::Drawing::Bitmap rgbImage(width, height);

	//Reconstruct the new image 
	//by replacing each pixel's RGP values with those of its cluster centroid
	for (int i = 0; i < dims; i++) {
		System::Drawing::Color c = System::Drawing::Color::FromArgb(\
			centers[clusters[i]], \
			centers[clusters[i]], \
			centers[clusters[i]]);

		int row = i / width;
		int col = i % width;
		rgbImage.SetPixel(col, row, c);
	}

	//Save the reconstructed image
	rgbImage.Save("..//..//..//Data//Output//outputRes" + "grayscale" + ".png");
	cout << "Grayscale Image Saved " << endl;

	stop_s = clock();
	TotalTime += (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000;

	cout << "Grayscale time: " << TotalTime << endl;

	delete[] means;
	delete[] clusters;
	delete[] centers;
	delete[] cluster_size;
}

int main(int argc, char** argv)
{

	int ImageWidth = 4, ImageHeight = 4;
	int cluster_num = 3;
	int thread_num = 8;
	System::String^ imagePath;
	std::string path;
	//Grayscale_Segmentation.jpg
	string img = "Grayscale_Segmentation.jpg";

	if (argc >= 4) {
		thread_num = atoi(argv[1]);
		img = argv[2];
		cluster_num = atoi(argv[3]);
	}
	else if (argc == 3) {
		thread_num = atoi(argv[1]);
		img = argv[2];
	}

	path = "..//..//..//Data//Input//" + img;
	imagePath = marshal_as<System::String^>(path);
	imageData imgData = inputImage(&ImageWidth, &ImageHeight, imagePath);

	createImage(imgData, ImageWidth, ImageHeight, thread_num, cluster_num);
	createImageGrayscale(imgData, ImageWidth, ImageHeight, thread_num, cluster_num);

	delete[] imgData.red;
	delete[] imgData.green;
	delete[] imgData.blue;
	delete[] imgData.gray;

	return 0;
}