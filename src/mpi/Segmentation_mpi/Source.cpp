#include <iostream>
#include <math.h>
#include <stdlib.h>
#include<string.h>
#include<msclr\marshal_cppstd.h>
#include <ctime>// include this header 
#include<mpi.h>
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


void createImage(imageData imgData, int width, int height, int cluster_num)
{

	// Get the number of processors
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	// Get the rank of the processor
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//Declaring the variables

	int flag = 1;
	int channels = 3;
	int dims = width * height;
	int share = dims / size;
	int thresh = cluster_num * channels;
	int remainder = dims % size;
	int start_s, stop_s, TotalTime = 0;

	int* centers = new int[cluster_num * channels];
	int* cluster_size = new int[cluster_num];
	int* temp = new int[cluster_num * channels];
	int* clusters = new int[cluster_num * dims];
	int* sendcounts = new int[size];
	int* displs = new int[size];
	int* init_centers = new int[cluster_num * channels];

	//Calculate the sendcounts and displacements
	for (int i = 0; i < size - 1; i++) {
		sendcounts[i] = share;
		displs[i] = i * share;
	}

	//The last process takes care of the remaining pixels 
	//in case the image dimensions are not divisible by the number of processes
	sendcounts[size - 1] = share + remainder; 
	displs[size - 1] = (size - 1) * share;

	int* subRed = new int[sendcounts[rank]];
	int* subGreen = new int[sendcounts[rank]];
	int* subBlue = new int[sendcounts[rank]];

	int* local_centers = new int[cluster_num * channels];
	int* local_cluster_size = new int[cluster_num];
	int* subclusters = new int[sendcounts[rank]];

	//Initialize the centers
	if (rank == 0) {

		int point;
		srand(time(NULL));
		start_s = clock();
		for (int i = 0; i < cluster_num; i++)
		{
			point = rand() % (dims);
			init_centers[i * channels] = imgData.red[point];
			init_centers[i * channels + 1] = imgData.green[point];
			init_centers[i * channels + 2] = imgData.blue[point];
		}

		for (int i = 0; i < cluster_num * channels; i++)//temp is used later to check if the centers have significantly shifted
		{
			temp[i] = init_centers[i];
		}
	}

	//Scatter the image data to all processors
	MPI_Scatterv(imgData.red, sendcounts, displs, MPI_INT, subRed, sendcounts[rank], MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatterv(imgData.green, sendcounts, displs, MPI_INT, subGreen, sendcounts[rank], MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatterv(imgData.blue, sendcounts, displs, MPI_INT, subBlue, sendcounts[rank], MPI_INT, 0, MPI_COMM_WORLD);

	//Broadcast the initial centers to all processors
	MPI_Bcast(init_centers, cluster_num * channels, MPI_INT, 0, MPI_COMM_WORLD);


	for (int i = 0; i < cluster_num * channels; i++)
	{
		centers[i] = init_centers[i];
	}

	while (flag) //Iterate until the centers stop moving
	{
		//Reset the local cluster data
		for (int i = 0; i < cluster_num; i++) {
			local_centers[i * channels] = 0;
			local_centers[i * channels + 1] = 0;
			local_centers[i * channels + 2] = 0;
			local_cluster_size[i] = 0;
		}

		//Calculate each pixel's distance from each center and assign it to the closest center
		for (int i = 0; i < sendcounts[rank]; i++) {
			int min_dist = INT_MAX;
			int min_index = 0;

			for (int j = 0; j < cluster_num; j++) {
				int dist = sqrt(pow(subRed[i] - centers[j * channels], 2) + \
					pow(subGreen[i] - centers[j * channels + 1], 2) + \
					pow(subBlue[i] - centers[j * channels + 2], 2));

				if (dist < min_dist) {
					min_dist = dist;
					min_index = j;
				}
			}

			//Assign the pixel to the closest center and update the local cluster data
			subclusters[i] = min_index;
			local_cluster_size[min_index]++;
			local_centers[min_index * channels] += subRed[i];
			local_centers[min_index * channels + 1] += subGreen[i];
			local_centers[min_index * channels + 2] += subBlue[i];
		}

		//Reset the cluster data and calculate the new data by reducing the data of local clusters
		for (int i = 0; i < cluster_num; i++) {
			centers[i * channels] = 0;
			centers[i * channels + 1] = 0;
			centers[i * channels + 2] = 0;
			cluster_size[i] = 0;
		}

		MPI_Allreduce(local_cluster_size, cluster_size, cluster_num, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(local_centers, centers, cluster_num * channels, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		//Gather the pixel categorization from all processors
		MPI_Gatherv(subclusters, sendcounts[rank], MPI_INT, clusters, sendcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);

		//Calculate the new centers by averaging the pixels assigned to each center
		for (int i = 0; i < cluster_num; i++) {
			centers[i * channels] /= cluster_size[i];
			centers[i * channels + 1] /= cluster_size[i];
			centers[i * channels + 2] /= cluster_size[i];
		}

		if (rank == 0) {
			//Check if the centers have significantly shifted
			int diff = 0;
			for (int i = 0; i < cluster_num; i++) {
				diff += abs(centers[i * channels] - temp[i * channels]) + \
					abs(centers[i * channels + 1] - temp[i * channels + 1]) + \
					abs(centers[i * channels + 2] - temp[i * channels + 2]);
			}

			//If the centers have not shifted significantly, stop iterating; otherwise, update the temp variable
			if (diff < thresh) {
				flag = 0;
			}
			else {
				for (int i = 0; i < cluster_num; i++) {
					temp[i * channels] = centers[i * channels];
					temp[i * channels + 1] = centers[i * channels + 1];
					temp[i * channels + 2] = centers[i * channels + 2];
				}
			}
		}


		//Broadcast the flag to all processors
		MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0) {
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

		MyNewImage.Save("..//..//..//..//Data//Output//outputResrgb" + ".png");
		cout << "RGB Image Saved " << endl;

		stop_s = clock();
		TotalTime += (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000;

		cout << "RGB time: " << TotalTime << endl;
	}

	delete[] subRed;
	delete[] subGreen;
	delete[] subBlue;
	delete[] clusters;
	delete[] subclusters;
	delete[] temp;
	delete[] sendcounts;
	delete[] displs;
	delete[] local_centers;
	delete[] local_cluster_size;
	delete[] centers;
	delete[] cluster_size;
}

void createImageGrayscale(imageData imgData, int width, int height, int cluster_num)
{

	// Get the number of processors
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	// Get the rank of the processor
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//Declaring the variables

	int flag = 1;
	int dims = width * height;
	int share = dims / size;
	int thresh = cluster_num;
	int remainder = dims % size;
	int start_s, stop_s, TotalTime = 0;

	int* centers = new int[cluster_num];
	int* cluster_size = new int[cluster_num];
	int* temp = new int[cluster_num];
	int* clusters = new int[cluster_num * dims];
	int* sendcounts = new int[size];
	int* displs = new int[size];
	int* init_centers = new int[cluster_num];

	//Calculate the sendcounts and displacements
	for (int i = 0; i < size - 1; i++) {
		sendcounts[i] = share;
		displs[i] = i * share;
	}

	//The last process takes care of the remaining pixels 
	//in case the image dimensions are not divisible by the number of processes
	sendcounts[size - 1] = share + remainder;
	displs[size - 1] = (size - 1) * share;

	int* subGray = new int[sendcounts[rank]];

	int* local_centers = new int[cluster_num];
	int* local_cluster_size = new int[cluster_num];
	int* subclusters = new int[sendcounts[rank]];

	//Initialize the centers
	if (rank == 0) {

		int point;
		srand(time(NULL));
		start_s = clock();
		for (int i = 0; i < cluster_num; i++)
		{
			point = rand() % (dims);
			init_centers[i] = imgData.gray[point];
		}

		for (int i = 0; i < cluster_num; i++)//temp is used later to check if the centers have significantly shifted
		{
			temp[i] = init_centers[i];
		}
	}

	//Scatter the image data to all processors
	MPI_Scatterv(imgData.gray, sendcounts, displs, MPI_INT, subGray, sendcounts[rank], MPI_INT, 0, MPI_COMM_WORLD);

	//Broadcast the initial centers to all processors
	MPI_Bcast(init_centers, cluster_num, MPI_INT, 0, MPI_COMM_WORLD);


	for (int i = 0; i < cluster_num; i++)
	{
		centers[i] = init_centers[i];
	}

	while (flag) //Iterate until the centers stop moving
	{
		//Reset the local cluster data
		for (int i = 0; i < cluster_num; i++) {
			local_centers[i] = 0;
			local_cluster_size[i] = 0;
		}

		//Calculate each pixel's distance from each center and assign it to the closest center
		for (int i = 0; i < sendcounts[rank]; i++) {
			int min_dist = INT_MAX;
			int min_index = 0;

			for (int j = 0; j < cluster_num; j++) {
				int dist = abs(subGray[i] - centers[j]);

				if (dist < min_dist) {
					min_dist = dist;
					min_index = j;
				}
			}

			//Assign the pixel to the closest center and update the local cluster data
			subclusters[i] = min_index;
			local_cluster_size[min_index]++;
			local_centers[min_index] += subGray[i];
		}

		//Reset the cluster data and calculate the new data by reducing the data of local clusters
		for (int i = 0; i < cluster_num; i++) {
			centers[i] = 0;
			cluster_size[i] = 0;
		}

		MPI_Allreduce(local_cluster_size, cluster_size, cluster_num, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(local_centers, centers, cluster_num, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		//Gather the pixel categorization from all processors
		MPI_Gatherv(subclusters, sendcounts[rank], MPI_INT, clusters, sendcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);

		//Calculate the new centers by averaging the pixels assigned to each center
		for (int i = 0; i < cluster_num; i++) {

			if (cluster_size[i] == 0) { //In case a cluster has no points
				continue;
			}

			centers[i] /= cluster_size[i];
		}

		if (rank == 0) {
			//Check if the centers have significantly shifted
			int diff = 0;
			for (int i = 0; i < cluster_num; i++) {
				diff += centers[i] - temp[i];
			}

			//If the centers have not shifted significantly, stop iterating; otherwise, update the temp variable
			if (diff < thresh) {
				flag = 0;
			}
			else {
				for (int i = 0; i < cluster_num; i++) {
					temp[i] = centers[i];
				}
			}
		}


		//Broadcast the flag to all processors
		MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0) {
		System::Drawing::Bitmap MyNewImage(width, height);
		for (int i = 0; i < dims; i++) {
			System::Drawing::Color c = System::Drawing::Color::FromArgb(\
				centers[clusters[i]], \
				centers[clusters[i]], \
				centers[clusters[i]]);
			int row = i / width;
			int col = i % width;
			MyNewImage.SetPixel(col, row, c);
		}

		MyNewImage.Save("..//..//..//..//Data//Output//outputRes" + "Grayscale" + ".png");
		cout << "Grayscale Image Saved " << endl;

		stop_s = clock();
		TotalTime += (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000;

		cout << "Grayscale time: " << TotalTime << endl;
	}

	delete[] clusters;
	delete[] subclusters;
	delete[] temp;
	delete[] sendcounts;
	delete[] displs;
	delete[] local_centers;
	delete[] local_cluster_size;
	delete[] centers;
	delete[] cluster_size;
}

int main(int argc, char** argv)
{
	int ImageWidth = 4, ImageHeight = 4;
	int cluster_num = 3;
	System::String^ imagePath;
	std::string path;
	//Grayscale_Segmentation.jpg
	string img = "Grayscale_Segmentation.jpg";

	if (argc >= 3) {
		img = argv[1];
		cluster_num = atoi(argv[2]);
	}
	else if (argc == 2) {
		img = argv[1];
	}

	path = "..//..//..//..//Data//Input//" + img;
	imagePath = marshal_as<System::String^>(path);
	imageData imgData = inputImage(&ImageWidth, &ImageHeight, imagePath);

	MPI_Init(NULL, NULL);
	createImage(imgData, ImageWidth, ImageHeight, cluster_num);
	createImageGrayscale(imgData, ImageWidth, ImageHeight, cluster_num);
	MPI_Finalize();

	delete[] imgData.red;
	delete[] imgData.green;
	delete[] imgData.blue;
	delete[] imgData.gray;

	return 0;
}