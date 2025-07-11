// LBG_algo.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include<iostream>
#include<fstream>
#include <string>
#include<cmath>
#include<Windows.h>
#include "vector"
#include <string>
#include <iterator>
#include <iostream>
#include <cassert>
#include <sstream>
# define ROW           7125
#define E 0.03
#define THRESHOLD 0.0001
using namespace std;

float code_book[32][13];
float universe[7125][13];
float cluster[32][7125][13];
int codebook_size = 1;
float current_dist, D = 11111.11;
int weight[13] = { 0, 1, 3, 7, 13, 19, 22, 25, 33, 42, 50, 56, 61 };


void split_code_book(){
	float temp_centroid1[13];
	float temp_centroid2[13];
	float** temp = new float*[2 * codebook_size];
	for (int i = 0; i < 2 * codebook_size; ++i)
		temp[i] = new float[13];

	if (codebook_size <= 32){
		int i, j, k;
		float centroid[13];
		for (i = 0; i < codebook_size; i++){
			for (j = 1; j <= 12; j++){
				centroid[j] = code_book[i][j];
			}
			for (j = 1; j <= 12; j++){
				temp_centroid1[j] = centroid[j] * (1 - E);

				temp[2 * i][j] = temp_centroid1[j];
				//cout << temp[2 * i][j] << " ";
			}

			for (j = 1; j <= 12; j++){
				temp_centroid2[j] = centroid[j] * (1 + E);
				temp[(2 * i) + 1][j] = temp_centroid2[j];
				//cout << temp[2 * i][j] << " ";
			}


		}
		codebook_size = codebook_size * 2;
		for (i = 0; i < codebook_size; i++){
			for (j = 1; j < 13; j++){
				code_book[i][j] = temp[i][j];

			}

		}

	}

}

void classify(){
	float min_temp;
	int i, j, k;
	float sum = 0.0;
	int index = -1;
	for (i = 0; i < codebook_size; i++){
		code_book[i][0] = 0;
	}
	for (i = 0; i < ROW; i++)
	{
		min_temp = 111111111;
		for (j = 0; j < codebook_size; j++)
		{
			sum = 0.0;
			for (k = 1; k <= 12; k++)
			{

				sum = sum + weight[k] * (universe[i][k] - code_book[j][k]) * (universe[i][k] - code_book[j][k]);
			}
			if (sum < min_temp)
			{
				index = j;
				min_temp = sum;

			}
		}


		for (j = 1; j <= 12; j++)
		{
			cluster[index][(int)code_book[index][0]][j] = universe[i][j];
		}
		code_book[index][0]++;
	}


}

void find_centroid(){
	int i, j, k;
	float centroid[13];
	for (i = 0; i < codebook_size; i++)
	{
		for (j = 1; j <= 12; j++)
			centroid[j] = 0;

		for (j = 0; j < code_book[i][0]; j++)
		{
			for (k = 1; k <= 12; k++)
			{
				centroid[k] = centroid[k] + cluster[i][j][k];
			}
		}
		for (j = 1; j <= 12; j++)
		{
			code_book[i][j] = centroid[j] / code_book[i][0];
		}
	}
}

void distortion(){
	float cluster_dist = 0.0;
	int i, j, k;
	float sum;
	for (i = 0; i < codebook_size; i++){
		cluster_dist = 0.0;
		for (j = 0; j < code_book[i][0]; j++)
		{
			sum = 0.0;
			for (k = 1; k <= 12; k++)
			{

				sum = sum + weight[k] * (code_book[i][k] - cluster[i][j][k]) * (code_book[i][k] - cluster[i][j][k]);
			}
			cluster_dist = cluster_dist + sum;
		}
		cluster_dist = cluster_dist / code_book[i][0];
		current_dist = current_dist + cluster_dist;
	}
	current_dist = current_dist / codebook_size;
}

void lbg()
{
	int i, j, k;
	codebook_size = 1;
	float old_distortion = 11111111.111;
	string in_path = "./write/ciFeatures.txt";
	std::fstream universe_file(in_path, std::ios_base::in);
	//making an array of the universe file
	for (i = 0; i < ROW; i++){
		for (j = 1; j < 13; j++){
			universe_file >> universe[i][j];
		}
	}
	//calculate centroid for the whole universe
	float centroid[13];
	for (i = 0; i < 13; i++){
		centroid[i] = 0;
	}
	for (i = 1; i < 13; i++){
		for (j = 0; j < ROW; j++){
			centroid[i] = centroid[i] + universe[j][i];
		}
		centroid[i] = centroid[i] / ROW;
	}

	for (i = 0; i < codebook_size; i++){
		for (j = 1; j <= 12; j++){
			code_book[i][j] = centroid[j];
			//cout << centroid[j] << " ";
		}
		//cout << endl;
	}

	//start of LBG algo
	ofstream report("./write/LBG_log.txt");
/*	cout << "codebook size" << codebook_size << endl;
	cout << "*****************************************\n";
	report << "codebook size" << codebook_size << endl;
	report << "*****************************************\n"; 
	for (i = 0; i < codebook_size; i++){
		for (j = 1; j < 13; j++){
			cout << code_book[i][j] << " ";
			report << code_book[i][j] << " ";
		}
		cout << endl;
		report << endl;
	}
	cout << endl;
	*/
	while (codebook_size<32){

		//split the codebook	
		split_code_book();

	/*	cout << "codebook size" << codebook_size << endl;
		cout << "*****************************************\n";
		report << "codebook size" << codebook_size << endl;
		report << "*****************************************\n";
		for (i = 0; i < codebook_size; i++){
			for (j = 1; j < 13; j++){
				cout << code_book[i][j] << " ";
				report << code_book[i][j] << " ";
			}
			cout << endl;
			report << endl;
		}
		cout << endl;*/

		while (D >= THRESHOLD){
			//classify vector
			classify();
			//find the centroids
			find_centroid();
			//compute distortion 
			distortion();
			cout << current_dist << endl;
			D = abs(current_dist - old_distortion);
			old_distortion = current_dist;

		}
		D = current_dist;

	}
	//cout << D;
	for (i = 0; i < codebook_size; i++){
		for (j = 1; j < 13; j++){
			cout << code_book[i][j] << " ";
			report << code_book[i][j] << " ";
		}
		cout << endl;
		report << endl;
	}
	cout << "********************************************************\n\n";
	cout << "Check the file LBG_log.txt file to see the iterations \n";
	cout << "********************************************************\n";

	report.close();
	report.clear();
	
	return;
}

