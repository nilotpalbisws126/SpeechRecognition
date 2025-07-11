#include "stdafx.h"
#include "stdafx.h"
#include<iostream>
#include<fstream>
#include <string>
#include<cmath>
#include<Windows.h>
#include "vector"
#include <sstream>
#include <algorithm>
#include<iomanip>
#define M 32
#define P 12
#define T 180 //obervation sequence
#define N 5
#define no_of_training_data 20
#define no_of_train_sample 20


using namespace std;

class HMMmodel{
public:
	long double pi[N];
	long double A[N][N];
	long double B[N][M];
};

class Codebook
{
public:
	int index;
	double ci[P];
};

class Zeta{
public:
	long double val[N][N];
};


class StateProbability{
public:
	long double alpha[N];
	long double beta[N];
	long double gamma[N];
};

vector <Codebook> codebook;
vector <double> temp_input;
vector <double> autoCorr;
vector <int> ObservationSeq;
vector <StateProbability> alpha_beta;
vector <long double>test_alpha;
vector <long double> test_probabilities;
vector <double> R;
vector <Zeta> zeta;
int window_size = 320;
int corr_size;
long int size;
double alpha[P];
long double pi[N];
long double pi_0[N] = { 1, 0, 0, 0, 0 };
long double A[N][N];
long double A_0[N][N] = { { 0.8, 0.2, 0, 0, 0 }, { 0, 0.8, 0.2, 0, 0 }, { 0, 0, 0.8, 0.2, 0 }, { 0, 0, 0, 0.8, 0.2 }, { 0, 0, 0, 0, 1 } };
long double B[N][M];
long double B_0[N][M] = { { 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000 },
{ 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000 },
{ 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000 },
{ 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000 },
{ 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000, 0.03125000 } };


void dc_norm()
{
	int i;
	size = temp_input.size();				// Length of input speech signal
	ofstream normalize;
	normalize.open("normalized_file.txt");

	std::fstream silent_file("silence.txt", std::ios_base::in);

	int a;
	double silence_sum = 0, silence_count = 0, dc_shift;
	while (silent_file >> a)
	{
		silence_sum = silence_sum + a;
		silence_count = silence_count + 1;

	}
	dc_shift = silence_sum / silence_count; //dc shift calculation
	//cout << "DC shift is :" << dc_shift << endl;
	silent_file.close();
	silent_file.clear();

	for (i = 0; i<size; i++)
	{
		temp_input[i] = temp_input[i] - dc_shift;
		// Mean adjustment to remove DC offset from the given speech signal 		
	}

	// Mean adjustment and normalization of "infile" vector
	double max_val = *max_element(temp_input.begin(), temp_input.end());
	double min_val = *min_element(temp_input.begin(), temp_input.end());
	max_val = abs(max_val); min_val = abs(min_val);
	if (min_val > max_val)
		max_val = min_val;
	for (i = 0; i < size; i++)
	{
		temp_input[i] = temp_input[i] * 5000 / max_val;
		normalize << temp_input[i] << endl;

	}

	normalize.close();

}

void crop(){
	std::fstream normalized_file("normalized_file.txt", std::ios_base::in);
	float c, sum;
	int count1 = 0, count2 = 0, i;
	int flag = 0;
	int end_ptr_close = 0;
	int zcr;
	bool initial_postion = false;
	double info_start = 0, info_content = 0, start_ptr = 0, end_ptr = 0, temp_end_ptr = 0;
	int ct = 0;
	while (!normalized_file.eof()){

		count1 = 0;
		sum = 0;
		zcr = 0;
		while (count1<100 && !normalized_file.eof()){
			normalized_file >> c;
			//ZCR calculation
			if (initial_postion && c < 0)
			{
				zcr++;
				initial_postion = false;
			}
			else if (!initial_postion && c>0){
				zcr++;
				initial_postion = true;
			}

			sum = sum + abs(c);
			count1++;
			count2++;

		}

		sum = sum / count1;
		//cout << sum << endl;
		//cout << zcr << endl << endl;


		if (flag == 1)
		{
			if ((sum > 380 || (sum<380 && zcr>25)) && (end_ptr_close == 0)){
				if (ct<1){
					end_ptr = end_ptr + 100;
					temp_end_ptr = end_ptr;
				}
				else{
					end_ptr = end_ptr + ct * 100;
					ct = 0;
				}

			}
			else{
				ct += 1;
				if (ct > 40){
					end_ptr = temp_end_ptr;
					end_ptr_close = 1;

				}
			}
		}
		else{
			if (sum>380) //if sum is greater for first frame (100 signal) shift the start marker to the end
			{
				flag = 1;
				start_ptr = count2 - 1;
			}
			else //else set both start and end marker to en of the frame
			{
				start_ptr = start_ptr + 100;
				end_ptr = start_ptr;
			}
		}

	}
	normalized_file.close();
	normalized_file.clear();
	ofstream cropped;
	cropped.open("cropped_file.txt");
	std::fstream input("normalized_file.txt", std::ios_base::in);

	string line;
	i = 0;
	//cout << info_start << endl;
	//cout << info_content << endl;
	while (i<start_ptr)
	{
		getline(input, line);
		i++;
	}
	while (i<(end_ptr))
	{

		input >> c;
		cropped << c << endl;
		i++;
	}

	cropped.close();
	input.close();
}

void crop1(int &folder_name, int &file_name){
	std::fstream normalized_file("normalized_file.txt", std::ios_base::in);
	float c, sum;
	int count1 = 0, count2 = 0, i;
	int flag = 0;
	int end_ptr_close = 0;
	int zcr;
	bool initial_postion = false;
	double info_start = 0, info_content = 0, start_ptr = 0, end_ptr = 0, temp_end_ptr = 0;
	int ct = 0;
	while (!normalized_file.eof()){

		count1 = 0;
		sum = 0;
		zcr = 0;
		while (count1<100 && !normalized_file.eof()){
			normalized_file >> c;
			//ZCR calculation
			if (initial_postion && c < 0)
			{
				zcr++;
				initial_postion = false;
			}
			else if (!initial_postion && c>0){
				zcr++;
				initial_postion = true;
			}

			sum = sum + abs(c);
			count1++;
			count2++;

		}

		sum = sum / count1;
		//cout << sum << endl;
		//cout << zcr << endl << endl;


		if (flag == 1)
		{
			if ((sum > 380 || (sum<380 && zcr>25)) && (end_ptr_close == 0)){
				if (ct<1){
					end_ptr = end_ptr + 100;
					temp_end_ptr = end_ptr;
				}
				else{
					end_ptr = end_ptr + ct * 100;
					ct = 0;
				}

			}
			else{
				ct += 1;
				if (ct > 40){
					end_ptr = temp_end_ptr;
					end_ptr_close = 1;

				}
			}
		}
		else{
			if (sum>380) //if sum is greater for first frame (100 signal) shift the start marker to the end
			{
				flag = 1;
				start_ptr = count2 - 1;
			}
			else //else set both start and end marker to en of the frame
			{
				start_ptr = start_ptr + 100;
				end_ptr = start_ptr;
			}
		}

	}
	normalized_file.close();
	normalized_file.clear();
	ostringstream in_name;
	in_name << "./cropped_train/" << folder_name << "/" << file_name << ".txt";
	ofstream cropped;
	string in_filename = in_name.str();
	cropped.open(in_filename);
	std::fstream input("normalized_file.txt", std::ios_base::in);

	string line;
	i = 0;
	//cout << info_start << endl;
	//cout << info_content << endl;
	while (i<start_ptr)
	{
		getline(input, line);
		i++;
	}
	while (i<(end_ptr))
	{

		input >> c;
		cropped << c << endl;
		i++;
	}

	cropped.close();
	input.close();
}


void autocorelation(vector<double> &frame)
{
	 long int segmentSize = window_size;
	vector <double> tempSegment;
	long int m, n;
	// voice speech in reverse order
	for (n = 0; n<segmentSize; n++)
		tempSegment.push_back(frame[segmentSize - 1 - n]);
	
	corr_size = segmentSize + segmentSize - 1;
	for (n = 0; n<corr_size; n++)
	{
		double temp = 0;
		for (m = 0; m<segmentSize; m++)
		{
			if (n - m >= 0 && n - m < segmentSize)
			{
				temp = temp + frame[n - m] * tempSegment[m];
			}
		}
		autoCorr.push_back(temp);
	}

	for (n = 0; n <= P; n++)
	{
		R.push_back(autoCorr[segmentSize - 1 + n]); //Saving the Ri s
		

	}
	tempSegment.clear();
	autoCorr.clear();
}

void calculate_alpha_i()
{
	//durbin's algo for calculating alpha
	int i, j;
	double K, E, temp;
	double k[P];
	for (i = 0; i<P; i++)
	{
		alpha[i] = 0;
		k[i] = 0;
	}

	E = R[0];
	for (i = 1; i <= P; i++)
	{
		temp = 0;
		for (j = 1; j <= i - 1; j++)
		{
			temp = temp + alpha[j - 1] * R[i - j];
		}
		K = (R[i] - temp) / E;
		k[i - 1] = K;
		for (j = 1; j <= i - 1; j++)
		{
			k[j - 1] = alpha[j - 1] - K*alpha[i - j - 1];
		}
		for (j = 0; j<P; j++)
		{
			alpha[j] = k[j];
		}
		E = (1 - K*K) * E;
	}
}


void calculate_c_i(ofstream &to_write)
{
	// Computing Ri and Ci
	int k = 0;
	double temp, c_i[P];
	std::fstream cropped_file("cropped_file.txt", std::ios_base::in);
	while (!cropped_file.eof()){
		cropped_file >> temp_input[k];
		k++;

	}
	int frame_size = window_size * 1 / 4;
	vector <double> frame;
	long int fileLen = temp_input.size();
	int count = 0;
	long int i, m, n;
	for (int i = 0; i < fileLen - window_size; i = i + frame_size)
	{
		for (n = 0; n<window_size; n++)
		{
			double win_val = 0.54 - 0.46 * cos(2 * 3.141*n / (window_size - 1));  //hamming window
			frame.push_back(win_val * temp_input[i + n]);
		}
		autocorelation(frame); //calling autocorelation for getting Ri values
		if (R[0] < 100)
		{
			R.clear();
			frame.clear();
			continue;
		}
		count = count + 1;
		double max_R = R[0];
		for (n = 0; n < R.size(); n++)
		{
			R[n] = R[n] / max_R;			// Normalizing the Ri values

		}
		calculate_alpha_i();	//calculating alpha using durbins method
		frame.clear();
		R.clear();
		// Calculating the Ci
		for (n = 1; n <= P; n++)
		{
			temp = 0;
			for (m = 1; m <= n - 1; m++)
			{
				temp = temp + (m / n)*c_i[m - 1] * alpha[n - m - 1];
			}
			c_i[n - 1] = alpha[n - 1] + temp;
			c_i[n - 1] = c_i[n - 1] * (1 + sin(3.141*(n - 1) / (P - 1))*P / 2);
		}
		for (n = 0; n<P; n++)
		{
			to_write << c_i[n] << " ";
			//cout << c_i[n] << endl;
		}
		to_write << endl;
	}

}
void loadCodebook()
{
	string path = "./write/LBG_log.txt";
	string line;
	double temp;
	ifstream codebookfile;
	codebookfile.open(path);
	int indx = 0;
	while (!codebookfile.eof())
	{
		indx = indx + 1;
		getline(codebookfile, line);
		stringstream string_s(line);
		Codebook var_temp;
		int i = 0;
		while (string_s >> temp)
		{
			var_temp.ci[i] = temp;
			i = i + 1;
		}
		var_temp.index = indx;									 // storing all ci
		codebook.push_back(var_temp);						 
	}
	codebookfile.close();
}




void train()
{
	double temp;
	cout << " ******************** Feature extraction ******************" << endl << endl;
	ostringstream out_file;
	out_file << "./write/ciFeatures.txt";
	// write the cepstral feature values in the file, name stored in out_name
	string out_filename = out_file.str();
	ofstream to_write;
	to_write.open(out_filename);

	for (int i = 0; i<4; i++)
	{
		for (int j = 1; j <= no_of_train_sample; j++)
		{
			cout << "Digit " << i << " from file " << j << endl;
			ostringstream in_name;
			in_name << "./digit_train/" << i << "/" << j << ".txt";
			string in_filename = in_name.str();
			ifstream infile_obj;
			infile_obj.open(in_filename);
			while (!infile_obj.eof())
			{
				infile_obj >> temp;
				temp_input.push_back(temp);
			}
			infile_obj.close();			
			dc_norm();				// dc shift and normalization
			crop1(i,j);				//remove unvoiced part
			calculate_c_i(to_write);	// calculate LPC and cepstral coefficients					
			temp_input.clear();
		}
	}
	to_write.close();
	cout << " Yepiiiiiiii you are done !!!!!!!!!!! " << endl;
	return;
}

long double alphabetasolution()
{
	// Initialization
	StateProbability P_O;
	for (int i = 0; i < N; i++)
	{
		
		P_O.alpha[i] = pi[i] * B[i][ObservationSeq[0] - 1];
		P_O.beta[i] = 1;
		P_O.gamma[i] = 0;
	}
	alpha_beta.push_back(P_O);

	// Induction	
	int length, k;
	if (ObservationSeq.size() <= T)
		length = ObservationSeq.size();
	else
		length = T + 1;


	for (k = 1; k < length; k++)
	{
		// Forward probability	
		for (int j = 0; j < N; j++)
		{
			long double temp_a = 0;
			for (int i = 0; i<N; i++)
			{
				temp_a = temp_a + alpha_beta[k - 1].alpha[i] * A[i][j];
			}
			P_O.alpha[j] = temp_a * B[j][ObservationSeq[k] - 1];
		}

		// Backward probability
		int indx = ObservationSeq[length - k] - 1;
		for (int i = 0; i < N; i++)
		{
			long double temp_b = 0;
			for (int j = 0; j<N; j++)
			{
				temp_b = temp_b + alpha_beta[k - 1].beta[j] * A[i][j] * B[j][indx];
			}
			P_O.beta[i] = temp_b;
			P_O.gamma[i] = 0;
		}
		alpha_beta.push_back(P_O);
	}

	// //The termintion step
	long double prob1 = 0, prob2 = 0;
	int temp = ObservationSeq[0] - 1;
	for (int i = 0; i < N; i++)
	{
		prob1 = prob1 + alpha_beta[k - 1].alpha[i];
		prob2 = prob2 + alpha_beta[k - 1].beta[i] * B[i][temp] * pi[i];
	}
//	cout << "\n\nP(O|lambda) is :" << prob1 << endl;

	return prob1;
}

void zetasolution3()
{
	Zeta zeta_temp;
	long double temp_sum = 0;
	int length;
	if (ObservationSeq.size() <= T) //taking maximum T observation sequence 
		length = ObservationSeq.size();
	else
		length = T ;
	for (int t = 1; t < length; t++)
	{
		temp_sum = 0;
		for (int i = 0; i < N; i++)
		{
			int indx = ObservationSeq[t] - 1;
			for (int j = 0; j < N; j++)
			{								 															
				zeta_temp.val[i][j] = alpha_beta[t - 1].alpha[i] * A[i][j] * B[j][indx] * alpha_beta[length - t].beta[j]; // Beta  are stored in in reverse order
				temp_sum = temp_sum + zeta_temp.val[i][j];
			}
		}

		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				zeta_temp.val[i][j] = zeta_temp.val[i][j] / temp_sum;
			}
		}
		zeta.push_back(zeta_temp);
	}
}


void ReEstimation() //Solution 3
{
	

	int length;
	if (ObservationSeq.size() <= T)
		length = ObservationSeq.size();
	else
		length = T + 1;
	vector <StateProbability> gama;
	for (int t = 0; t < length - 1; t++)
	{
		StateProbability tempGama;
		for (int i = 0; i < N; i++)
		{
			long double temp = 0;
			for (int j = 0; j < N; j++)
			{
				temp = temp + zeta[t].val[i][j];
			}
			tempGama.gamma[i] = temp;
			tempGama.alpha[i] = 0;
			tempGama.beta[i] = 0;
			if (t == 0)
			{
				pi[i] = temp;
				//cout <<pi[i]<<endl;
			}
		}
		gama.push_back(tempGama);
	}

	//	New A matrix
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			long double numerator = 0, denominator = 0;
			for (int t = 0; t < length - 1; t++)
			{
				numerator = numerator + zeta[t].val[i][j];
				denominator = denominator + gama[t].gamma[i];
			}
			A[i][j] = numerator / denominator;

		}

	}

	//New B model
	for (int j = 0; j < N; j++)
	{
		for (int k = 0; k < M; k++)
		{
			long double numerator = 0, denominator = 0;

			for (int t = 0; t < length - 1; t++)
			{
				int indx = ObservationSeq[t] - 1;
				if (k == indx)
					numerator = numerator + gama[t].gamma[j];
				denominator = denominator + gama[t].gamma[j];
			}
			B[j][k] = numerator / denominator;
			if (B[j][k] == 0)
				B[j][k] = 10e-80;

		}

	}
	gama.clear();


	for (int i = 0; i <N; i++)
	{
		long double colSum = 0;
		long double tempMax = 0, substract = colSum - 1;

		for (int j = 0; j < M; j++)
		{
			colSum = colSum + B[i][j];
		}


		if (substract > 0)
		{
			int maxIndex = 0;
			for (int j = 0; j < M; j++)
			{
				if (B[i][j] > tempMax)
				{
					tempMax = B[i][j];
					maxIndex = j;
				}
			}
			B[i][maxIndex] = B[i][maxIndex] - substract; 
			
		}



	}
}

void createObservationSequence(ofstream &obsSeq_obj, double cepstral[], int len)
{
	double tokura_wts[12] = { 1, 3, 7, 13, 19, 22, 25, 33, 42, 50, 56, 61 };
	vector <double> distance_from_cluster;
	for (int i = 0; i < codebook.size(); i++)
	{
		double temp_d = 0;
		for (int l = 0; l < len; l++)
		{
			double  temp_prod = (codebook[i].ci[l] - cepstral[l]) *
				(codebook[i].ci[l] - cepstral[l]);
			temp_d = temp_d + temp_prod * tokura_wts[l];
		}
		distance_from_cluster.push_back(temp_d);				// Stores distances of each training word from all centroids	
	}
	int index = distance(distance_from_cluster.begin(), min_element(distance_from_cluster.begin(), distance_from_cluster.end()));
	
	obsSeq_obj << codebook[index].index << " ";
	
	ObservationSeq.push_back(codebook[index].index);
	
	distance_from_cluster.clear();
}



void calculate_coeff(ofstream &obsSeq_obj)
{
	//obsSeq_obj << "going to write";
	// Calculating Ci
	int k = 0;
	vector <double> frame;
	std::fstream cropped_file("cropped_file.txt", std::ios_base::in);
	while (!cropped_file.eof()){
		cropped_file >> temp_input[k];
		k++;
	}
	long int i, m, n;
	double temp, c_i[P];
	int frame_size = window_size * 1 / 4;
	long int fileLen = temp_input.size();
	int count = 0;
	for (int i = 0; i < fileLen - window_size; i = i + frame_size)
	{
		for (n = 0; n<window_size; n++)
		{
			double win_val = 0.54 - 0.46 * cos(2 * 3.141*n / (window_size - 1));  //hamming window
			frame.push_back(win_val * temp_input[i + n]);
		}
		autocorelation(frame); //calculating autocorelation 
		if (R[0] < 100)
		{
			R.clear();
			frame.clear();
			continue;
		}
		double max_R = R[0];
		count = count + 1;
		
		for (n = 0; n < R.size(); n++)
		{
			R[n] = R[n] / max_R;			// Normalizing the autocorrelation values

		}
		calculate_alpha_i(); //calculating alpha using durbins algo
		frame.clear();
		R.clear();
		// calculating Ci

		for (n = 1; n <= P; n++)
		{
			temp = 0;
			for (m = 1; m <= n - 1; m++)
			{
				temp = temp + (m / n)*c_i[m - 1] * alpha[n - m - 1]; 
			}
			c_i[n - 1] = alpha[n - 1] + temp;
			c_i[n - 1] = c_i[n - 1] * (1 + sin(3.141*(n - 1) / (P - 1))*P / 2);
		}
		count = count + 1;
		createObservationSequence(obsSeq_obj, c_i, P); //This assigns the Ci to codebook and makes the observation sequence
		
	}
	
	//cout << "Observation sequence size " << count << endl;
}



void buildHMM(){
	cout << "I appreciate your patience !!!"<<endl;
	cout << "Let me build the model :)" <<endl<<endl;
	codebook.clear();
	temp_input.clear();
	autoCorr.clear();
	R.clear();
	temp_input.clear();
	double temp;
	loadCodebook();
	
	vector <HMMmodel> HMM2;
	HMMmodel modelTemp;
	ostringstream ObservationSeq_string;
	ObservationSeq_string << "./write/observationSeq.txt";
	string ObservationSeq_filename = ObservationSeq_string.str();
	ofstream ObservationSeq_objct;
	ObservationSeq_objct.open(ObservationSeq_filename);

		for (int i = 0; i<4; i++)
		{
			vector <HMMmodel> resultHMM;
			ObservationSeq_objct << "Number " << i << endl;

			for (int j = 1; j <= no_of_training_data; j++)
			{

				cout << "Digit " << i << " from file no " << j << endl;
				ostringstream in_name;
				in_name << "./digit_train/" << i << "/" << j << ".txt";
				string in_filename = in_name.str();
				ifstream infile_obj;
				infile_obj.open(in_filename);
				while (!infile_obj.eof())
				{
					infile_obj >> temp;
					temp_input.push_back(temp);
				}
				infile_obj.close();						
				dc_norm();				// DC shift, Normalization on the signal	
				crop();
				calculate_coeff(ObservationSeq_objct);			// calculate ci													
				ObservationSeq_objct << endl;
				
					for (int n = 0; n<N; n++)
					{
						pi[n] = pi_0[n];
						for (int m = 0; m<N; m++)
							A[n][m] = A_0[n][m];
						for (int m = 0; m<M; m++)
							B[n][m] = B_0[n][m];
					}
				
				
					//cout << "paaluu";

				long double prob_old;
				long int count = 0;
				while (1)
				{
					
					long double prob_new = alphabetasolution();

					zetasolution3();
					ReEstimation();
					alpha_beta.clear();
					zeta.clear();

					if (count == 0)
					{
						prob_old = prob_new;
						count = count + 1;
					}
					else
					{
						//cout << (prob_new - prob_old) * 100 / prob_old << endl;
						//If change is less than the threshold then no need for reestimating
						if ((prob_new - prob_old) * 100 / prob_old < 10e-5)
						{
							break;
						}
						prob_old = prob_new;
						count = count + 1;
					}
					
				}
				temp_input.clear();
				ObservationSeq.clear();
				for (int n = 0; n< N; n++)
				{
					modelTemp.pi[n] = pi[n];

					for (int m = 0; m < N; m++){
						modelTemp.A[n][m] = A[n][m];
						//cout << modelTemp.A[n][m]<<" ";
					}
					//cout << endl<< endl;
					for (int m = 0; m < M; m++){
						modelTemp.B[n][m] = B[n][m];
						
					}
						
				}
				resultHMM.push_back(modelTemp);
			}
			
			for (int n = 0; n<N; n++)
			{
				modelTemp.pi[n] = 0;
				for (int l = 0; l<no_of_training_data; l++)
					modelTemp.pi[n] = modelTemp.pi[n] + resultHMM[l].pi[n];
				modelTemp.pi[n] = modelTemp.pi[n] / no_of_training_data;

				for (int m = 0; m<N; m++)
				{
					modelTemp.A[n][m] = 0;
					for (int l = 0; l<no_of_training_data; l++)
						modelTemp.A[n][m] = modelTemp.A[n][m] + resultHMM[l].A[n][m];
					modelTemp.A[n][m] = modelTemp.A[n][m] / no_of_training_data;
				}
				for (int m = 0; m<M; m++)
				{
					modelTemp.B[n][m] = 0;
					for (int l = 0; l<no_of_training_data; l++)
						modelTemp.B[n][m] = modelTemp.B[n][m] + resultHMM[l].B[n][m];
					modelTemp.B[n][m] = modelTemp.B[n][m] / no_of_training_data;
				}
			}

			HMM2.push_back(modelTemp);
			resultHMM.clear();
		}
	
	cout << " Model building done !!!!!!! Hurraaah " << endl;
	ObservationSeq_objct.close();

	ostringstream out_file_path;
	out_file_path << "./write/finalHMM.txt";
	string finalhmm = out_file_path.str();
	ofstream hmm_objects;
	hmm_objects.open(finalhmm);

	for (int i = 0; i < 4; i++)
	{
		ostringstream out_file_path2;
		out_file_path2 << "./write/HMM_training/" << i << ".txt";
		string trainedhmm = out_file_path2.str();
		ofstream trainedhmm_obj;
		trainedhmm_obj.open(trainedhmm);

		hmm_objects << "Pi Values for " << i << endl;
		for (int n = 0; n< N; n++)
		{
			hmm_objects << scientific << setprecision(20) << HMM2[i].pi[n] << " ";
			trainedhmm_obj << scientific << setprecision(20) << HMM2[i].pi[n] << endl;
		}
		hmm_objects << endl << "A Matrix for " << i << endl;
		for (int n = 0; n < N; n++)
		{
			for (int m = 0; m < N; m++)
			{
				hmm_objects << scientific << setprecision(20) << HMM2[i].A[n][m] << " ";
				trainedhmm_obj << scientific << setprecision(20) << HMM2[i].A[n][m] << endl;
			}
			hmm_objects << endl;
		}
		hmm_objects << endl << "B Matrix for " << i << endl;
		for (int n = 0; n < N; n++)
		{
			for (int m = 0; m < M; m++)
			{
				hmm_objects << scientific << setprecision(20) << HMM2[i].B[n][m] << " ";
				trainedhmm_obj << scientific << setprecision(20) << HMM2[i].B[n][m] << endl;
			}
			hmm_objects << endl;
		}
		hmm_objects << endl;
		trainedhmm_obj.close();
	}
	hmm_objects.close();



	return;
}

void getProbabilities(){

	for (int i = 0; i < N; i++)			// Initialization
	{
		int indx = ObservationSeq[0] - 1;
		test_alpha.push_back(pi[i] * B[i][indx]);
	}

	vector <long double> alpha1;		// Induction			
	int length,k;
	if (ObservationSeq.size() <= T)
		length = ObservationSeq.size();
	else
		length = T+1;
	// Calculation of forward probability
	for ( k = 1; k < length; k++)		
	{
		for (int j = 0; j < N; j++)
		{
			long double temp_a = 0;
			for (int i = 0; i<N; i++)
				temp_a = temp_a + test_alpha[i] * A[i][j];
			int indx = ObservationSeq[k] - 1;
			alpha1.push_back(temp_a * B[j][indx]);
		}
		test_alpha.clear();
		for (int i = 0; i < N; i++)
			test_alpha.push_back(alpha1[i]);
		alpha1.clear();
	}
	// Termination
	long double prob1 = 0;
	int temp = ObservationSeq[0] - 1;
	for (int i = 0; i < N; i++)
		prob1 = prob1 + test_alpha[i];
	test_probabilities.push_back(prob1);
	test_alpha.clear();
}


void compare_trained_hmm(){
	
	for (int i = 0; i < 4; i++)
	{
		long double temp;
		ostringstream training_file;
		training_file << "./write/HMM_training/" << i << ".txt";
		string filename = training_file.str();
		ifstream infile_obj;
		infile_obj.open(filename);
		int indx = 0;
		while (!infile_obj.eof())
		{
			infile_obj >> temp;
			if (indx >= 0 && indx <= 4)	//first five values are PI and next
				pi[indx] = temp;
			if (indx >= 5 && indx <= 29)
				A[(indx - 5) / N][(indx - 5) % N] = temp;
			if (indx >= 30 && indx <= 189)
				B[(indx - 30) / M][(indx - 30) % M] = temp;
			indx = indx + 1;
		}
		
		infile_obj.close();
		getProbabilities();
	}
}

int predictDigit()
{
	int index;
	/*for(int i = 0; i < test_probabilities.size(); i++){
		cout << test_probabilities[i] << endl;
	}*/
	index = distance(test_probabilities.begin(), max_element(test_probabilities.begin(), test_probabilities.end()));
	for (int i = 0; i < test_probabilities.size(); i++){
		cout << test_probabilities[i] << endl;
	}
	return index;

	
	
}

void testing(){
	cout << "Wait I am trying to analyse what you said....................." << endl<<endl;
	codebook.clear();
	ObservationSeq.clear();
	temp_input.clear();
	autoCorr.clear();
	R.clear();
	temp_input.clear();
	alpha_beta.clear();
	double temp;
	ostringstream in_name;
	in_name << "test.txt";
	string in_filename = in_name.str();
	ifstream infile_obj;
	infile_obj.open(in_filename);

	loadCodebook();
	ostringstream ObservationSeqtest_string;
	ObservationSeqtest_string << "./write/observationSeqtest.txt";
	string ObservationSeq_filename = ObservationSeqtest_string.str();
	ofstream ObservationSeq_objct_test;
	ObservationSeq_objct_test.open(ObservationSeq_filename);

	while (!infile_obj.eof())
	{
		infile_obj >> temp;
		temp_input.push_back(temp);
	}
	dc_norm();				// DC shift and normalization
	crop();					//cut the unvoiced part
	infile_obj.close();
	ObservationSeq_objct_test << temp;
	calculate_coeff(ObservationSeq_objct_test);			// calculate Ci for the testing data
	
	compare_trained_hmm();
	int predicted_number = predictDigit();

	cout << "You just said " << predicted_number << "   ?????" << endl<<endl;
	ObservationSeq_objct_test.close();
	temp_input.clear();
	ObservationSeq.clear();
	test_probabilities.clear();
	return;
}

