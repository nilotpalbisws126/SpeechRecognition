#include "stdafx.h"
#include "stdafx.h"
#include<iostream>
#include<fstream>
#include <string>
#include<cmath>
#include<Windows.h>
#include "vector"
#include <sstream>
#include <time.h>
#include <iomanip>
#include <algorithm>
#include "feature_ext.h"
#include "make_codebook.h"


using namespace std;

int _tmain(int argc, _TCHAR* argv[])
{
	int opt;
	cout << "****************************************************************" << endl;
	cout << "********************** Digit recognition system ****************" << endl;
	cout << "****************************************************************" << endl;
	cout << "Choose a option" << endl;
	cout << "1. Train the model with your data locate in digit_train folder and build the model (Caution: It takes too much time)" << endl;
	cout << "2. Test the system with live voice that will be recored now and I will try to detect what you said " << endl;
	cin >> opt;
	switch (opt){
	case 1: 
		train();
		lbg();
		buildHMM();
		break;
	case 2: 
		cout << "Now record the Input and Press input \n\n\n";

		system("Recording_Module.exe 2 test.wav test.txt");
		testing();
		break;
	default:
		cout << "You have chosen a wrong option. Choose 1 or 2" << endl;

	}
	
	//*******************************************************************************************************************************

	
	system("pause");
	return 0;
}

