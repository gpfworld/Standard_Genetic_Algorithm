#include <iostream>
#include "genetic.h"
#include "genetic.cpp"

using namespace std;

int main (int argc, char *argv[])
{
//	if (argc == 1){
//		cout << "Error, must give a file name for parameters" << endl;
//		exit(1);
//	}
//	char *file = argv[1];
//	if (argc == 2){
//		cout << "Error, must give a file name for outuput" << endl;
//		exit(1);
//	}
//	char *o_file = argv[2];
//
//	int seed = 100;
//	if (argc == 4)
//		seed = atoi(argv[3]);

	string parameters_path = "/Users/wintic/Documents/gpf/code/clion_code/standard-genetic-algorithm/parameters.txt";
	string run_path = "/Users/wintic/Documents/gpf/code/clion_code/standard-genetic-algorithm/run.txt";
	Genetic *genetic_solver = new Genetic(parameters_path.c_str(),run_path.c_str(), 100);
	genetic_solver->solve();
	//char *bit_string = new char[22];
	//double value;		
	//genetic_solver->bit_encode(-0.999999,-1.000000,2.000000,6,bit_string);
	//cout<<bit_string<<endl;
	//genetic_solver->bit_decode("1011010110101011010111", -1.00000, 2.000000, 6, value);
	//cout<<value<<endl;

	return 0;
}
