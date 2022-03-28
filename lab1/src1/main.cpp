#include <iostream>
#include <fstream>
#include <sys/time.h>
#include "utility.h"
#include "algebra.h"

int main(){
	//Utility_setRandomSeed();

	std::ifstream file("./data.txt");

	int N = 2300;
	struct timespec start, end;

	clock_gettime(CLOCK_MONOTONIC_RAW, &start);
	double epsilon = 1e-5;

	double* A = Algebra_allocateMatrix(N);
	double* b = Algebra_allocateVector(N);

	for (int i = 0; i < N; ++i){
		for (int j = 0; j < N; ++j){
			file >> A[i * N + j];
		}
	}
	for (int i = 0; i < N; ++i){
		file >> b[i];
	}
	double* x = Algebra_allocateVector(N);

	int nSteps = 0;
	x = Algebra_solveSystem(A, x, b, epsilon, N, &nSteps);

	double error = Algebra_getErrorOfSolving(A, x, b, N);

	std::cout << "Error: " << error << "\nSteps: " << nSteps << std::endl;

	free(A);
	free(x);
	free(b);
	
	clock_gettime(CLOCK_MONOTONIC_RAW, &end);

	int time_ms = (end.tv_sec - start.tv_sec) * 1000 + (end.tv_nsec - start.tv_nsec) / 1000000;
	printf("Program execution took %d.%d seconds\n", time_ms / 1000, time_ms % 1000);

	return 0;
}
