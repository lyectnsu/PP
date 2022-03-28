#include "utility.h"

void Utility_setRandomSeed(){
	auto nowTime = std::chrono::system_clock::now().time_since_epoch();
    auto timeSinceEpochSec = std::chrono::duration_cast<std::chrono::seconds>(nowTime).count();
	srand(timeSinceEpochSec);
}

double Utility_getRandomDouble(double dMin, double dMax){
    double d = (double)rand() / RAND_MAX;
    return dMin + d * (dMax - dMin);
}

void Utility_swap(double** a, double** b){
    double* tmp = *a;
    *a = *b;
    *b = tmp;
}

void Utility_printVector(double* vector, int N){
    std::cerr << "[ ";
    for (int i = 0; i < N; ++i){
        std::cerr << vector[i] << " ";
    }
    std::cerr << "]" << std::endl;
}

void Utility_printMatrix(double* matrix, int N){
    
    for (int i = 0; i < N; ++i){
        std::cerr << "| ";
        for (int j = 0; j < N; ++j){
            std::cerr << std::setw(5) << matrix[i*N + j] << " ";
        }
        std::cerr << "|" << std::endl;
    }
}