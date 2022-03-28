#pragma once

#include <iostream>
#include <iomanip>
#include <chrono>

void Utility_setRandomSeed();

double Utility_getRandomDouble(double dMin, double dMax);

void Utility_swap(double** a, double** b);

void Utility_printVector(double* vector, int N);

void Utility_printMatrix(double* matrix, int N);