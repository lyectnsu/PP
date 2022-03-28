#pragma once

#include <iostream>
#include <iomanip>
#include <chrono>

#ifndef ROOT_PROC
#define ROOT_PROC 0
#endif

void Utility_setRandomSeed();

double Utility_getRandomDouble(double dMin, double dMax);

void Utility_swap(double** a, double** b);

void Utility_printVector(double* vector, int N);

void Utility_printMatrix(double* matrix, int N);

void Utility_initScatter(int* sendCounts, int* sendCountsN, int* displacements, int* displacementsN, int M, int N, int nProcs);