#pragma once

#include <mpi.h>
#include <cmath>
#include <cstring>

#include "utility.h"

#ifndef ROOT_PROC
#define ROOT_PROC 0
#endif

double* Algebra_allocateMatrix(int M, int N);

double* Algebra_allocateVector(int N);

void Algebra_fillMatrix(double* matrix, int M, int N, double dMin, double dMax);

void Algebra_fillVector(double* vector, int N, double dMin, double dMax);

void Algebra_multiplyMatrixOnVector(double* matrix, double* vector, double* result, int M, int N);

void Algebra_multiplyScalarOnVector(double scalar, double* vector, double* result, int N);

double Algebra_scalarProduct(double* vec1, double* vec2, int N);

void Algebra_sum(double* vec1, double* vec2, double* result, int N);

void Algebra_difference(double* vec1, double* vec2, double* result, int N);

double Algebra_norm(double* vector, int N);

void Algebra_copyVector(double* src, double* dest, int N);

double* Algebra_solveSystem(
	double* part_A, double* part_x, double* part_b,
	int M, int N,
	int* sendCounts, int* sendCountsN, int* displacements, int* displacementsN,
	double EPS, int MAX_STEPS,
	int rank,
	int* nSteps
);

double Algebra_getErrorOfSolving(
	double* part_A, double* x, double* b,
	int M, int N,
	int* sendCounts, int* sendCountsN, int* displacements, int* displacementsN,
	int rank
);