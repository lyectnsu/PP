#pragma once

#include <cmath>

#include "utility.h"

double* Algebra_allocateMatrix(int N);

double* Algebra_allocateVector(int N);

void Algebra_fillMatrix(double* matrix, int N, double dMin, double dMax);

void Algebra_fillVector(double* vector, int N, double dMin, double dMax);

void Algebra_multiplyMatrixOnVector(double* matrix, double* vector, double* result, int N);

void Algebra_multiplyScalarOnVector(double scalar, double* vector, double* result, int N);

double Algebra_scalarProduct(double* vec1, double* vec2, int N);

void Algebra_sum(double* vec1, double* vec2, double* result, int N);

void Algebra_difference(double* vec1, double* vec2, double* result, int N);

double Algebra_norm(double* vector, int N);

void Algebra_copyVector(double* src, double* dest, int N);

double* Algebra_solveSystem(double* A, double* xPrev, double* b, double epsilon, int N, int* nSteps);

double Algebra_getErrorOfSolving(double* A, double* x, double* b, int N);