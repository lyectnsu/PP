#include "algebra.h"

double* Algebra_allocateMatrix(int N){
	double* matrixToReturn = NULL;
	while (matrixToReturn == NULL){
		matrixToReturn = (double*)malloc(N * N * sizeof(double));
	}
	for (int i = 0; i < N; ++i){
		for (int j = 0; j < N; ++j){
			matrixToReturn[i*N+j] = 0.0;
		}
	}

	return matrixToReturn;
}

double* Algebra_allocateVector(int N){
	double* vectorToReturn = NULL;
	while (vectorToReturn == NULL){
		vectorToReturn = (double*)malloc(N * sizeof(double));
	}
	for (int i = 0; i < N; ++i){
		vectorToReturn[i] = 0.0;
	}
	return vectorToReturn;
}

// filling matrix symmetrically with random doubles
// and adding weigths to the principal diagonal
//
// should not be optimized
//
void Algebra_fillMatrix(double* matrix, int N, double dMin, double dMax){
	for (int i = 0; i < N; ++i){
		for (int j = i; j < N; ++j){
			if (i == j){
				matrix[i * N + j] += dMax * 2;
			}
			double randomDouble = Utility_getRandomDouble(dMin, dMax);
			matrix[i * N + j] = randomDouble;
			matrix[j * N + i] = randomDouble;
		}
	}
}

// filling vector with random doubles
//
// should not be optimized
//
void Algebra_fillVector(double* vector, int N, double dMin, double dMax){
	for (int i = 0; i < N; ++i){
		vector[i] = Utility_getRandomDouble(dMin, dMax);
	}
}

void Algebra_multiplyMatrixOnVector(double* matrix, double* vector, double* result, int N){
	for (int i = 0; i < N; ++i){
		result[i] = 0;
		for (int j = 0; j < N; ++j){
			result[i] += matrix[i*N + j] * vector[j];
		}
	}
}

void Algebra_multiplyScalarOnVector(double scalar, double* vector, double* result, int N){
	for (int i = 0; i < N; ++i){
		result[i] = scalar * vector[i];
	}
}

double Algebra_scalarProduct(double* vec1, double* vec2, int N){
	double result = 0;
	for (int i = 0; i < N; ++i){
		result += vec1[i] * vec2[i];
	}
	return result;
}

void Algebra_sum(double* vec1, double* vec2, double* result, int N){
	for (int i = 0; i < N; ++i){
		result[i] = vec1[i] + vec2[i];
	}
}

void Algebra_difference(double* vec1, double* vec2, double* result, int N){
	for (int i = 0; i < N; ++i){
		result[i] = vec1[i] - vec2[i];
	}
}

double Algebra_norm(double* vector, int N){
	return std::sqrt(Algebra_scalarProduct(vector, vector, N));
}

void Algebra_copyVector(double* src, double* dest, int N){
	for (int i = 0; i < N; ++i){
		dest[i] = src[i];
	}
}

double* Algebra_solveSystem(double* A, double* x, double* b, double epsilon, int N, int* nSteps){
	double rr;
	double aNext, bNext;
	double* zPrev = Algebra_allocateVector(N);
	double* zNext = Algebra_allocateVector(N);
	double* xPrev = Algebra_allocateVector(N);
	double* xNext = Algebra_allocateVector(N);
	double* rPrev = Algebra_allocateVector(N);	
	double* rNext = Algebra_allocateVector(N);
	double* Az    = Algebra_allocateVector(N);
	double* az    = Algebra_allocateVector(N);
	double* bz    = Algebra_allocateVector(N);
	double* aAz   = Algebra_allocateVector(N);

	Algebra_copyVector(x, xPrev, N);

	// starting point
	// Az is used as a temp container. It should be Ax, but it will cause
	// redundant memory allocation
	Algebra_multiplyMatrixOnVector(A, xPrev, Az, N);   // Ax_0
	Algebra_difference(b, Az, rPrev, N);               // r_0 = b - Ax_0
	Algebra_copyVector(rPrev, zPrev, N);               // z_0 = r_0

	double bNorm = Algebra_norm(b, N);

	int streak = 0;
	while (streak != 5){
		double error = Algebra_norm(rPrev, N) / bNorm;
		//printf("%.7lf\n", error);
		if (error < epsilon){

			if (error == 0){
				break;
			}
			streak++;
		}
		else {
			streak = 0;
		}
		Algebra_multiplyMatrixOnVector(A, zPrev, Az, N);   // Az = A * z_n
		rr = Algebra_scalarProduct(rPrev, rPrev, N);       // (r_n, r_n)

		// a_{n+1} = (r_n, r_n) / (Az_n, z_n)
		aNext = rr / Algebra_scalarProduct(Az, zPrev, N);

		// x_{n+1} = x_n + a_{n+1}z_n
		Algebra_multiplyScalarOnVector(aNext, zPrev, az, N);
		Algebra_sum(xPrev, az, xNext, N);

		// r_{n+1} = r_n - a_{n+1}Az_n
		Algebra_multiplyScalarOnVector(aNext, Az, aAz, N);
		Algebra_difference(rPrev, aAz, rNext, N);

		// b_{n+1} = (r_{n+1}, r_{n+1}) / (r_n, r_n);
		bNext = Algebra_scalarProduct(rNext, rNext, N) / rr;

		// z_{n+1} = r_{n+1} + b_{n+1}z_n
		Algebra_multiplyScalarOnVector(bNext, zPrev, bz, N);
		Algebra_sum(rNext, bz, zNext, N);

		Utility_swap(&xPrev, &xNext);
		Utility_swap(&rPrev, &rNext);
		Utility_swap(&zPrev, &zNext);
		(*nSteps)++;
	}

	free(zPrev);
	free(zNext);
	free(xNext);
	free(Az);
	free(az);
	free(bz);
	free(aAz);
	free(rPrev);	
	free(rNext);
	free(x);

	return xPrev;
}

double Algebra_getErrorOfSolving(double* A, double* x, double* b, int N){
	double* Ax = Algebra_allocateVector(N);
	Algebra_multiplyMatrixOnVector(A, x, Ax, N);

	Algebra_difference(Ax, b, Ax, N);
	double error = Algebra_norm(Ax, N);

	free(Ax);
	return error;
}