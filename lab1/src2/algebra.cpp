#include "algebra.h"

double* Algebra_allocateMatrix(int M, int N){
	double* matrixToReturn = NULL;
	while (matrixToReturn == NULL){
		matrixToReturn = (double*)malloc(M * N * sizeof(double));
	}
	memset(matrixToReturn, 0.0, M * N);

	return matrixToReturn;
}

double* Algebra_allocateVector(int N){
	double* vectorToReturn = NULL;
	while (vectorToReturn == NULL){
		vectorToReturn = (double*)malloc(N * sizeof(double));
	}
	memset(vectorToReturn, 0.0, N);

	return vectorToReturn;
}

// filling matrix symmetrically with random doubles
// and adding weigths to the principal diagonal
//
// should not be optimized
//
void Algebra_fillMatrix(double* matrix, int M, int N, double dMin, double dMax){
	for (int i = 0; i < M; ++i){
		for (int j = i; j < N; ++j){
			if (i == j){
				matrix[i * N + j] += dMax * 2;
			}
			double randomDouble = Utility_getRandomDouble(dMin, dMax);
			matrix[i * N + j] = randomDouble;
			matrix[j * M + i] = randomDouble;
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

void Algebra_multiplyMatrixOnVector(double* matrix, double* vector, double* result, int M, int N){
	for (int i = 0; i < M; ++i){
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

double* Algebra_solveSystem(
	double* part_A, double* x, double* b,
	int M, int N,
	int* sendCounts, int* sendCountsN, int* displacements, int* displacementsN,
	double EPS, int MAX_STEPS,
	int rank,
	int* nSteps
){

	double aNext, bNext;
	double rr, Azz, rNrN, bNorm;
	double part_rr, part_Azz, part_rNrN, part_bNorm;

	double* zPrev      = Algebra_allocateVector(N);

	double* part_zPrev = Algebra_allocateVector(sendCounts[rank]);
	double* part_zNext = Algebra_allocateVector(sendCounts[rank]);
	double* part_xPrev = Algebra_allocateVector(sendCounts[rank]);
	double* part_xNext = Algebra_allocateVector(sendCounts[rank]);
	double* part_rPrev = Algebra_allocateVector(sendCounts[rank]);
	double* part_rNext = Algebra_allocateVector(sendCounts[rank]);

	double* part_Az    = Algebra_allocateVector(sendCounts[rank]);
	double* part_aAz   = Algebra_allocateVector(sendCounts[rank]);
	double* part_b     = Algebra_allocateVector(sendCounts[rank]);
	double* part_bz    = Algebra_allocateVector(sendCounts[rank]);
	double* part_az    = Algebra_allocateVector(sendCounts[rank]);
	
	MPI_Scatterv(x, sendCounts, displacements, MPI_DOUBLE, part_xPrev, sendCounts[rank], MPI_DOUBLE, ROOT_PROC, MPI_COMM_WORLD);
	MPI_Scatterv(b, sendCounts, displacements, MPI_DOUBLE,     part_b, sendCounts[rank], MPI_DOUBLE, ROOT_PROC, MPI_COMM_WORLD);

	// part_Ax = part_A * x
	Algebra_multiplyMatrixOnVector(part_A, x, part_Az, sendCounts[rank], N);
	// part_rPrev = part_b - part_Ax
	Algebra_difference(part_b, part_Az, part_rPrev, sendCounts[rank]);

	// part_zPrev = part_rPrev
	Algebra_copyVector(part_rPrev, part_zPrev, sendCounts[rank]);

	// gather parts of zPrev together
	MPI_Allgatherv(part_zPrev, sendCounts[rank], MPI_DOUBLE, zPrev, sendCounts, displacements, MPI_DOUBLE, MPI_COMM_WORLD);
	
	// calculate bNorm
	part_bNorm = Algebra_scalarProduct(part_b, part_b, sendCounts[rank]);
	MPI_Allreduce(&part_bNorm, &bNorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	// rr = (rPrev, rPrev)
	part_rr = Algebra_scalarProduct(part_rPrev, part_rPrev, sendCounts[rank]);
	MPI_Allreduce(&part_rr, &rr, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	int stepCount = 0;
	int streak = 0;

	// All we need to have as a whole at each step is a vector zPrev
	// Another vectors such as xPrev, rPrev should be calculated in parts
	while (streak != 5 && stepCount < MAX_STEPS){
		double error = 0, part_error = 0;
		
		part_error = Algebra_scalarProduct(part_rPrev, part_rPrev, sendCounts[rank]);
		MPI_Allreduce(&part_error, &error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		error = error / bNorm;

		if (error < EPS){
			if (error == 0){
				break;
			}
			streak++;
		}
		else {
			streak = 0;
		}

		// part_Az = part_A * zPrev
		Algebra_multiplyMatrixOnVector(part_A, zPrev, part_Az, sendCounts[rank], N);
		// Azz = (A * zPrev, zPrev)
		part_Azz = Algebra_scalarProduct(part_Az, part_zPrev, sendCounts[rank]);
		MPI_Allreduce(&part_Azz, &Azz, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);		

		aNext = rr / Azz;

		// part_az = aNext * part_zPrev
		Algebra_multiplyScalarOnVector(aNext, part_zPrev, part_az, sendCounts[rank]);
		// part_xNext = part_xPrev + part_az
		Algebra_sum(part_xPrev, part_az, part_xNext, sendCounts[rank]);

		// part_aAz = aNext * part_Az
		Algebra_multiplyScalarOnVector(aNext, part_Az, part_aAz, sendCounts[rank]);
		// part_rNext = part_rPrev - part_aAz
		Algebra_difference(part_rPrev, part_aAz, part_rNext, sendCounts[rank]);

		// rNrN = (rNext, rNext)
		part_rNrN = Algebra_scalarProduct(part_rNext, part_rNext, sendCounts[rank]);
		MPI_Allreduce(&part_rNrN, &rNrN, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		
		bNext = rNrN / rr;

		// part_bz = bNExt * part_zPrev
		Algebra_multiplyScalarOnVector(bNext, part_zPrev, part_bz, sendCounts[rank]);
		Algebra_sum(part_rNext, part_bz, part_zNext, sendCounts[rank]);

		// make zPrev, xPrev, rPrev point to zNext, xNext, rPrev
		std::swap(part_zPrev, part_zNext);
		std::swap(part_xPrev, part_xNext);
		std::swap(part_rPrev, part_rNext);
		std::swap(rr, rNrN);

		// gather parts of zPrev together and give it to each process
		MPI_Allgatherv(part_zPrev, sendCounts[rank], MPI_DOUBLE, zPrev, sendCounts, displacements, MPI_DOUBLE, MPI_COMM_WORLD);
		++stepCount;
	}
	
	// resulting number of steps is max stepCount among all processes
	MPI_Allreduce(&stepCount, nSteps, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

	// gather parts of xPrev (the result) together and give it to each process

	free(zPrev);

	free(part_zPrev);
	free(part_zNext);
	free(part_xNext);
	free(part_rPrev);
	free(part_rNext);

	free(part_Az);
	free(part_aAz);
	free(part_b);
	free(part_bz);
	free(part_az);

	return part_xPrev;
}

double Algebra_getErrorOfSolving(
	double* part_A, double* x, double* b,
	int M, int N,
	int* sendCounts, int* sendCountsN, int* displacements, int* displacementsN,
	int rank
){

	double* part_Ax  = Algebra_allocateVector(sendCounts[rank]);
	double* part_b   = Algebra_allocateVector(sendCounts[rank]);
	double* part_bAx = Algebra_allocateVector(sendCounts[rank]);

	MPI_Scatterv(b,  sendCounts,  displacements, MPI_DOUBLE, part_b,  sendCounts[rank], MPI_DOUBLE, ROOT_PROC, MPI_COMM_WORLD);

	Algebra_multiplyMatrixOnVector(part_A, x, part_Ax, sendCounts[rank], N);
	Algebra_difference(part_Ax, part_b, part_bAx, sendCounts[rank]);

	double part_error = Algebra_scalarProduct(part_bAx, part_bAx, sendCounts[rank]);
	
	free(part_Ax);
	free(part_b);
	free(part_bAx);
	return part_error;
}