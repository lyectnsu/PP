#include <iostream>
#include <fstream>
#include <mpi.h>

#include "utility.h"
#include "algebra.h"

#ifndef ROOT_PROC
#define ROOT_PROC 0
#endif

int main(int argc, char* argv[]){

	int M = 100;
	int N = 100;
	int MAX_STEPS = 10000;
	double DMIN = -100.0;
	double DMAX =  100.0;
	double EPS = 1e-10; // EPS*EPS

	int rank = 0, nProcs = 0;
	double timeStart, timeEnd;

	MPI_Init(&argc, &argv);
	timeStart = MPI_Wtime();
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int *sendCounts     = (int*)malloc(nProcs * sizeof(int));
	int *sendCountsN    = (int*)malloc(nProcs * sizeof(int));

	int *displacements  = (int*)malloc(nProcs * sizeof(int));
	int *displacementsN = (int*)malloc(nProcs * sizeof(int));

	memset(    sendCounts, 0, nProcs);
	memset(   sendCountsN, 0, nProcs);
	memset( displacements, 0, nProcs);
	memset(displacementsN, 0, nProcs);

	Utility_initScatter(sendCounts, sendCountsN, displacements, displacementsN, M, N, nProcs);

	double* A = NULL;
	double* b = Algebra_allocateVector(N);
	double* x = Algebra_allocateVector(N);

	double* part_A = Algebra_allocateMatrix(sendCounts[rank], N);
	double* part_x = Algebra_allocateVector(sendCounts[rank]);

	if (rank == ROOT_PROC){
		//Utility_setRandomSeed();

		A = Algebra_allocateMatrix(M, N);
		std::ifstream file("./data.txt");
		for (int i = 0; i < N; ++i){
			for (int j = 0; j < N; ++j){
				file >> A[i * N + j];
			}
		}
		for (int i = 0; i < N; ++i){
			file >> b[i];
		}
	}

	MPI_Scatterv(A, sendCountsN, displacementsN, MPI_DOUBLE, part_A, sendCountsN[rank], MPI_DOUBLE, ROOT_PROC, MPI_COMM_WORLD);
	MPI_Bcast(x, N, MPI_DOUBLE, ROOT_PROC, MPI_COMM_WORLD);
	MPI_Bcast(b, N, MPI_DOUBLE, ROOT_PROC, MPI_COMM_WORLD);
	
	// solving system of equations
	int nSteps = 0;
	part_x = Algebra_solveSystem(
		part_A, x, b,
		M, N,
		sendCounts, sendCountsN, displacements, displacementsN,
		EPS, MAX_STEPS,
		rank,
		&nSteps
	);

	MPI_Allgatherv(part_x, sendCounts[rank], MPI_DOUBLE, x, sendCounts, displacements, MPI_DOUBLE, MPI_COMM_WORLD);

	// calculating error
	double error = 0, part_error = 0;
	part_error = Algebra_getErrorOfSolving(
		part_A, x, b,
		M, N,
		sendCounts, sendCountsN, displacements, displacementsN,
		rank
	);

	MPI_Allreduce(&part_error, &error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	if (rank == ROOT_PROC){
		error = std::sqrt(error);
		printf("Error: %lf\nSteps taken: %d\n", error, nSteps);
	}

	free(A);
	free(x);
	free(b);
	free(part_A);
	free(part_x);
	free(sendCounts);
	free(sendCountsN);
	free(displacements);
	free(displacementsN);

	timeEnd = MPI_Wtime();

	if (rank == ROOT_PROC){
		printf("Program execution took %.2lf seconds\n", timeEnd - timeStart);
	}

	MPI_Finalize();

	return 0;
}
