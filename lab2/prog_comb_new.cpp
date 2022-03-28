#include <iostream>
#include <omp.h>

void fillMatrixByZerosSequence(float* matrix, int N, int M){
	for (int i = 0; i < N; ++i){
		for (int j = 0; j < M; j += 4){
			asm volatile (
				"xorps %%xmm0, %%xmm0;"
				"movaps %%xmm0, (%0);"
				:
				: "r"(matrix + (i*M + j))
				: "%xmm0", "memory"
			);
		}
	}
}

void fillMatrixByZeros(float* matrix, int N, int M){
	#pragma omp for nowait collapse(2)
	for (int i = 0; i < N; ++i){
		for (int j = 0; j < M; j += 4){
			asm volatile (
				"xorps %%xmm0, %%xmm0;"
				"movaps %%xmm0, (%0);"
				:
				: "r"(matrix + (i*M + j))
				: "%xmm0", "memory"
			);
		}
	}
}

float* allocateMatrix(int N, int M){
	float* matrixToReturn = (float*)aligned_alloc(16, N * M * sizeof(float));
	fillMatrixByZerosSequence(matrixToReturn, N, M);
	return matrixToReturn;	
}

void fillMatrixAsUnit(float* matrix, int N, int M){
	fillMatrixByZeros(matrix, N, M);

	#pragma omp for nowait
	for (int i = 0; i < N; ++i){
		matrix[i*N + i] = 1.0;
	}
}

void transposeMatrix(float* matrix, float* result, int N, int M){
	#pragma omp for nowait collapse(2)
	for (int i = 0; i < N; ++i){
		for (int j = 0; j < M; ++j){
			result[j * N + i] = matrix[i * M + j];
		}
	}
}

void addMatrices(float* lMatrix, float* rMatrix, float* result, int N, int M){
	#pragma omp for nowait collapse(2)
	for (int i = 0; i < N; ++i){
		for (int j = 0; j < M; j += 4){
			asm volatile (
				"movaps (%0), %%xmm0;"
				"movaps (%1), %%xmm1;"
				"addps %%xmm0, %%xmm1;"
				"movaps %%xmm1, (%2);"
				:
				: "r"(lMatrix + (i*M + j)), "r"(rMatrix + (i*M + j)), "r"(result + (i*M + j))
				: "%xmm0", "%xmm1", "memory"
			);
		}
	}
}

void substractMatrices(float* lMatrix, float* rMatrix, float* result, int N, int M){
	#pragma omp for nowait collapse(2)
	for (int i = 0; i < N; ++i){
		for (int j = 0; j < M; j += 4){
			asm volatile (
				"movaps (%0), %%xmm0;"
				"movaps (%1), %%xmm1;"
				"subps %%xmm1, %%xmm0;"
				"movaps %%xmm0, (%2);"
				:
				: "r"(lMatrix + (i*M + j)), "r"(rMatrix + (i*M + j)), "r"(result + (i*M + j))
				: "%xmm0", "%xmm1", "memory"
			);
		}
	}
}

void divideMatrixByScalar(float* matrix, float scalar, float* result, int N, int M){
	#pragma omp for nowait collapse(2)
	for (int i = 0; i < N; ++i){
		for (int j = 0; j < M; j += 4){
			asm volatile (
				"movaps (%0), %%xmm0;"
				"movss %1, %%xmm1;"
				"shufps $0x00, %%xmm1, %%xmm1;"
				"divps %%xmm1, %%xmm0;"
				"movaps %%xmm0, (%2);"
				:
				: "r"(matrix + (i*M + j)), "m"(scalar), "r"(result + (i*M + j))
				: "%xmm0", "%xmm1", "memory"
			);
		}	
	}
}

void multiplyMatrices(float* lMatrix, float* rMatrix, float* result, int N, int M, int K){
	fillMatrixByZeros(result, N, K);
	#pragma omp for collapse(3)
	for (int i = 0; i < N; ++i){
		for (int j = 0; j < M; ++j){
			for (int k = 0; k < K; k += 4){
				asm volatile (
					"movaps (%2), %%xmm2;"
					"movss (%0), %%xmm0;"
					"shufps $0x00, %%xmm0, %%xmm0;"
					"movaps (%1), %%xmm1;"
					"mulps %%xmm0, %%xmm1;"
					"addps %%xmm1, %%xmm2;"
					"movaps %%xmm2, (%2);"
					: 
					: "r"(lMatrix + (i*M + j)), "r"(rMatrix + (j*K + k)), "r"(result + (i*K + k))
					: "%xmm0", "%xmm1", "%xmm2", "memory"
				);
			}
		}
	}
}

void findMaxColSum(float* matrix, int N, int M, float& maxColSum){
	#pragma omp for reduction(max:maxColSum)
	for (int j = 0; j < M; ++j){
		float curColSum = 0;
		for (int i = 0; i < N; ++i){
			curColSum += ((matrix[i*M + j] < 0) ? -matrix[i*M + j] : matrix[i*M + j]);
		}
		if (curColSum > maxColSum){
			maxColSum = curColSum;
		}
	}
}

void findMaxRowSum(float* matrix, int N, int M, float& maxRowSum){
	#pragma omp for reduction(max:maxRowSum)
	for (int j = 0; j < M; ++j){
		float curRowSum = 0;
		for (int i = 0; i < N; ++i){
			curRowSum += ((matrix[i*M + j] < 0) ? -matrix[i*M + j] : matrix[i*M + j]);
		}
		if (curRowSum > maxRowSum){
			maxRowSum = curRowSum;
		}
	}
}

double getTimeDelta(timespec &start, timespec &end){
    return (end.tv_sec - start.tv_sec) + 0.000000001*(end.tv_nsec - start.tv_nsec);
}

void findInverse(float* A, float* Av, int N, int M, int NSTEPS){
	float* I  = allocateMatrix(N, M);
	float* B  = allocateMatrix(N, M);
	float* R  = allocateMatrix(N, M);
	float* AT = allocateMatrix(N, M);
	float* BA = allocateMatrix(N, M);
	float* Rr = allocateMatrix(N, M);
	float* Rn = allocateMatrix(N, M);
	float Ai = -1, Aj = -1;

	#pragma omp parallel
	{
		fillMatrixAsUnit(I, N, M);  // I is used as buffer further
		fillMatrixAsUnit(Rr, N, M);
		fillMatrixAsUnit(Rn, N, M);

		transposeMatrix(A, AT, N, M); // AT = A^T
		findMaxRowSum(A, N, M, Ai);
		findMaxColSum(A, N, M, Aj);

		divideMatrixByScalar(AT, Ai*Aj, B, N, M); // B = AT / (Ai*Aj)
		multiplyMatrices(B, A, BA, N, M, N); // BA = B * A
		substractMatrices(I, BA, R, N, M);   // R = I - B*A
	
		for (int step = 1; step <= NSTEPS - 1; ++step){
				multiplyMatrices(Rn, R, I, N, M, N);

				#pragma omp single
				std::swap(I, Rn);

				addMatrices(Rr, Rn, Rr, N, M);
		}
		multiplyMatrices(Rr, B, Av, N, M, N);
	}

	free(I );
	free(B );
	free(R );
	free(AT);
	free(BA);
	free(Rr);
	free(Rn);
}

void printMatrix(float* matrix, int N, int M){

	int N_old = N;
	int M_old = M;
	int MAX_ROW_SHOW = 5;
	int MAX_COL_SHOW = 5;

	N = (N > MAX_ROW_SHOW) ? MAX_ROW_SHOW : N;
	M = (M > MAX_COL_SHOW) ? MAX_COL_SHOW : M;
	for (int i = 0; i < N; ++i){
		for (int j = 0; j < M; ++j){
			printf("%+.3f ", matrix[i*M_old + j]);
		}
		if (N == MAX_ROW_SHOW){
			printf(". . .");
		}
		printf("\n");
	}
	
	if (M == MAX_COL_SHOW){
		for (int i = 0; i < 3; ++i){
			for (int j = 0; j < N; ++j){
				printf("   .   ");
			}
			for (int j = 0; j < 2 * i; ++j){
				printf(" ");
			}
			printf(".\n");
		}
	}
	
	printf("\n");
}


int main(){

	int N = 1000;
	int M = 1000;
	int NSTEPS = 14;

	float* A   = allocateMatrix(N, M);
	float* Av  = allocateMatrix(N, M);

	for (int i =0 ; i < N; ++i){
		for (int j = 0; j < M; ++j){
			A[i*M + j] = (float)rand()/rand();
		}
	}

	double omp_start = omp_get_wtime();
	findInverse(A, Av, N, M, NSTEPS);
	double omp_end = omp_get_wtime();

	/*
	float* B = allocateMatrix(N, M);
	multiplyMatrices(A, Av, B, N, M, N);
	printMatrix(B, N, M);
	*/
	printf("SIMD + OpenMP program took %lfs. to compute inversed matrix\n", omp_end - omp_start);

	//free(B);
	free(A);
	free(Av);
	return 0;
}
