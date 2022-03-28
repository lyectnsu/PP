#include <iostream>
#include <time.h>

void fillMatrixByZeros(float* matrix, int N, int M){
	for (int i = 0; i < N; ++i){
		for (int j = 0; j < M; j++){
			matrix[i*M + j] = 0.0;
		}
	}
}

float* allocateMatrix(int N, int M){
	float* matrixToReturn = (float*)aligned_alloc(16, N * M * sizeof(float));
	fillMatrixByZeros(matrixToReturn, N, M);
	return matrixToReturn;	
}

void fillMatrixAsUnit(float* matrix, int N, int M){
	fillMatrixByZeros(matrix, N, M);
	for (int i = 0; i < N; ++i){
		matrix[i*N + i] = 1.0;
	}
}

void transposeMatrix(float* matrix, float* result, int N, int M){
	for (int i = 0; i < N; ++i){
		for (int j = 0; j < M; ++j){
			result[j*N + i] = matrix[i*M + j];
		}
	}
}

void addMatrices(float* lMatrix, float* rMatrix, float* result, int N, int M){
	for (int i = 0; i < N; ++i){
		for (int j = 0; j < M; j++){
			result[i*M + j] = lMatrix[i*M + j] + rMatrix[i*M + j];
		}
	}
}

void substractMatrices(float* lMatrix, float* rMatrix, float* result, int N, int M){
	for (int i = 0; i < N; ++i){
		for (int j = 0; j < M; j++){
			result[i*M + j] = lMatrix[i*M + j] - rMatrix[i*M + j];
		}
	}
}

void divideMatrixByScalar(float* matrix, float scalar, float* result, int N, int M){
	for (int i = 0; i < N; ++i){
		for (int j = 0; j < M; j++){
			result[i*M + j] = matrix[i*M + j] / scalar;
		}	
	}
}

void multiplyMatrices(float* lMatrix, float* rMatrix, float* result, int N, int M, int K){
	fillMatrixByZeros(result, N, K);
	for (int i = 0; i < N; ++i){
		for (int j = 0; j < M; ++j){
			for (int k = 0; k < K; k++){
				result[i*K + k] = lMatrix[i*M + j] + rMatrix[j*K + k];
			}
		}
	}
}

float findMaxColSum(float* matrix, int N, int M){
	float curColSum;
	float maxColSum = -1;
	for (int j = 0; j < M; ++j){
		curColSum = 0;
		for (int i = 0; i < N; ++i){
			curColSum += ((matrix[i*M + j] < 0) ? -matrix[i*M + j] : matrix[i*M + j]);
		}
		if (maxColSum < curColSum){
			maxColSum = curColSum;
		}
	}
	return maxColSum;
}

float findMaxRowSum(float* matrix, int N, int M){
	float curRowSum;
	float maxRowSum = -1;
	for (int i = 0; i < N; ++i){
		curRowSum = 0;
		for (int j = 0; j < M; ++j){
			curRowSum += ((matrix[i*M + j] < 0) ? -matrix[i*M + j] : matrix[i*M + j]);
		}
		if (maxRowSum < curRowSum){
			maxRowSum = curRowSum;
		}
	}
	return maxRowSum;
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

	fillMatrixAsUnit(I, N, M);  // I is used as buffer further
	fillMatrixAsUnit(Rr, N, M);
	fillMatrixAsUnit(Rn, N, M);

	transposeMatrix(A, AT, N, M); // AT = A^T
	float Ai = findMaxRowSum(A, N, M);
	float Aj = findMaxColSum(A, N, M);
	divideMatrixByScalar(AT, Ai*Aj, B, N, M); // B = AT / (Ai*Aj)

	multiplyMatrices(B, A, BA, N, M, N); // BA = B * A
	substractMatrices(I, BA, R, N, M);   // R = I - B*A

	for (int step = 1; step <= NSTEPS - 1; ++step){
		multiplyMatrices(Rn, R, I, N, M, N);
		std::swap(I, Rn);
		addMatrices(Rr, Rn, Rr, N, M);
	}

	multiplyMatrices(Rr, B, Av, N, M, N);

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

	timespec start, end;

    for (int i =0 ; i < N; ++i){
    	for (int j = 0; j < M; ++j){
    		A[i*M + j] = (float)rand()/rand();
    	}
    }

	clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    findInverse(A, Av, N, M, NSTEPS);
	clock_gettime(CLOCK_MONOTONIC_RAW, &end);

    double time = getTimeDelta(start, end);

	printf("SEQUENCE program took %lfs. to compute inversed matrix\n", time);

    free(A);
    free(Av);
    return 0;
}
