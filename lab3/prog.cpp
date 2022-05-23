#include <cstdio>
#include <mpi.h>

double* allocateMatrix(int nRows, int nCols){
	double* matrixToReturn = (double*)malloc(nRows * nCols * sizeof(double));
	memset(matrixToReturn, 0.0, nRows * nCols);
	return matrixToReturn;
}

double getRandomDouble(double dMin, double dMax){
    double d = (double)rand() / RAND_MAX;
    return dMin + d * (dMax - dMin);
}

// Matrix is being filled using getRandomDouble function.
// This function uses rand() to generate pseudo-random numbers.
// Since rand() is pseudo-random generator, it is deterministic.
// Hence every call of this function will fill matrix in the same numbers.
void fillMatrix(double* matrix, int rows, int cols, double dMin, double dMax){
	for (int i = 0; i < rows; ++i){
		for (int j = 0; j < cols; ++j){
			double randomDouble = getRandomDouble(dMin, dMax);
			matrix[i * cols + j] = randomDouble;
		}
	}
}

// lMatrix is NxM, rMatrix is MxK, res is NxK
void multiplyMatrices(double* lMatrix, double* rMatrix, double* res, int N, int M, int K){
	for (int i = 0; i < N; ++i){
		for (int j = 0; j < M; ++j){
			for (int k = 0; k < K; ++k){
				res[i*K + k] += lMatrix[i*M + j] * rMatrix[j*K + k]; 
			}
		}
	}
}

int main(int argc, char* argv[]){

	int ndims = 2;
	int xdim = x_dim_number;
	int ydim = y_dim_number;

	// A has dimensions (N1, N2), B has dimensions (N2, N3), C has dimensions (N1, N3);
	// N1 and N3 must be divisible by 
	// gcd(ydim1, ydim2, ..., ydimn)=24 and gcd(xdim1, xdim2, ..., xdimn)=24
	int N1 = 24 * 1357; // ydim * rows_per_proc; 
	int N2 = 333; // every number is acceptable
	int N3 = 24 * 1246; // xdim * cols_per_proc;

	int rows_per_proc = N1 / ydim;
	int cols_per_proc = N3 / xdim;

	MPI_Init(NULL, NULL);

	int world_rank;
	int nProcs;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

	/* A_ROW datatype:
	 * 
	 *    +--------+--------+--------+       +--------+
	 * A: | DOUBLE | DOUBLE | DOUBLE |  ...  | DOUBLE |  ...
	 *    +--------+--------+--------+       +--------+
	 *     \___________  ___________________/ \__  ___________
	 *                 \/                        \/
	 *            First A_ROW                Next A_ROW
	*/
	MPI_Datatype A_ROW;
	MPI_Type_vector(N2, 1, 1, MPI_DOUBLE, &A_ROW);
	MPI_Type_commit(&A_ROW);

	/* B_COL datatype:
	 * 
	 * 2D representation of 1D matrix B:
	 *
	 *
     *             B_COL start   B_COL start   B_COL start   B_COL start           B_COL start
     *                  |             |             |             |                     |
	 *                  V             V             V             V                     V
	 *           +-------------+-------------+-------------+-------------+       +-------------+
	 *         / |           0 |           1 |           2 |           3 |  ...  |        N2-1 |
	 *        |  +-------------+-------------+-------------+-------------+       +-------------+
	 *        |  |      1*N2+0 |      1*N2+1 |      1*N2+2 |      1*N2+3 |  ...  |    2*(N2-1) |
	 *        |  +-------------+-------------+-------------+-------------+       +-------------+
	 *  One   |  |      2*N2+0 |      2*N2+1 |      2*N2+2 |      2*N2+3 |  ...  |    3*(N2-1) |
	 * B_COL {   +-------------+-------------+-------------+-------------+       +-------------+
	 *        |  |      3*N2+0 |      3*N2+1 |      3*N2+2 |      3*N2+3 |  ...  |    4*(N2-1) |
	 *        |  +-------------+-------------+-------------+-------------+       +-------------+
	 *        |        ...           ...           ...           ...         \         ...
	 *        |  +-------------+-------------+-------------+-------------+       +-------------+
	 *         \ | (N1-1)*N2+0 | (N1-1)*N2+1 | (N1-1)*N2+2 | (N1-1)*N2+3 |  ...  |(N1-1)*(N2-1)|
	 *           +-------------+-------------+-------------+-------------+       +-------------+
	 *                  ^             ^             ^             ^                     ^
	 *                  |             |             |             |                     |
	 *              B_COL end     B_COL end     B_COL end     B_COL end             B_COL end
	*/
	MPI_Datatype B_COL_NO_RESIZE, B_COL;
	MPI_Type_vector(N2, 1, N3, MPI_DOUBLE, &B_COL_NO_RESIZE);
	MPI_Type_commit(&B_COL_NO_RESIZE);
	MPI_Type_create_resized(B_COL_NO_RESIZE, 0, 1*sizeof(double), &B_COL);
	MPI_Type_commit(&B_COL);

	/* C_BLOCK datatype
	 *
	 * 2D representation of 1D matrix C:
	 *
	 * 					One block of C_BLOCK_NO_RESIZE      C_BLOCK line ends
	 *  C_BLOCK line        (contains some doubles)             |       
	 *	starts     _________________/\____________________      |
     *     |      /                                       \     |
	 *     |     +-------------+-------------+-------------+----|--------+       +-------------+
	 *   0 +---> |           0 |           1 |           2 | <--+      3 |  ...  |        N3-1 |
	 *     |     +-------------+-------------+-------------+----|--------+       +-------------+
	 *   1 +---> |      1*N3+0 |      1*N3+1 |      1*N3+2 | <--+ 1*N3+3 |  ...  |    2*(N3-1) |
	 *     |     +-------------+-------------+-------------+----|--------+       +-------------+
	 *   2 +---> |      2*N3+0 |      2*N3+1 |      2*N3+2 | <--+ 2*N3+3 |  ...  |    3*(N3-1) |
	 *     |     +-------------+-------------+-------------+----|--------+       +-------------+
	 *   3 +---> |      3*N3+0 |      3*N3+1 |      3*N3+2 | <--+ 3*N3+3 |  ...  |    4*(N3-1) |
	 *     |     +-------------+-------------+-------------+----|--------+       +-------------+
	 *    ...          ...           ...           ...           ...         \         ...
	 *     |     +-------------+-------------+-------------+-------------+       +-------------+
	 *     +---> | (N1-1)*N3+0 | (N1-1)*N3+1 | (N1-1)*N3+2 | (N1-1)*N3+3 |  ...  |(N1-1)*(N3-1)|
	 *           +-------------+-------------+-------------+-------------+       +-------------+
	 *            \________________  _________________________________________________________/    
	 *                             \/
	 *                        One stride
	 *
	 * First C_BLOCK starts at index "0", next starts from index "3" and so on.
	 * The block following the last one in the first row must start from index "rows_per_proc*N3"
	 *
	 * This datatype used only for gathering the C matrix. The principle of collecting the C matrix
	 * is shown below, after matrix multiplication.
	*/

	MPI_Datatype C_BLOCK_NO_RESIZE, C_BLOCK;
	MPI_Type_vector(rows_per_proc, cols_per_proc, N3, MPI_DOUBLE, &C_BLOCK_NO_RESIZE);
	MPI_Type_commit(&C_BLOCK_NO_RESIZE);
	MPI_Type_create_resized(C_BLOCK_NO_RESIZE, 0, cols_per_proc*sizeof(double), &C_BLOCK);
	MPI_Type_commit(&C_BLOCK);

	/* Cartesian grid creation
	 *
	 *               xdim
     *              __/\________________________________________________________________
	 *             /                                                                    \
	 *          +--------------------------------------------------------------------------> X coord
	 *          | +----------+   +----------+   +----------+   +----------+   +----------+
	 *        / | |(0, 0),  0|---|(1, 0),  1|---|(2, 0),  2|---|(3, 0),  3|---|(4, 0),  4|
	 *       |  | +----------+   +----------+   +----------+   +----------+   +----------+
	 * ydim {   |       |              |              |              |              |
	 *       |  | +----------+   +----------+   +----------+   +----------+   +----------+
	 *       |  | |(0, 1),  5|---|(1, 1),  6|---|(2, 1),  7|---|(3, 1),  8|---|(4, 1),  9|
	 *       |  | +----------+   +----------+   +----------+   +----------+   +----------+
	 *       |  |       |              |              |              |              |
	 *       |  | +----------+   +----------+   +----------+   +----------+   +----------+
	 *        \ | |(0, 2), 10|---|(1, 2), 11|---|(2, 2), 12|---|(3, 2), 13|---|(4, 2), 14|
	 *          | +----------+   +----------+   +----------+   +----------+   +----------+
	 *          |
	 *  Y coord V
	 *
	 *	Cell structure:
	 *	+---------------------------------------------+
	 *	| (cart_coords[0], cart_coords[1]), cart_rank |
	 *	+---------------------------------------------+
	*/
	int dims[ndims] = {xdim, ydim};
	int periods[ndims] = {false, false};
	int reorder = true;
	MPI_Comm COMM_CART, COMM_ROW, COMM_COL;
	
	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &COMM_CART);

	/* Get subgrids of COMM_CART
	 * MPI_Cart_sub(..., {true, false}, ...) makes new communicators:
	 *
	 * @==========================================================@
	 * #                                                          #
	 * #  +-----+   +-----+   +-----+   +-----+   +-----+         #         First
	 * #  |   0 |---|   1 |---|   2 |---|   3 |---|   4 |--- ...  #   <---- row
	 * #  +-----+   +-----+   +-----+   +-----+   +-----+         #         communicator
	 * #     |         |         |         |         |            #
	 * @=====|=========|=========|=========|=========|============@
	 *       |         |         |         |         |
	 * @=====|=========|=========|=========|=========|============@
	 * #     |         |         |         |         |            #
	 * #  +-----+   +-----+   +-----+   +-----+   +-----+         #         Second
	 * #  |   0 |---|   1 |---|   2 |---|   3 |---|   4 |--- ...  #   <---- row
	 * #  +-----+   +-----+   +-----+   +-----+   +-----+         #         communicator
	 * #     |         |         |         |         |            #
	 * @=====|=========|=========|=========|=========|============@
	 *       |         |         |         |         |
	 *      ...       ...       ...       ...       ...               <---- and so on
	 *
	 * Cells illustrated with ranks which belong to new communicators.
	 * Every communicator has its own ranking (but not a ranking rule)!
	 * @@@ be careful it's a trap! @@@
	 *
	 * MPI_Cart_sub(..., {false, true}, ...) separates COMM_CART along the X axis
	 * 		and makes new "vertical" communicators as like as in the picture shown above.
	*/
	int save_xaxis[2] = {true, false};
	int save_yaxis[2] = {false, true};
	MPI_Cart_sub(COMM_CART, save_xaxis, &COMM_ROW);
	MPI_Cart_sub(COMM_CART, save_yaxis, &COMM_COL);


	// Get the cart_rank and the coordinates in a grid
	int cart_rank;
	MPI_Comm_rank(COMM_CART, &cart_rank);

	int cart_coords[ndims];
	MPI_Cart_coords(COMM_CART, cart_rank, ndims, cart_coords);

	/* Matrix creating and filling
	 * 
	 * Assume that first dimension of A matrix is divisible by ydim.
	 * Also assume that second dimension of B matrix is divisible by xdim.
	 * 
	 * Matrices will be created and filled only in root process (process with cart_rank == 0),
	 * then parts will be broadcasted to another processes.
	*/
	double* A = NULL;
	double* B = NULL;
	double* C = NULL;

	double* A_part = allocateMatrix(rows_per_proc, N2);
	double* B_part = allocateMatrix(N2, cols_per_proc);
	double* C_part = allocateMatrix(rows_per_proc, cols_per_proc);

	if (cart_rank == 0){
		A = allocateMatrix(N1, N2);
		B = allocateMatrix(N2, N3);
		C = allocateMatrix(N1, N3);
		fillMatrix(A, N1, N2, -100.0, 100.0);
		fillMatrix(B, N2, N3, -100.0, 100.0);
	}

	timeStart = MPI_Wtime();

	// Scatter and broadcast A parts
	if (cart_coords[0] == 0){
		MPI_Scatter(A, rows_per_proc, A_ROW, A_part, rows_per_proc*N2, MPI_DOUBLE, 0, COMM_COL);
	}
	MPI_Bcast(A_part, rows_per_proc*N2, MPI_DOUBLE, 0, COMM_ROW);

	// Scatter and broadcast B parts
	if (cart_coords[1] == 0){
		MPI_Scatter(B, cols_per_proc, B_COL, B_part, N2*cols_per_proc, MPI_DOUBLE, 0, COMM_ROW);
	}
	MPI_Bcast(B_part, N2*cols_per_proc, MPI_DOUBLE, 0, COMM_COL);

	// Now every process has its own part of A matrix and part of B matrix.
	// After mulpiplication every process will have a small block of C matrix.
	multiplyMatrices(A_part, B_part, C_part, rows_per_proc, N2, cols_per_proc);

	
	/* Gathering C matrix
	 *
	 * Every process will send only one block to root process. So, proc_recievecount
	 * must be 1.
	 *
	 * Consider matrix made out of 2x3 blocks. Every cell of this diagram is a single double.
	 *
	 *       0   1   2   3   4   5
	 *   +---------------------------> X coord
	 *   | +---+---+---+---+---+---+
	 * 0 | |   |   |///|///|###|###|
	 *   | +---+---+---+---+---+---+
	 * 1 | |   |   |///|///|###|###|
	 *   | +---+---+---+---+---+---+
	 * 2 | |   |   |///|///|###|###|
	 *   | +---+---+---+---+---+---+
	 * 3 | |@@@|@@@|***|***|||||||||
	 *   | +---+---+---+---+---+---+
	 * 4 | |@@@|@@@|***|***|||||||||
	 *   | +---+---+---+---+---+---+
	 * 5 | |@@@|@@@|***|***|||||||||
	 *   | +---+---+---+---+---+---+
	 *   V
	 *   Y coord
	 *
	 * Pair of 2D coordinates can be mapped to a 1D index using formula:
	 *		Y * 6 + X
	 *
	 * C_BLOCK starting indexes: {0, 2, 4, 18, 20, 22}. 
	 * This indexes can be calculated using formula:
	 * 		Y * rows_per_proc * xdim + X
	 * where rows_per_proc - constant height of a block, 
	 *       xdim - X dimension of a grid.
	*/

	// Calculate proc_recievecount and proc_displacement locally.
	int proc_recievecount = 1;
	int proc_displacement = (cart_coords[1] * rows_per_proc * xdim) + cart_coords[0];
	int rc_C[nProcs];
	int dp_C[nProcs];

	// Gather rc and dp from all processes to local arrays.
	MPI_Allgather(&proc_recievecount, 1, MPI_INT, rc_C, 1, MPI_INT, COMM_CART);
	MPI_Allgather(&proc_displacement, 1, MPI_INT, dp_C, 1, MPI_INT, COMM_CART);
	
	// Gather C matrix
	MPI_Gatherv(C_part, rows_per_proc*cols_per_proc, MPI_DOUBLE, C, rc_C, dp_C, C_BLOCK, 0, COMM_CART);
	
	timeEnd = MPI_Wtime();
	
	// freeing used memory
	free(A);
	free(B);
	free(C);
	free(A_part);
	free(B_part);
	free(C_part);
	
	MPI_Comm_free(&COMM_CART);
	MPI_Comm_free(&COMM_ROW);
	MPI_Comm_free(&COMM_COL);
	MPI_Type_free(&A_ROW);
	MPI_Type_free(&B_COL_NO_RESIZE);
	MPI_Type_free(&B_COL);
	MPI_Type_free(&C_BLOCK_NO_RESIZE);
	MPI_Type_free(&C_BLOCK);
	
	if (cart_rank == 0){
		printf("Program execution took %.2lf seconds\n", timeEnd - timeStart);
	}

	MPI_Finalize();
	return 0;
}