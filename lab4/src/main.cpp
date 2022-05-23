#include <cstdio>
#include <mpi.h>

#include "field.h"
#include "procdata.h"
#include "settings.h"
#include "utility.h"

void findCycle(Field** part_fields, char* part_cycle_found_arr, int cur_step){
	for (int i = 0; i < cur_step; ++i){
		part_cycle_found_arr[i] = 0;
		if (Field_equal(part_fields[i], part_fields[cur_step])){
			part_cycle_found_arr[i] = 1;
		}
	}
}

int main(int argc, char* argv[]){

	MPI_Init(&argc, &argv);

	ProcData* proc_data = ProcData_create();

	ProcData_init(proc_data);

	Field* starting_field = NULL;

	if (ProcData_getRank(proc_data) == ROOT_PROC){
		starting_field = Field_create(N, M);
		//Field_placeGlider(starting_field, 0, 0);
		Field_fillRandom(starting_field);
	}

	MPI_Scatterv(
		Field_getDataPointer(starting_field, 0), 
		ProcData_getMultipliedSendCounts(proc_data),
		ProcData_getMultipliedDisplacements(proc_data),
		MPI_CHAR,
		Field_getDataPointer(ProcData_getPartField(proc_data, 0), M),
		ProcData_getProcessMultipliedSendCounts(proc_data),
		MPI_CHAR,
		ROOT_PROC, 
		MPI_COMM_WORLD
	);

	MPI_Request req_prev_recv;
	MPI_Request req_next_recv;
	MPI_Request req_prev_send;
	MPI_Request req_next_send;
	MPI_Request req_cycle_found;
	MPI_Status  req_status;

	int step_number = 1;
	char cycle_found = 0;

	char part_cycle_found_arr[MAX_STEPS];
	char cycle_found_arr[MAX_STEPS];

	double start_time = MPI_Wtime();

	while (step_number < MAX_STEPS){

		if (PRINT_FIELD){
			MPI_Barrier(MPI_COMM_WORLD);
			for (int i = 0; i < ProcData_getNumberOfProcs(proc_data); ++i){
				if (ProcData_getRank(proc_data) == i){
					Field_print(ProcData_getPartField(proc_data, step_number - 1));
				}
				MPI_Barrier(MPI_COMM_WORLD);
			}
			MPI_Barrier(MPI_COMM_WORLD);
			system("sleep 0.15");
			system("clear");
		}
		Field* cur_part_field = ProcData_getPartField(proc_data, step_number);
		Field* prev_part_field = ProcData_getPartField(proc_data, step_number - 1);
		int n_procs = ProcData_getNumberOfProcs(proc_data);
		int rank = ProcData_getRank(proc_data);

		// Send first row to previous proc
		MPI_Isend(
			Field_getFirstRow(prev_part_field), M, MPI_CHAR,
			Utility_getPrevIndex(rank, 0, n_procs - 1), 0,
			MPI_COMM_WORLD, &req_prev_send
		);
		
		// Send last row to next proc
		MPI_Isend(
			Field_getLastRow(prev_part_field), M, MPI_CHAR,
			Utility_getNextIndex(rank, 0, n_procs - 1), 1,
			MPI_COMM_WORLD, &req_next_send
		);
		
		// Receive pre-first row from previous proc
		MPI_Irecv(
			Field_getTFirstRow(prev_part_field), M, MPI_CHAR,
			Utility_getPrevIndex(rank, 0, n_procs - 1), 1,
			MPI_COMM_WORLD, &req_prev_recv 
		);
		
		// Receive after-last row from next proc
		MPI_Irecv(
			Field_getTLastRow(prev_part_field), M, MPI_CHAR,
			Utility_getNextIndex(rank, 0, n_procs - 1), 0,
			MPI_COMM_WORLD, &req_next_recv 
		);
		
		findCycle(ProcData_getPartFields(proc_data), part_cycle_found_arr, step_number - 1);
		MPI_Iallreduce(part_cycle_found_arr, cycle_found_arr, step_number - 1, MPI_CHAR, MPI_LAND, MPI_COMM_WORLD, &req_cycle_found);

		Field_updateMainPart(prev_part_field, cur_part_field);

		MPI_Wait(&req_prev_send, &req_status);
		MPI_Wait(&req_prev_recv, &req_status);
		Field_updateFirstRow(prev_part_field, cur_part_field);

		MPI_Wait(&req_next_send, &req_status);
		MPI_Wait(&req_next_recv, &req_status);
		Field_updateLastRow(prev_part_field, cur_part_field);

		MPI_Wait(&req_cycle_found, &req_status);
		for (int i = 0; i < step_number - 1; ++i){
			if (cycle_found_arr[i] == 1){
				cycle_found = 1;
				break;
			}
		}
		if (cycle_found){
			break;
		}
		step_number++;
	}

	double end_time = MPI_Wtime();

	if (ProcData_getRank(proc_data) == ROOT_PROC){
		if (cycle_found){
			printf("FOUND CYCLE AFTER %d STEPS.\n", step_number);
		}
		else {
			printf("RAN OUT OF STEPS.\n");
		}
		printf("Time of execution: %.5lfs.\n", end_time - start_time);
	}

	MPI_Finalize();

	return 0;
}