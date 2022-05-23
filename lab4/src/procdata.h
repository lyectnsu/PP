#pragma once

#include <cstring>
#include <cstdlib>
#include <mpi.h> 

#include "settings.h"
#include "utility.h"
#include "field.h"

typedef struct ProcData_s {
	int n_procs;
	int rank;

	Field** part_fields; // saved fields

	int* sc; // send_counts
	int* dp; // displacements
	int* scM; // multiplied by M send_counts
	int* dpM; // multiplied by M displacements
	int proc_sc; // sc[rank]
	int proc_dp; // dp[rank]
	int proc_scM; // scM[rank]
	int proc_dpM; // dpM[rank]
} ProcData;

ProcData* ProcData_create(){
	ProcData* proc_data = (ProcData*)malloc(sizeof(ProcData));
	
	proc_data->n_procs = -1;
	proc_data->rank = -1;

	proc_data->part_fields = NULL;

	proc_data->sc = NULL;
	proc_data->dp = NULL;
	proc_data->scM = NULL;
	proc_data->dpM = NULL;
	proc_data->proc_sc = 0;
	proc_data->proc_dp = 0;
	proc_data->proc_scM = 0;
	proc_data->proc_dpM = 0;

	return proc_data;	
}

void ProcData_init(ProcData* proc_data){
	MPI_Comm_size(MPI_COMM_WORLD, &proc_data->n_procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_data->rank);

	proc_data->sc = Utility_allocateIntArray(proc_data->n_procs);
	proc_data->dp = Utility_allocateIntArray(proc_data->n_procs);
	proc_data->scM = Utility_allocateIntArray(proc_data->n_procs);
	proc_data->dpM = Utility_allocateIntArray(proc_data->n_procs);

	// compute scatter arrays 
	int remaining = N % proc_data->n_procs;
    int sum = 0;

    for (int i = 0; i < proc_data->n_procs; i++) {
        proc_data->sc[i] = N / proc_data->n_procs;

        if (remaining > 0) {
            proc_data->sc[i]++;
            remaining--;
        }

        proc_data->dp[i] = sum;
        sum += proc_data->sc[i];

        proc_data->scM[i] = proc_data->sc[i] * M;
        proc_data->dpM[i] = proc_data->dp[i] * M;
    }

    proc_data->proc_sc = proc_data->sc[proc_data->rank];
	proc_data->proc_dp = proc_data->dp[proc_data->rank];
	proc_data->proc_scM = proc_data->scM[proc_data->rank];
	proc_data->proc_dpM = proc_data->dpM[proc_data->rank];

    // create fields
	proc_data->part_fields = (Field**)malloc(MAX_STEPS * sizeof(Field*));
	for (int i = 0; i < MAX_STEPS; ++i){
		proc_data->part_fields[i] = Field_create(proc_data->proc_sc + 2, M);
	}
}

int ProcData_getNumberOfProcs(ProcData* proc_data){
	return proc_data->n_procs;
}

int ProcData_getRank(ProcData* proc_data){
	return proc_data->rank;
}

Field** ProcData_getPartFields(ProcData* proc_data){
	return proc_data->part_fields;
}

Field* ProcData_getPartField(ProcData* proc_data, int idx){
	return proc_data->part_fields[idx];
}

int* ProcData_getMultipliedSendCounts(ProcData* proc_data){
	return proc_data->scM;
}

int* ProcData_getMultipliedDisplacements(ProcData* proc_data){
	return proc_data->dpM;
}

int* ProcData_getSendCounts(ProcData* proc_data){
	return proc_data->sc;
}

int* ProcData_getDisplacements(ProcData* proc_data){
	return proc_data->dp;
}

int ProcData_getProcessMultipliedSendCounts(ProcData* proc_data){
	return proc_data->proc_scM;
}

int ProcData_getProcessMultipliedDisplacements(ProcData* proc_data){
	return proc_data->proc_dpM;
}

int ProcData_getProcessSendCounts(ProcData* proc_data){
	return proc_data->proc_sc;
}

int ProcData_getProcessDisplacements(ProcData* proc_data){
	return proc_data->proc_dp;
}