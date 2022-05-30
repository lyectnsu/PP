#include <cstdio>
#include <mpi.h>
#include <pthread.h>

#include "task.h"
#include "settings.h"

#define TASK_REQUEST_TAG 1
#define TASK_ANSWER_TAG 2

typedef struct ThreadArgs_s {
	int start;
	int end;
	int rank;
	int n_procs;
	Task** task_arr;
	double* tasks_result;
	int* tasks_done;
	int* tasks_sent;
	double exec_time;
} ThreadArgs;

pthread_t thrs[1];

pthread_mutex_t execution_mutex;  

ThreadArgs* ThreadArgs_createThreadArgs(
	int start, int end, int rank, int n_procs, Task** task_arr, double* tasks_result, int* tasks_done, int* tasks_sent
){
	ThreadArgs* ta = (ThreadArgs*)malloc(sizeof(ThreadArgs));
	ta->start = start;
	ta->end = end;
	ta->rank = rank;
	ta->n_procs = n_procs;
	ta->task_arr = task_arr;
	ta->tasks_result = tasks_result;
	ta->tasks_done = tasks_done;
	ta->tasks_sent = tasks_sent;
	ta->exec_time = 0;
	return ta;
}

double calculate_time(struct timespec start, struct timespec end){
	return end.tv_sec - start.tv_sec + 0.000000001 * (double)(end.tv_nsec - start.tv_nsec);
}

void* executor_function(void* args_s){
	ThreadArgs* args = (ThreadArgs*)args_s;

	struct timespec start, end;
	

	bool task_taken;
	for (int i = args->start; i < args->end; ++i){
		task_taken = false;
		pthread_mutex_lock(&execution_mutex);
		Task* cur_task = args->task_arr[i];
		if (Task_getStatus(cur_task) == NOT_EXECUTED){
			Task_setStatus(cur_task, EXECUTED);
			task_taken = true;
		}
		pthread_mutex_unlock(&execution_mutex);
		if (task_taken){
			clock_gettime(CLOCK_MONOTONIC_RAW, &start);
			*args->tasks_result += Task_execute(cur_task);
			clock_gettime(CLOCK_MONOTONIC_RAW, &end);
			*args->exec_time += calculate_time(start, end);
			*args->tasks_done += 1;
		}
	}

	bool may_have_task[args->n_procs];
	int empty_procs = 0;
	
	for (int i = 0; i < args->n_procs; ++i){
		may_have_task[i] = true;
	}
	may_have_task[args->rank] = false;

	int new_task_num;
	while (true){
		if (empty_procs == args->n_procs - 1){
			return NULL;
		}
		for (int proc_num = 0; proc_num < args->n_procs; ++proc_num){
			if (may_have_task[proc_num]){
				MPI_Send(&args->rank, 1, MPI_INT, proc_num, TASK_REQUEST_TAG, MPI_COMM_WORLD); // send "give me task!"
				MPI_Recv(&new_task_num, 1, MPI_INT, proc_num, TASK_ANSWER_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // receive task number to execute
				if (new_task_num == -1){
					may_have_task[proc_num] = false;
					empty_procs++;
				}
				else {
					Task* cur_task = args->task_arr[new_task_num];
					clock_gettime(CLOCK_MONOTONIC_RAW, &start);
					*args->tasks_result += Task_execute(cur_task);
					clock_gettime(CLOCK_MONOTONIC_RAW, &end);
					*args->exec_time += calculate_time(start, end);
					*args->tasks_done += 1;
				}
			}
		}
	}
}

void* manager_function(void* args_s){
	ThreadArgs* args = (ThreadArgs*)args_s;

	int requester;
	int rejected = 0;
	bool proc_is_empty = false;
	while (true){
		if (rejected == args->n_procs - 1){
			pthread_exit(0);
		}
		
		// wait for signal from executor from another proc
		MPI_Recv(&requester, 1, MPI_INT, MPI_ANY_SOURCE, TASK_REQUEST_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		if (!proc_is_empty){
			bool task_sent = false;
			for (int i = args->start; i < args->end; ++i){
				pthread_mutex_lock(&execution_mutex);
				Task* cur_task = args->task_arr[i];
				if (Task_getStatus(cur_task) == NOT_EXECUTED){
					Task_setStatus(cur_task, EXECUTED);
					task_sent = true;
				}
				pthread_mutex_unlock(&execution_mutex);
				if (task_sent){
					*args->tasks_sent += 1;
					MPI_Send(&i, 1, MPI_INT, requester, TASK_ANSWER_TAG, MPI_COMM_WORLD); // send task num to requester
					break;
				}
			}
			
			if (!task_sent){
				proc_is_empty = true;
			}
		}
		if (proc_is_empty){
			int no_task = -1;
			MPI_Send(&no_task, 1, MPI_INT, requester, TASK_ANSWER_TAG, MPI_COMM_WORLD);
			rejected++;
		}
	}
}

int main(int argc, char* argv[]){
	int rank;
	int n_procs;

	int provided_level;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided_level);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &n_procs);

	int total_tasks = n_procs * TASKS_PER_PROC;

	Task** task_arr = Task_createTaskArray(total_tasks);
	pthread_mutex_init(&execution_mutex, NULL);

	int iter_count = 0;

	double total_time_start = MPI_Wtime();
	while (iter_count < TOTAL_ITERATIONS){
		// reset (set values by zeros) array and then fill it
		Task_resetTaskArray(task_arr, total_tasks);
		for (int i = 0; i < n_procs; ++i){
			for (int j = 0; j < TASKS_PER_PROC; ++j){
				Task_setNewRepeat(task_arr[i * TASKS_PER_PROC + j], i, iter_count, n_procs);
			}
		}

		// init attributes of a thread which will be created
		pthread_attr_t attrs;
		pthread_attr_init(&attrs);
		pthread_attr_setdetachstate(&attrs, PTHREAD_CREATE_JOINABLE);

		double tasks_result = 0;
		int tasks_done = 0;
		int tasks_sent = 0;

		ThreadArgs* thread_args = ThreadArgs_createThreadArgs(
			TASKS_PER_PROC * rank,
			TASKS_PER_PROC * (rank + 1),
			rank,
			n_procs,
			task_arr,
			&tasks_result,
			&tasks_done,
			&tasks_sent
		);

		// create task manager
		pthread_create(&thrs[0], &attrs, manager_function, thread_args);
		// while another thread is managing tasks, this thread will be executing tasks
		executor_function(thread_args);

		// join manager thread to executor thread
		pthread_join(thrs[0], NULL);
		// free attributes since they won't be used anymore
		pthread_attr_destroy(&attrs);
		pthread_mutex_destroy(&execution_mutex);

		double min_time, max_time;
		MPI_Allreduce(&thread_args->exec_time, &min_time, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		MPI_Allreduce(&thread_args->exec_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

		if (rank == 0){
			printf("Iteration #%d:\n", iter_count + 1);
			printf("\tDisbalance: %.5lf %%\n", (max_time - min_time) / max_time * 100);
		}

		for (int i = 0; i < n_procs; ++i){
			MPI_Barrier(MPI_COMM_WORLD);
			if (rank == i){
				printf("\tPROCESS WITH RANK #%d\n", rank);
				printf("\t\tExecution time: %lf s.\n", thread_args->exec_time);
				printf("\t\tTasks result:   %lf\n", tasks_result);
				printf("\t\tTasks done:     %d\n", tasks_done);
				printf("\t\tTasks sent:     %d\n", tasks_sent);
			}
			fflush(stdout);
			MPI_Barrier(MPI_COMM_WORLD);
		}

		free(thread_args);

		++iter_count;
	}
	double total_time_end = MPI_Wtime();

	if (rank == 0){
		printf("Total time: %.5lf s.\n", total_time_end - total_time_start);
	}

	Task_freeTaskArray(task_arr, total_tasks);

	MPI_Finalize();
	return 0;
}