#include <cstdio>
#include <mpi.h>

#define TASK_REQUEST_TAG 0

pthread_t thrs[1];


/*

thread 1 (executor):


thread 2 (manager):


*/

double executor_function(void* args_s){
    ThreadArgs* args = (ThreadArgs*)args_s;

    // summary result of all executed by this process tasks
    double global_res = 0;

    for (int i = args->start; i < args->end; ++i){
        pthread_mutex_lock(&execution_mutex);
        Task* cur_task = args->task_arr[i];
        if (Task_getStatus(cur_task) == NOT_EXECUTED){
            Task_setStatus(cur_task, EXECUTED);
        }
        pthread_mutex_unlock(&execution_mutex);
        global_res += Task_execute(cur_task);
    }

    while (true){
        for (int proc_num = 0; i < args->n_procs; ++i){
            if (args->rank != proc_num){
                MPI_Send(&args->rank, 1, MPI_INT, proc_num, TASK_REQUEST_TAG, MPI_COMM_WORLD); // send "give me task!"
                MPI_Recv(&new_task_num, 1, MPI_INT, proc_num, TASK_ANSWER_TAG, MPI_COMM_WORLD, ); // receive task number to execute
                if (){

                }
            }
        }
    }

}

void manager_function(void* args_s){
    ThreadArgs* args = (ThreadArgs*)args_s;

    while (true){
        // wait for signal from executor
        MPI_Recv() // wait for signal from proc (rank == requester)
        pthread_mutex_lock(&execution_mutex);
        for (int i = args->start; i < args->end; ++i){
            Task* cur_task = args->task_arr[i];
            if (Task_getStatus(cur_task) == NOT_EXECUTED){
                Task_setStatus(cur_task, EXECUTED);
                MPI_Send() // send task num to requester
                break;
            }
        }
        pthread_mutex_unlock(&execution_mutex);
    }
}

int main(int argc, char* argv[]){
    int rank;
    int n_procs;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);

    int total_tasks = n_procs * TASKS_PER_PROC;

    Task** task_arr = Task_createTaskArray(total_tasks);

    int iter_count = 0;
    while (iter_count < TOTAL_ITERATIONS){
        // reset (set values by zeros) array and then fill it
        Task_resetTaskArray(task_arr, total_tasks);
        for (int i = 0; i < total_tasks; ++i){
            Task_setNewRepeat(task_arr[i], rank, iter_count, n_procs);
        }

        // init attributes of a thread which will be created
        pthread_attr_t attrs;
        pthread_attr_init(&attrs);
        pthread_attr_setdetachstate(&attrs, PTHREAD_CREATE_JOINABLE);

        ThreadArgs executor_args = ThreadArgs_createThreadArgs(
            task_arr,
            rank,
            n_procs,
        );

        // create task manager
        pthread_create(&thrs[i], &attrs, &manager_function, &thread_args);

        // while another thread is managing tasks, this thread will be executing tasks
        executor_function(thread_args);

        // free attributes since they won't be used anymore
        pthread_attr_destroy(&attrs);
        // join manager thread to executor thread
        pthread_join(thrs[0], NULL);

        ++iter_count;
    }

    return 0;
}