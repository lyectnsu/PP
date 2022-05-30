#pragma once

#include <cstdlib>
#include <cmath>

#include "settings.h"

enum TaskStatus {
    NOT_EXECUTED,
    EXECUTED
};

typedef struct Task_s{
    int repeat;
    TaskStatus status;
} Task;

void Task_resetTask(Task* task){
    task->repeat = 0;
    task->status = NOT_EXECUTED;
}

void Task_resetTaskArray(Task** task_arr, int size){
    for (int i = 0; i < size; ++i){
        Task_resetTask(task_arr[i]);
    }
}

Task* Task_createTask(){
    Task* task = (Task*)malloc(sizeof(Task));
    Task_resetTask(task);
    return task;
}

Task** Task_createTaskArray(int size){
    Task** task_arr = (Task**)malloc(size * sizeof(Task*));
    for (int i = 0; i < size; ++i){
        task_arr[i] = Task_createTask();
    }
    return task_arr;
}

void Task_freeTaskArray(Task** task_arr, int size){
    for (int i = 0; i < size; ++i){
        free(task_arr[i]);
    }
    free(task_arr);
}

void Task_setNewRepeat(Task* task, int rank, int iter_count, int n_procs){
    task->repeat = TASK_QUANT_TIME * (1 + abs(rank - (iter_count % n_procs)));
}

int Task_getRepeat(Task* task){
    return task->repeat;
}

TaskStatus Task_getStatus(Task* task){
    return task->status;
}

void Task_setStatus(Task* task, TaskStatus new_status){
    task->status = new_status;
}

double Task_execute(Task* task){
    double local_res = 0;
    for (int i = 0; i < task->repeat; ++i){
        local_res += sin(i);
    }
    return local_res;
}
