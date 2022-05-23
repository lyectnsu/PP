#pragma once

#include <cstdlib>
#include <cstring>

int* Utility_allocateIntArray(int size){
	int* array = (int*)malloc(size * sizeof(int));
	memset(array, 0, size);
	return array;
}

char* Utility_allocateCharArray(int size){
	char* array = (char*)malloc(size * sizeof(char));
	memset(array, 0, size);
	return array;
}

int Utility_getPrevIndex(int icur, int imin, int imax){
	--icur;
	if (icur == imin - 1){
		icur = imax;
	}
	return icur;
}

int Utility_getNextIndex(int icur, int imin, int imax){
	++icur;
	if (icur == imax + 1){
		icur = imin;
	}
	return icur;	
}

void Utility_printCharArray(char* arr, int size){
	for (int i = 0; i < size; ++i){
		printf("%d ", arr[i]);
	}
	printf("\n");
}