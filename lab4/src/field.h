#pragma once

#include <cstdlib>

#include "utility.h"
#include "settings.h"

typedef struct Field_s {
	int h;
	int w;
	char* data;
} Field;

Field* Field_create(int h, int w){
	Field* field = (Field*)malloc(sizeof(Field));
	field->h = h;
	field->w = w;
	field->data = Utility_allocateCharArray(h * w);
	return field;
}

// get poiter to data with offset.
// Equal to &(field->data[offset])
char* Field_getDataPointer(Field* field, int offset){
	if (field == NULL){
		return NULL;
	}
	return field->data + offset;
}

// Get first significant row. This row belong to process which field
// is passed to the function
char* Field_getFirstRow(Field* field){
	return Field_getDataPointer(field, field->w);
}

// Get last significant row
char* Field_getLastRow(Field* field){
	return Field_getDataPointer(field, field->w * (field->h - 2));
}

// Get absolute first row. This row will be filled by values received 
// from another process. This row does not belong to process which field
// is passed to the function
char* Field_getTFirstRow(Field* field){
	return Field_getDataPointer(field, 0);
}

// Get absolute last row
char* Field_getTLastRow(Field* field){
	return Field_getDataPointer(field, field->w * (field->h - 1));
}

char Field_getCellState(Field* field, int i, int j){
	return field->data[i * field->w + j];
}

void Field_setCellState(Field* field, int i, int j, char cell_state){
	field->data[i * field->w + j] = cell_state;
}

char Field_getNextCellState(Field* field, int i, int j){

	char prev_cell_state = Field_getCellState(field, i, j);

	int is[8], js[8];

	/*
	 *  0 | 1 | 2
	 * ---+---+---
	 *  3 |///| 4 
	 * ---+---+---
	 *  5 | 6 | 7
	*/

	is[0] = Utility_getPrevIndex(i, 0, field->h - 1);
	is[1] = Utility_getPrevIndex(i, 0, field->h - 1);
	is[2] = Utility_getPrevIndex(i, 0, field->h - 1);
	is[3] = i;
	is[4] = i;
	is[5] = Utility_getNextIndex(i, 0, field->h - 1);
	is[6] = Utility_getNextIndex(i, 0, field->h - 1);
	is[7] = Utility_getNextIndex(i, 0, field->h - 1);

	js[0] = Utility_getPrevIndex(j, 0, field->w - 1);
	js[3] = Utility_getPrevIndex(j, 0, field->w - 1);
	js[5] = Utility_getPrevIndex(j, 0, field->w - 1);
	js[1] = j;
	js[6] = j;
	js[2] = Utility_getNextIndex(j, 0, field->w - 1);
	js[4] = Utility_getNextIndex(j, 0, field->w - 1);
	js[7] = Utility_getNextIndex(j, 0, field->w - 1);

	int alive_around = 0;
	for (int k = 0; k < 8; ++k){
		int neighbor_state = Field_getCellState(field, is[k], js[k]);
		if (neighbor_state == 1) {
			++alive_around;
		}
	}

	if (prev_cell_state == 0){
		if (alive_around == 3){
			return 1;
		}
	}
	else {
		if (alive_around == 2 || alive_around == 3){
			return 1;
		}
	}
	return 0;
}

// "s" stands for "start" and "e" stands for "end"
void Field_update(Field* prev_field, Field* cur_field, int srow, int scol, int erow, int ecol){
	for (int i = srow; i < erow; ++i){
		for (int j = scol; j < ecol; ++j){
			char cell_state = Field_getNextCellState(prev_field, i, j);
			Field_setCellState(cur_field, i, j, cell_state);
		}
	}
}

// In all functions below consider heights and widths of the fields are the same
void Field_updateMainPart(Field* prev_field, Field* cur_field){
	Field_update(prev_field, cur_field, 2, 0, cur_field->h - 2, cur_field->w);
}

void Field_updateFirstRow(Field* prev_field, Field* cur_field){
	Field_update(prev_field, cur_field, 1, 0, 2, cur_field->w);
	
}

void Field_updateLastRow(Field* prev_field, Field* cur_field){
	Field_update(prev_field, cur_field, cur_field->h - 2, 0, cur_field->h - 1, cur_field->w);
}

char Field_equal(Field* f1, Field* f2){
	for (int i = 1; i < f1->h - 1; ++i){
		for (int j = 0; j < f1->w; ++j){
			if (Field_getCellState(f1, i, j) != Field_getCellState(f2, i, j)){
				return 0;
			}
		}
	}
	return 1;
}

void Field_fillRandom(Field* field){
	for (int i = 1; i < field->h - 1; ++i){
		for (int j = 0; j < field->w; ++j){
			char cell_state = rand() % 2;
			Field_setCellState(field, i, j, cell_state);
		}
	}
}

void Field_placeGlider(Field* field, int i, int j){
	Field_setCellState(field,     i,     j, 0);
	Field_setCellState(field,     i, j + 1, 1);
	Field_setCellState(field,     i, j + 2, 0);
	Field_setCellState(field, i + 1,     j, 0);
	Field_setCellState(field, i + 1, j + 1, 0);
	Field_setCellState(field, i + 1, j + 2, 1);
	Field_setCellState(field, i + 2,     j, 1);
	Field_setCellState(field, i + 2, j + 1, 1);
	Field_setCellState(field, i + 2, j + 2, 1);
}

void Field_print(Field* field){
	for (int i = 1; i < field->h - 1; ++i){
		for (int j = 0; j < field->w; ++j){
			char cell_state = Field_getCellState(field, i, j);
			if (cell_state == 1) printf("@");
			if (cell_state == 0) printf(".");
		}
		printf("\n");
	}
}