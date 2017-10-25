#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "atom.h"

/*小さいから大きい順を並び替え関数*/
/*DdataとswapByが同じ構造体のなかにあるとだめです*/
void swap_single(double *data, int val_i, int val_j) {
	double swap;

	swap = data[val_i];
	data[val_i] = data[val_j];
	data[val_j] = swap;
}

void swap_3d(position_t *data, int val_i, int val_j) {
	position_t swap;

	swap = data[val_i];
	data[val_i] = data[val_j];
	data[val_j] = swap;
}

int partition(position_t *data, double *swpby, int left, int right, double pivot) {

	int leftpoint = left - 1;
	int rightpoint = right;

	while (true)
	{
		while (swpby[++leftpoint] < pivot);
		while (rightpoint > 0 && swpby[--rightpoint]>pivot);

		if (leftpoint >= rightpoint) break;
		else {
			swap_single(swpby, leftpoint, rightpoint);
			swap_3d(data, leftpoint, rightpoint);
		}
	}

	swap_single(swpby, leftpoint, right);
	swap_3d(data, leftpoint, right);
	return leftpoint;
}

void quicksort(position_t *data, double *swpby, int left, int right) {
	if (right - left <= 0) {
		return;
	}
	else {
		double pivot = swpby[right];
		int partitionPoint = partition(data, swpby, left, right, pivot);
		quicksort(data, swpby, left, partitionPoint - 1);
		quicksort(data, swpby, partitionPoint + 1, right);
	}
}

void swaping_small_big(position_t Ddata[], double swapBy[], int Ndata, const char *fname, double adjust){
	FILE *ff;
	char ffname[30];
	strcpy(ffname, fname);
	strcat(ffname, ".csv");
	ff = fopen(ffname, "w");
	position_t swap;
	double dataSwap;

	if (adjust == 0.0) adjust = 1.0;

	quicksort(Ddata, swapBy, 0, Ndata - 1);

	for (int i = 0; i < Ndata; i++) {
		fprintf(ff, "%g,%g,%g,%g\n", Ddata[i].x/adjust, Ddata[i].y/adjust, Ddata[i].z/adjust, swapBy[i]/adjust);
	}
	fclose(ff);
}