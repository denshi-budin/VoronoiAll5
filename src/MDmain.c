#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "atom.h"

void fail(const char *error){
    printf("Error occurred: ");
    printf("\"%s\"\n", error);
    exit(EXIT_FAILURE);
}

int main(int argc, char *argv[])
{
    double timer = 0, timer_past = 0;
    int mins = 0, hour = 0, day = 0;
    int pointuse = 0;
    time_t t;
    /* 構造体宣言 */
    settingParameter ();
    atoms = (position_t *) calloc (RNum, sizeof(position_t));
    adis = (double *) calloc (RNum, sizeof(double));
    nearAtm = (position_t *) calloc (N*2, sizeof(position_t));
    if (remake == 0){
        mpoint = (mpoint_t *) calloc (mN+1, sizeof(mpoint_t));
        datapoint = (dpoint_t *) calloc (mN+1, sizeof(dpoint_t));
    }

    for (int i = 0; i < N*2 ; i++) nearAtm[i].x = nearAtm[i].y = nearAtm[i].z = -9999;

    StartTimer();
    if (argc >= 2){
        randseed = atoi(argv[1]);
    }
    else if (remake != 0); 
    else randseed = (unsigned) time(&t);
    if (pointVor != 1) printf("Rand seed: %d\n", randseed);
    srand(randseed);
	/* 初期化 */
    RandomBox(mpoint, datapoint);
    copyPoint(mpoint, datapoint);
    eulervalue(datapoint);
    timer = GetTimer();
    if (pointVor != 1) {
        printf("Point arrange at ");
        PrintTimer(timer - timer_past);
    }
    timer_past = timer;

	/* ボロノイ分割のメインループ */
    for (int i = 0; i < mN; i++) {
        if (!inSide(mpoint[i].pos, latorigin, latticeLength)) continue;
        pointuse++;
    }
    printf("Total voronoi point used: %d\n", pointuse);
    pointuse = 0;

	for(int v=0;v<mN;v++){
        if (!inSide(mpoint[v].pos, latorigin, latticeLength)) continue;
        
        MakeRotateEle(atoms, mpoint, datapoint, v);

        timer = GetTimer();
        if (pointVor != 1) {
            printf("Done arrange voronoi point %d at ", ++pointuse);
            PrintTimer(timer - timer_past);
        }
        timer_past = timer;
	}

    relocate(mpoint);
    //checkrange(mpoint);
	writeFile(mpoint, datapoint);
    printf("Time taken: ");
    PrintTimer(GetTimer());
	return 0;
}
