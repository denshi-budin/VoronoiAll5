#include <stdio.h>
#include <math.h>
#include "atom.h"

#define RAD(v) (v)*(3.14159265359 / 180.0)

position_t arotate3d(position_t atoms, position_t cent, dpoint_t data)
{
	position_t dummy = {0., 0., 0.};
	position_t old;

	old.x = atoms.x - cent.x;
	old.y = atoms.y - cent.y;
	old.z = atoms.z - cent.z;

	dummy.x = (old.x * data.A[0][0]) + (old.y * data.A[0][1]) + (old.z * data.A[0][2]) + cent.x;//111012
	dummy.y = (old.x * data.A[1][0]) + (old.y * data.A[1][1]) + (old.z * data.A[1][2]) + cent.y;//010002
	dummy.z = (old.x * data.A[2][0]) + (old.y * data.A[2][1]) + (old.z * data.A[2][2]) + cent.z;//212022

	atoms.x = dummy.x;
	atoms.y = dummy.y;
	atoms.z = dummy.z;

	return atoms;
}

position_t arotate2d(position_t atoms, position_t cent, dpoint_t data)
{

	position_t dummy = {0., 0., 0.};

	double oldX = atoms.x - cent.x;
	double oldY = atoms.y - cent.y;

		/* 回転させて(0,0,0) -> (latticeLengthX/2,latticeLengthY/2,latticeLengthZ/2) だけ平行移動 */
	dummy.x = oldX * data.A[0][0] + oldY * data.A[0][1] + cent.x;
	dummy.y = oldX * data.A[1][0] + oldY * data.A[1][1] + cent.y;

	atoms.x = dummy.x;
	atoms.y = dummy.y;

	return atoms;
}

position_t rotate(position_t atoms, position_t cent, dpoint_t data)
{
	if(pointVor == 1) return atoms;
	if(latticeSizeZ == 1) return arotate2d(atoms, cent, data);
	else return arotate3d(atoms, cent, data);
}

void eulervalue(dpoint_t datapoint[])
{
	int editval = 0;
	int pointN = 0;
	double ppsi = 0.0, ptheta = 0.0, pphi = 0.0;
	double space;
	if (remake == 0){
		if (latticeSizeZ == 1) space = 90.0 / (double)pointVor;
		else space = 180.0 / (double)pointVor;
		for (int i = 0; i < mN; ++i)
		{
			for (int j = 0; j < 3; ++j) for (int k = 0; k < 3; ++k) datapoint[i].A[j][k] = 0.0;

			if(datapoint[i].oriPoint != -1) continue;

			if(latticeSizeZ == 1){
				do {
					datapoint[i].theta = GetRandom(5, 90);
				} while (fabs(ptheta - datapoint[i].theta) <= space);
			}
			else{
				do {
					datapoint[i].psi = GetRandom(5, 180);
				} while (fabs(ppsi - datapoint[i].psi) <= space);
				do {
					datapoint[i].theta = GetRandom(5, 180);
				} while (fabs(ptheta - datapoint[i].theta) <= space);
				do {
					datapoint[i].phi = GetRandom(5, 180);
				} while  (fabs(pphi - datapoint[i].phi) <= space);
			}
			printf("Euler val point %d (psi:theta:phi): %g : %g : %g\n", i+1, datapoint[i].psi, datapoint[i].theta, datapoint[i].phi);
			ppsi = datapoint[i].psi; 
			ptheta = datapoint[i].theta;
			pphi = datapoint[i].phi;
		}

		if (pointVor != 1){
			printf("Edit euler val? 1:yes 0:no : ");
			scanf("%d", &editval);
			while(editval >= 1){
				do {
					printf("Input point that need to change: ");
					scanf("%d", &pointN);
				} while ((pointN < 1) || (pointN > pointVor)); 
				if(latticeSizeZ == 1){
					printf("Input new theta val for point %d:-\n", pointN);
					scanf("%lf", &datapoint[pointN-1].theta);
				}
				else{
					printf("Input new psi val for point %d:-\n", pointN);
					scanf("%lf", &datapoint[pointN-1].psi);
					printf("Input new theta val for point %d:-\n", pointN);
					scanf("%lf", &datapoint[pointN-1].theta);
					printf("Input new phi val for point %d:-\n", pointN);
					scanf("%lf", &datapoint[pointN-1].phi);
				}
				pointN = 0;
				printf("Edit other point? 1:yes 0:no : ");
				scanf("%d", &editval);
			}
		}
	}
	else for (int i = 0; i < mN; ++i) for (int j = 0; j < 3; ++j) for (int k = 0; k < 3; ++k) datapoint[i].A[j][k] = 0.0; 


	for (int i = 0; i < mN; ++i)
	{
		if (latticeSizeZ == 1){
			datapoint[i].A[0][0] = cos(RAD(datapoint[i].theta));//オイラー角行列公式
			datapoint[i].A[0][1] =-sin(RAD(datapoint[i].theta));
			datapoint[i].A[1][0] = sin(RAD(datapoint[i].theta));
			datapoint[i].A[1][1] = cos(RAD(datapoint[i].theta));
		}
		else {
			datapoint[i].A[0][0] = (cos(RAD(datapoint[i].psi)) * cos(RAD(datapoint[i].phi))) - (cos(RAD(datapoint[i].theta)) * sin(RAD(datapoint[i].phi)) * sin(RAD(datapoint[i].psi)));//オイラー角行列公式
			datapoint[i].A[0][1] = (cos(RAD(datapoint[i].psi)) * sin(RAD(datapoint[i].phi))) + (cos(RAD(datapoint[i].theta)) * cos(RAD(datapoint[i].phi)) * sin(RAD(datapoint[i].psi)));
			datapoint[i].A[0][2] = sin(RAD(datapoint[i].psi)) * sin(RAD(datapoint[i].theta));
			datapoint[i].A[1][0] = (-sin(RAD(datapoint[i].psi)) * cos(RAD(datapoint[i].phi))) - (cos(RAD(datapoint[i].theta)) * sin(RAD(datapoint[i].phi)) * cos(RAD(datapoint[i].psi)));
			datapoint[i].A[1][1] = (-sin(RAD(datapoint[i].psi)) * sin(RAD(datapoint[i].phi))) + (cos(RAD(datapoint[i].theta)) * cos(RAD(datapoint[i].phi)) * cos(RAD(datapoint[i].psi)));
			datapoint[i].A[1][2] = cos(RAD(datapoint[i].psi)) * sin(RAD(datapoint[i].theta));
			datapoint[i].A[2][0] = sin(RAD(datapoint[i].theta)) * sin(RAD(datapoint[i].phi));
			datapoint[i].A[2][1] = -sin(RAD(datapoint[i].theta)) * cos(RAD(datapoint[i].phi));
			datapoint[i].A[2][2] = cos(RAD(datapoint[i].theta));
		}

		if(datapoint[i].oriPoint == -1) continue;
		datapoint[i].psi = datapoint[datapoint[i].oriPoint].psi;
		datapoint[i].theta = datapoint[datapoint[i].oriPoint].theta;
		datapoint[i].phi = datapoint[datapoint[i].oriPoint].phi;
		for (int j = 0; j < 3; ++j) for (int k = 0; k < 3; ++k) 
			datapoint[i].A[j][k] = datapoint[datapoint[i].oriPoint].A[j][k];
	}
}
