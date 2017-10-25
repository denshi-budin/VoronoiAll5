#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "atom.h"

int subcell;
bool ev, od;

void Vregister(position_t atoms[], mpoint_t mpoint[], int total, int pt)
{
	static int Nv = 0;
	mpoint[pt].voronoi = (position_t *) calloc(total + 1, sizeof(position_t));
	mpoint[pt].dis = (double *) calloc(total + 1, sizeof(double));
	mpoint[pt].voronoiNum = total;

	for (int i = 0; i < total; ++i)
	{
			mpoint[pt].voronoi[i] = atoms[i];
			mpoint[pt].dis[i] = adis[i];
	}

	if(inSide(mpoint[pt].pos,simboxori,latticeLength)){
		avgsizeV = (avgsizeV * Nv) + sizeVoronoi(mpoint[pt].voronoi,mpoint[pt].voronoiNum, mpoint[pt].dis);
		avgsizeV /= (double)Nv + 1.0;
		Nv += 1;
	}
}

bool even(int a) {
	if (a % 2 == 0) return true;
	else return false;
}

bool insideb(int x, int y, int z, int latX, int latY, int latZ) {
	if (x >= 0 && x < latX && y >= 0 && y < latY && z >= 0 && z < latZ) 
		return true;
	else return false;
}

bool value(int num, int point, int pt, int count, int *k) {
	int val = 0;
	if (even(num)) {
		val = point + num;
		(*k) = val;
		if ((count == 0) && (pt < pointVor)) od = false;
		if (!ev) return false;
	}
	else {
		val = point - num;
		(*k) = val;
		if ((count == 0) && (pt < pointVor)) ev = false;
		if (!od) return false;
	}
	return true;
}

int FCC(position_t atoms[], mpoint_t mpoint[], dpoint_t datapoint[], 
	int pt, int num, int i, int j, int k, position_t lengthCenter)
{
	int bil = num;
	position_t originPlace;
	originPlace.x =  mpoint[pt].pos.x - lengthCenter.x;
	originPlace.y =  mpoint[pt].pos.y - lengthCenter.y;
	originPlace.z =  mpoint[pt].pos.z - lengthCenter.z;

	atoms[bil].x = originPlace.x + (latticeConst.x * i);
	atoms[bil].y = originPlace.y + (latticeConst.y * j);
	atoms[bil].z = originPlace.z + (latticeConst.z * k);
	atoms[bil] = rotate(atoms[bil], mpoint[pt].pos, datapoint[pt]);
	if(ainOrOut(atoms[bil], mpoint, &adis[bil], pt)) bil += 1;

	atoms[bil].x = originPlace.x + (latticeConst.x * (i + 0.5));
	atoms[bil].y = originPlace.y + (latticeConst.y * (j + 0.5));
	atoms[bil].z = originPlace.z + (latticeConst.z * k);
	atoms[bil] = rotate(atoms[bil], mpoint[pt].pos, datapoint[pt]);
	if(ainOrOut(atoms[bil], mpoint, &adis[bil], pt)) bil += 1;

	atoms[bil].x = originPlace.x + (latticeConst.x * i);
	atoms[bil].y = originPlace.y + (latticeConst.y * (j + 0.5));
	atoms[bil].z = originPlace.z + (latticeConst.z * (k + 0.5));
	atoms[bil] = rotate(atoms[bil], mpoint[pt].pos, datapoint[pt]);
	if(ainOrOut(atoms[bil], mpoint, &adis[bil], pt)) bil += 1;

	atoms[bil].x = originPlace.x + (latticeConst.x * (i + 0.5));
	atoms[bil].y = originPlace.y + (latticeConst.y * j);
	atoms[bil].z = originPlace.z + (latticeConst.z * (k + 0.5));
	atoms[bil] = rotate(atoms[bil], mpoint[pt].pos, datapoint[pt]);
	if(ainOrOut(atoms[bil], mpoint, &adis[bil], pt)) bil += 1;

	return bil;
}

int DIAMOND(position_t atoms[], mpoint_t mpoint[], dpoint_t datapoint[], 
	int pt, int num, int i, int j, int k, position_t lengthCenter)
{
	int bil = num;
	position_t originPlace;
	originPlace.x =  mpoint[pt].pos.x - lengthCenter.x;
	originPlace.y =  mpoint[pt].pos.y - lengthCenter.y;
	originPlace.z =  mpoint[pt].pos.z - lengthCenter.z;

	atoms[bil].x = originPlace.x + (latticeConst.x * i);
	atoms[bil].y = originPlace.y + (latticeConst.y * j);
	atoms[bil].z = originPlace.z + (latticeConst.z * k);
	atoms[bil] = rotate(atoms[bil], mpoint[pt].pos, datapoint[pt]);
	if(ainOrOut(atoms[bil], mpoint, &adis[bil], pt)) 
		bil += 1;
	atoms[bil].x = originPlace.x + (latticeConst.x * (i + 0.5));
	atoms[bil].y = originPlace.y + (latticeConst.y * (j + 0.5));
	atoms[bil].z = originPlace.z + (latticeConst.z * k);
	atoms[bil] = rotate(atoms[bil], mpoint[pt].pos, datapoint[pt]);
	if(ainOrOut(atoms[bil], mpoint, &adis[bil], pt)) 
		bil += 1;
	atoms[bil].x = originPlace.x + (latticeConst.x * i);
	atoms[bil].y = originPlace.y + (latticeConst.y * (j + 0.5));
	atoms[bil].z = originPlace.z + (latticeConst.z * (k + 0.5));
	atoms[bil] = rotate(atoms[bil], mpoint[pt].pos, datapoint[pt]);
	if(ainOrOut(atoms[bil], mpoint, &adis[bil], pt)) 
		bil += 1;
	atoms[bil].x = originPlace.x + (latticeConst.x * (i + 0.5));
	atoms[bil].y = originPlace.y + (latticeConst.y * j);
	atoms[bil].z = originPlace.z + (latticeConst.z * (k + 0.5));
	atoms[bil] = rotate(atoms[bil], mpoint[pt].pos, datapoint[pt]);
	if(ainOrOut(atoms[bil], mpoint, &adis[bil], pt)) 
		bil += 1;
	/* si */
	atoms[bil].x = originPlace.x + (latticeConst.x * (i + 0.25));
	atoms[bil].y = originPlace.y + (latticeConst.y * (j + 0.25));
	atoms[bil].z = originPlace.z + (latticeConst.z * (k + 0.25));
	atoms[bil] = rotate(atoms[bil], mpoint[pt].pos, datapoint[pt]);
	if(ainOrOut(atoms[bil], mpoint, &adis[bil], pt)) 
		bil += 1;
	atoms[bil].x = originPlace.x + (latticeConst.x * (i + 0.75));
	atoms[bil].y = originPlace.y + (latticeConst.y * (j + 0.75));
	atoms[bil].z = originPlace.z + (latticeConst.z * (k + 0.25));
	atoms[bil] = rotate(atoms[bil], mpoint[pt].pos, datapoint[pt]);
	if(ainOrOut(atoms[bil], mpoint, &adis[bil], pt)) 
		bil += 1;
	atoms[bil].x = originPlace.x + (latticeConst.x * (i + 0.75));
	atoms[bil].y = originPlace.y + (latticeConst.y * (j + 0.25));
	atoms[bil].z = originPlace.z + (latticeConst.z * (k + 0.75));
	atoms[bil] = rotate(atoms[bil], mpoint[pt].pos, datapoint[pt]);
	if(ainOrOut(atoms[bil], mpoint, &adis[bil], pt)) 
		bil += 1;
	atoms[bil].x = originPlace.x + (latticeConst.x * (i + 0.25));
	atoms[bil].y = originPlace.y + (latticeConst.y * (j + 0.75));
	atoms[bil].z = originPlace.z + (latticeConst.z * (k + 0.75));
	atoms[bil] = rotate(atoms[bil], mpoint[pt].pos, datapoint[pt]);
	if(ainOrOut(atoms[bil], mpoint, &adis[bil], pt)) 
		bil += 1;

	return bil;
}

int BCC(position_t atoms[], mpoint_t mpoint[], dpoint_t datapoint[], 
	int pt, int num, int i, int j, int k, position_t lengthCenter)
{
	int bil = num;
	position_t originPlace;
	originPlace.x =  mpoint[pt].pos.x - lengthCenter.x;
	originPlace.y =  mpoint[pt].pos.y - lengthCenter.y;
	originPlace.z =  mpoint[pt].pos.z - lengthCenter.z;

	atoms[bil].x = originPlace.x + (latticeConst.x * i);
	atoms[bil].y = originPlace.y + (latticeConst.y * j);
	atoms[bil].z = originPlace.z + (latticeConst.z * k);
	atoms[bil] = rotate(atoms[bil], mpoint[pt].pos, datapoint[pt]);
	if(ainOrOut(atoms[bil], mpoint, &adis[bil], pt)) 
		bil += 1;
	atoms[bil].x = originPlace.x + (latticeConst.x * (i + 0.5));
	atoms[bil].y = originPlace.y + (latticeConst.y * (j + 0.5));
	atoms[bil].z = originPlace.z + (latticeConst.z * (k + 0.5));
	atoms[bil] = rotate(atoms[bil], mpoint[pt].pos, datapoint[pt]);
	if(ainOrOut(atoms[bil], mpoint, &adis[bil], pt)) 
		bil += 1;

	return bil;
}

int HCP(position_t atoms[], mpoint_t mpoint[], dpoint_t datapoint[], 
	int pt, int num, int i, int j, int k, position_t lengthCenter)
{
	int bil = num;
	position_t originPlace;
	originPlace.x =  mpoint[pt].pos.x - lengthCenter.x;
	originPlace.y =  mpoint[pt].pos.y - lengthCenter.y;
	originPlace.z =  mpoint[pt].pos.z - lengthCenter.z;

	atoms[bil].x = originPlace.x + (latticeConst.x * i);
	atoms[bil].y = originPlace.y + (latticeConst.y * j);
	atoms[bil].z = originPlace.z + (latticeConst.z * k);
	atoms[bil] = rotate(atoms[bil], mpoint[pt].pos, datapoint[pt]);
	if(ainOrOut(atoms[bil], mpoint, &adis[bil], pt)) bil += 1;

	atoms[bil].x = originPlace.x + (latticeConst.x * (i + 0.5));
	atoms[bil].y = originPlace.y + (latticeConst.y * (j + 0.5));
	atoms[bil].z = originPlace.z + (latticeConst.z * k);
	atoms[bil] = rotate(atoms[bil], mpoint[pt].pos, datapoint[pt]);
	if(ainOrOut(atoms[bil], mpoint, &adis[bil], pt)) bil += 1;

	atoms[bil].x = originPlace.x + (latticeConst.x * i);
	atoms[bil].y = originPlace.y + (latticeConst.y * (j + 4. / 6.));
	atoms[bil].z = originPlace.z + (latticeConst.z * (k + 0.5));
	atoms[bil] = rotate(atoms[bil], mpoint[pt].pos, datapoint[pt]);
	if(ainOrOut(atoms[bil], mpoint, &adis[bil], pt)) bil += 1;

	atoms[bil].x = originPlace.x + (latticeConst.x * (i + 0.5));
	atoms[bil].y = originPlace.y + (latticeConst.y * (j + 1. / 6.));
	atoms[bil].z = originPlace.z + (latticeConst.z * (k + 0.5));
	atoms[bil] = rotate(atoms[bil], mpoint[pt].pos, datapoint[pt]);
	if(ainOrOut(atoms[bil], mpoint, &adis[bil], pt)) bil += 1;

	return bil;
}

int set_atom(char *setname, position_t atoms[], mpoint_t mpoint[], dpoint_t datapoint[], 
	int pt, int num, int i, int j, int k, position_t lengthCenter){

	if (strcmp("fcc", setname) == 0) return FCC(atoms, mpoint, datapoint, pt, num, i, j, k, lengthCenter);
	else if (strcmp("bcc", setname) == 0) return BCC(atoms, mpoint, datapoint, pt, num, i, j, k, lengthCenter);
	else if (strcmp("diamond", setname) == 0) return DIAMOND(atoms, mpoint, datapoint, pt, num, i, j, k, lengthCenter);
	else if(strcmp("hcp", setname) == 0) return HCP(atoms, mpoint, datapoint, pt, num, i, j, k, lengthCenter);
	else {
		char strgerr[100];
		sprintf(strgerr, "Unsupport structure: %s", setname);
		fail(strgerr);
	}
	return 0;
}

void MakeRotateEle(position_t atoms[], mpoint_t mpoint[], dpoint_t datapoint[], int pt)
{
	int num = 0, cnt = 0, cntP = 0, cntl = 0, numpast = 0;
	int stk, xfi, yfi, zfi;
	int i, j, k;
	//int breakval = 0;
	subcell = 1;
	position_t halflength;
	halflength.x = halflength.y = halflength.z = floor((SizeRotateBox+PLUS)/2) * latConstBox;
	xfi = yfi = zfi = (SizeRotateBox+PLUS);
	if (latticeSizeZ == 1) {
		zfi = (SizeRotateBoxZ+PLUS);
		halflength.z = floor((SizeRotateBoxZ+PLUS)/2) * latticeConst.z;
	}

	stk = (int)round(zfi/2);
	ev = od = true;

	for(int c=0; c<zfi; c++){
		i = (int)round(xfi/2);
		j = (int)round(yfi/2);
		//printf("c : %d  ", c);
		if(!value(c, stk, pt, subcell, &k)) {stk = k; continue;}
		//printf("od : %d ; ev %d ", od, ev );
		if ((!od) && (!ev)) break;
		stk = k;
		cnt = 0;
		int cntx = 0, cnty = 0;

		num = set_atom(lattice, atoms, mpoint, datapoint, pt, num, i, j, k, halflength);

		int subcnt = 0;
		while (cnt < (xfi + yfi - 1))
		{
			cntP = 0;
			if (even(cntl)) {
				if (even(cntx)) {
					cntx++;
					for (int a = 0; a < cntx; a++)
					{
						numpast = num;
						i--;
						if (insideb(i, j, k, xfi, yfi, zfi)) {
							num = set_atom(lattice, atoms, mpoint, datapoint, pt, num, i, j, k, halflength);
							cntP = 1;
							if(numpast != num)subcnt++;
						}
					}
				}
				else {
					cntx++;
					for (int a = 0; a < cntx; a++)
					{
						numpast = num;
						i++;
						if (insideb(i, j, k, xfi, yfi, zfi)) {
							num = set_atom(lattice, atoms, mpoint, datapoint, pt, num, i, j, k, halflength);
							cntP = 1;
							if(numpast != num)subcnt++;
						}
					}
				}
				//printf("Sub count x: %d\n", subcnt);
			}
			else {
				if (even(cnty)) {
					cnty++;
					for (int a = 0; a < cnty; a++)
					{
						numpast = num;
						j--;
						if (insideb(i, j, k, xfi, yfi, zfi)) {
							num = set_atom(lattice, atoms, mpoint, datapoint, pt, num, i, j, k, halflength);
							cntP = 1;
							if(numpast != num)subcnt++;
						}
					}
				}
				else {
					cnty++;
					for (int a = 0; a < cnty; a++)
					{
						numpast = num;
						j++;
						if (insideb(i, j, k, xfi, yfi, zfi)) {
							num = set_atom(lattice, atoms, mpoint, datapoint, pt, num, i, j, k, halflength);
							cntP = 1;
							if(numpast != num)subcnt++;
						}
					}
				}
				//printf("Sub count y: %d\n", subcnt);
			}

			cnt += cntP;
			cntl++;
			subcell = subcnt;
		}
		//printf("subcell %d, %d: %d\n", c, k, subcell );
	}

	Vregister(atoms, mpoint, num, pt);
	
}
