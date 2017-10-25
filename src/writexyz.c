#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "atom.h"

void writeSample(position_t data[], int sumVnum, const char *nmF){
	char target[100];
	FILE *fw;
	fpos_t positionf;
	int sampleSum = 0;
	strcpy(target,"output/");
	strcat(target,nmF);
	if((fw = fopen(target, "w")) == NULL) {
		fail("output file create error!");
		return;
	}
	fgetpos(fw, &positionf);
	fprintf(fw, "          \n");
	fprintf(fw, "Lattice=\" %f 0.0 0.0 0.0 %f 0.0 0.0 0.0 %f \"\n", latticeLength.x * 1.0e10, latticeLength.y * 1.0e10, latticeLength.z * 1.0e10);

	for (int i = 0; i < sumVnum; i++) {
		if (inSide(data[i], simboxori, latticeLength)) {
				fprintf(fw, "%s\t%f\t%f\t%f\n", atom1, data[i].x * 1.0e10, data[i].y * 1.0e10, data[i].z * 1.0e10);
				sampleSum++;
			}
	}
	fsetpos(fw, &positionf);
	fprintf(fw, "%d", sampleSum);
	fclose(fw);
	printf("Number of atoms in %s: %d\n", nmF, sampleSum);
}

void writeLammps(mpoint_t mpoint[], int sumVnum){
	FILE *lm;
	fpos_t positionf;
	char fname[200];
	int sampleSum = 0;
	sprintf(fname, "output/%sin.mddata", atom1);
	if((lm = fopen(fname, "w")) == NULL) {
		fail("output lammps create error!");
		return;
	}

	fprintf(lm, " Start File for LAMMPS Voronoi %s; S:%d%d%d vN:%d g:%f r:%f asV:%g rand:%d \n\n",
			atom1, latticeSizeX, latticeSizeY, latticeSizeZ, pointVor, gap, nearAtom, avgsizeV/ao, randseed);
	fgetpos(lm, &positionf);
	fprintf(lm, "                   \n");
	fprintf(lm, "\n 1 atom types\n");
	fprintf(lm, "\n 0.0  %lf xlo xhi", latticeLength.x * 1.0e10);
	fprintf(lm, "\n 0.0  %lf ylo yhi", latticeLength.y * 1.0e10);
	fprintf(lm, "\n 0.0  %lf zlo zhi\n", latticeLength.z * 1.0e10);

	fprintf(lm, "\n Masses\n\n 1  %g\n", massa);
	fprintf(lm, "\n Atoms\n");

	int count = 1;
	for (int i = 0; i < sumVnum; i++) {
		if (inSide(sample[i], simboxori, latticeLength)) {
				fprintf(lm, "\n %d\t1\t%f\t%f\t%f", count++, sample[i].x * 1.0e10, sample[i].y * 1.0e10, sample[i].z * 1.0e10);
				sampleSum++;
			}
	}
	fsetpos(lm, &positionf);
	fprintf(lm, " %d atoms", sampleSum);
	fclose(lm);
}

void writePoint(mpoint_t mpoint[], dpoint_t datapoint[]){
	FILE *fv;
	fpos_t positionf;
	int pt = 0;
	if((fv = fopen("output/VoronoiPoint.xyz", "w")) == NULL){
		fail("voronoi point create error!");
		return;
	}
	fgetpos(fv, &positionf);
	fprintf(fv, "              \n");
	fprintf(fv, "Lattice=\" %f 0.0 0.0 0.0 %f 0.0 0.0 0.0 %f \" rand= %d\n", 
		latticeLength.x * 1.0e10, latticeLength.y * 1.0e10, latticeLength.z * 1.0e10, randseed);

	for (int i = 0; i < mN; i++) {
		if (!inSide(mpoint[i].pos, latorigin, latticeLength)) continue;
		pt++;
		char vnum[10] = {'\0'};
		if (datapoint[i].oriPoint == -1)
			sprintf(vnum, "V%d", i);
		else sprintf(vnum, "V%d", datapoint[i].oriPoint);
		fprintf(fv, "%s\t%f\t%f\t%f\t0.5\t%g\t%g\t%g\n", vnum, mpoint[i].pos.x * 1.0e10, mpoint[i].pos.y * 1.0e10, mpoint[i].pos.z * 1.0e10, 
															datapoint[i].psi, datapoint[i].theta, datapoint[i].phi);
	}
	fsetpos(fv, &positionf);
	fprintf(fv, "%d", pt);
	fclose(fv);
}

void writeVoronoi(mpoint_t mpoint[], dpoint_t datapoint[], int sumVnum){
	FILE *fp;
	if((fp = fopen("output/VoronoiRegion.xyz", "w")) == NULL){
		fail("voronoi region create error!");
		return;
	}

	fprintf(fp, "%d\n", sumVnum);
	fprintf(fp, "Lattice=\" %f 0.0 0.0 0.0 %f 0.0 0.0 0.0 %f \"\n", latticeLength.x * 1.0e10, latticeLength.y * 1.0e10, latticeLength.z * 1.0e10);

	for(int i = 0; i < mN; i++){
		if (!inSide(mpoint[i].pos, latorigin, latticeLength)) continue;
		for(int j=0;j<mpoint[i].voronoiNum;j++){
			char num[10] = {'\0'};
			if (datapoint[i].oriPoint == -1)
				sprintf(num, "N%d", i);
			else sprintf(num, "N%d", datapoint[i].oriPoint);
			fprintf(fp, "%s\t%f\t%f\t%f\t0.3\t%f\n", num, mpoint[i].voronoi[j].x * 1.0e10, mpoint[i].voronoi[j].y * 1.0e10, mpoint[i].voronoi[j].z * 1.0e10,
													mpoint[i].dis[j]);
		}
	}
	fclose(fp);
}

void writeFile(mpoint_t mpoint[], dpoint_t datapoint[]){
	FILE *ls;
	char readdir[100];
	strcpy(readdir, "");
	system("ls output > list.txt");
	if((ls = fopen("list.txt", "r")) != NULL){
		int folder = 0;
		while(fscanf(ls, "%s", readdir) != EOF){
			if(strcmp(readdir, "back") == 0) {
				folder = 1;
				break;
			}
			strcpy(readdir, "");
		}
		if(folder == 0) {
			system("mkdir output/back");
		}
	}
	fclose(ls);
	system ("rm list.txt");

	printf("Write to file...\n");
	system("mv output/*.xyz output/back");
	system("mv output/*.mddata output/back");

	int totalVor = 0;
	for(int i=0;i<mN;i++){
		if (!inSide(mpoint[i].pos, latorigin, latticeLength)) continue;
		totalVor += mpoint[i].voronoiNum;
	}
	printf("Expected number of atoms: %d\n", N);
	printf("Total number of atoms: %d\n", totalVor);
	printf("Average from %d grains size: %g\n", pointVor, avgsizeV/ao);

	writePoint(mpoint, datapoint);
	writeVoronoi(mpoint, datapoint, totalVor);
	writeSample(sample, totalVor, "sample.xyz");
	writeSample(nearAtm, NA, "nearatom.xyz");
	writeLammps(mpoint, totalVor);
}
