#include <stdio.h>
#include <stdlib.h>
#include "atom.h"

void writeFile(mpoint_t mpoint[])
{
	FILE *fw;
	if((fw = fopen("output/sample.xyz", "w")) == NULL) {
		printf("output open error!\n");
		return;
	}
	FILE *fv;
	if((fv = fopen("output/VoronoiPoint.xyz", "w")) == NULL){
		printf("voronoi open error!\n");
		return;
	}
	FILE *fp;
	if((fp = fopen("output/VoronoiRegion.xyz", "w")) == NULL){
		printf("voronoi open error!\n");
		return;
	}
	int sampleSum = 0;
	fpos_t positionf;

	int sumVnum = 0;
	for(int i=0;i<mN;i++){
		sumVnum += mpoint[i].voronoiNum;
	}
	printf("予想原子数 : %d\n", N);
	printf("全原子数：%d\n", sumVnum);

	fprintf(fv, "%d\n", mN);
	fprintf(fv, "Lattice=\"%f 0.0 0.0 0.0 %f 0.0 0.0 0.0 %f\"\n", latticeLengthX * 1.0e10, latticeLengthY * 1.0e10, latticeLengthZ * 1.0e10);
	fprintf(fp, "%d\n", sumVnum);
	fprintf(fp, "Lattice=\"%f 0.0 0.0 0.0 %f 0.0 0.0 0.0 %f\"\n", latticeLengthX * 1.0e10, latticeLengthY * 1.0e10, latticeLengthZ * 1.0e10);
	fgetpos(fw, &positionf);
	fprintf(fw, "          \n");
	fprintf(fw, "Lattice=\"%f 0.0 0.0 0.0 %f 0.0 0.0 0.0 %f\"\n", latticeLengthX * 1.0e10, latticeLengthY * 1.0e10, latticeLengthZ * 1.0e10);
	int aNum = 0;
	for(int i=0;i<mN;i++){
		char vnum[10] = {'\0'};
		sprintf(vnum, "V%d", i);
		fprintf(fv, "%s\t%f\t%f\t%f\t0.7\n", vnum, mpoint[i].pos.x * 1.0e10, mpoint[i].pos.y * 1.0e10, mpoint[i].pos.z * 1.0e10);
		for(int j=0;j<mpoint[i].voronoiNum;j++){
			char num[10] = {'\0'};
			sprintf(num, "N%d", i);
			fprintf(fp, "%s\t%f\t%f\t%f\t0.7\n", num, mpoint[i].voronoi[j].x * 1.0e10, mpoint[i].voronoi[j].y * 1.0e10, mpoint[i].voronoi[j].z * 1.0e10);
			// fprintf(fw, "%s\t%f\t%f\t%f\n", atom1, mpoint[i].voronoi[j].x * 1.0e10, mpoint[i].voronoi[j].y * 1.0e10, mpoint[i].voronoi[j].z * 1.0e10);
			if ((sample[aNum].x >= 0 && sample[aNum].x <= latticeLengthX) &&
					(sample[aNum].y >= 0 && sample[aNum].y <= latticeLengthY) &&
				  (sample[aNum].z >= 0 && sample[aNum].z <= latticeLengthZ)){
						fprintf(fw, "%s\t%f\t%f\t%f\n", atom1, sample[aNum].x * 1.0e10, sample[aNum].y * 1.0e10, sample[aNum].z * 1.0e10);
						sampleSum++;
				}
			aNum++;
		}
	}
	fsetpos(fw, &positionf);
	fprintf(fw, "%d", sampleSum);
	printf("原子数 : %d\n", sampleSum);
	/*for(int i=0;i<N;i++){
		fprintf(fw, "%s\t%f\t%f\t%f\t0.8\n", atom1, atoms[i].x * 1.0e10, atoms[i].y * 1.0e10, atoms[i].z * 1.0e10);
	}*/
}
