#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atom.h"

void lower_string(char s[]) {
   int c = 0;
   while (s[c] != '\0') {
      if (s[c] >= 'A' && s[c] <= 'Z') {
         s[c] = s[c] + 32;
      }
      c++;
   }
}

int readpointdata(){
	FILE *fp, *ls;
	position_t *pointr;
	dpoint_t *dtr;
	char readdir[100];
	strcpy(readdir, "");
	system("ls > list.txt");
	if((ls = fopen("list.txt", "r")) != NULL){
		int folder = 0;
		while(fscanf(ls, "%s", readdir) != EOF){
			if(strcmp(readdir, "output") == 0) {
				folder = 1;
				break;
			}
			strcpy(readdir, "");
		}
		if(folder == 0) {
			system("mkdir output");
			system("mkdir output/back");
			fclose(ls);
			system ("rm list.txt");
			return 0;
		}
	}
	fclose(ls);
	system ("rm list.txt");
	

	remake = 0;
	if((fp = fopen("output/VoronoiPoint.xyz", "r")) == NULL){
		//printf("voronoi point file open error!\n");
		return 0;
	}
	else {
		printf("Remake from voronoi point? 1:yes 0:No : ");
		scanf("%d", &remake);

		if (remake == 0) {
			fclose(fp);
			return 0;
		}
	}

	int nP = 0, rP = 0;
	char str[20];
	pointVor = 0;
	fscanf(fp, "%d", &nP);
	printf("Total point in file: %d\n", nP);
	pointr = (position_t *) calloc (nP+1, sizeof(position_t));
    dtr = (dpoint_t *) calloc (nP+1, sizeof(dpoint_t));

	fscanf(fp, "%s %lf %s %s %s %lf %s %s %s %lf %s %s %d", 
		str, &latticeLength.x, str, str, str, &latticeLength.y, str, str, str, &latticeLength.z, str, str, &randseed);

	//printf("latticeLength: %g, %g, %g \n", latticeLength.x, latticeLength.y, latticeLength.z);

	latticeLength.x *= ao;
	latticeLength.y *= ao;
	latticeLength.z *= ao;

	latticeSizeX = (int)(latticeLength.x / latticeConst.x);
	latticeSizeY = (int)(latticeLength.y / latticeConst.y);
	latticeSizeZ = (int)(latticeLength.z / latticeConst.z);
	printf("latticeSize: %d, %d, %d \n", latticeSizeX, latticeSizeY, latticeSizeZ);

	for (int i = 0; i < nP; ++i)
	{
		fscanf(fp, "%s %lf %lf %lf %s %lf %lf %lf", 
			       str, &pointr[i].x, &pointr[i].y, &pointr[i].z, str, &dtr[i].psi, &dtr[i].theta, &dtr[i].phi);
		pointr[i].x *= ao;
		pointr[i].y *= ao;
		pointr[i].z *= ao;
		if (inSide(pointr[i], simboxori, latticeLength)){
			rP += 1;
			//printf("[%d] %lf %lf %lf %lf %lf %lf\n", i, pointr[i].x/ao, pointr[i].y/ao, pointr[i].z/ao, dtr[i].psi, dtr[i].theta, dtr[i].phi);
		}
	}

	printf("Number of ori point: %d\n", rP);
	pointVor = rP;
	mN = pointVor * 8;

	mpoint = (mpoint_t *) calloc (mN+1, sizeof(mpoint_t));
    datapoint = (dpoint_t *) calloc (mN+1, sizeof(dpoint_t));

    for (int i = 0; i < rP; ++i)
    {
    	mpoint[i].pos = pointr[i];
    	datapoint[i].oriPoint = -1;
    	datapoint[i].psi = dtr[i].psi;
    	datapoint[i].theta = dtr[i].theta;
    	datapoint[i].phi = dtr[i].phi;
    }

    fclose(fp);
    return 1;

}

void config(){
	FILE *cfg;
	char strg[50];
	if((cfg = fopen("Configuration.cfg", "r")) != NULL){
		fscanf(cfg, " %s %s", strg, lattice);
		lower_string(lattice);
		fscanf(cfg, "%s %lf", strg, &a0);
		fscanf(cfg, "%s %s", strg, atom1);
		//fscanf(cfg, "%s %d", strg, &ParticleNum);
		fscanf(cfg, "%s %lf", strg, &massa);
		fclose(cfg);
	}
	else {
		if((cfg = fopen("Configuration.cfg", "w")) == NULL){
			fail("Configuration file error");
		}
		printf("Enter lattice structure: ");
		scanf("%s", lattice);
		fprintf(cfg, "Lattice: %s\n", lattice);
		lower_string(lattice);
		printf("Size per 1 unit lattice: ");
		scanf("%lf", &a0);
		fprintf(cfg, "Scale: %g\n", a0);
		printf("Element name: ");
		scanf("%s", atom1);
		fprintf(cfg, "Element: %s\n", atom1);
		// printf("ParticalNum: ");
		// scanf("%d", &ParticleNum);
		// fprintf(cfg, "ParticalNum: %d\n", ParticleNum);
		printf("Mass: ");
		scanf("%lf", &massa);
		fprintf(cfg, "Mass: %g\n", massa);
		fclose(cfg);
	}

	if (strcmp("fcc", lattice) == 0) ParticleNum = 4;
	else if (strcmp("bcc", lattice) == 0) ParticleNum = 2;
	else if (strcmp("diamond", lattice) == 0) ParticleNum = 8;
	else if (strcmp("hcp", lattice) == 0) ParticleNum = 4;

	if(strcmp("hcp", lattice) == 0) {
		latticeConst.x = a0 * ao;
		latticeConst.y = a0 * sqrt(3.) * ao;
		latticeConst.z = a0 * 2. * sqrt(2.) / sqrt(3.) * ao;
	}
	else {
		latticeConst.x = a0 * ao;
		latticeConst.y = a0 * ao;
		latticeConst.z = a0 * ao;
	}
}