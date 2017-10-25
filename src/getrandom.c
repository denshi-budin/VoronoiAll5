#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atom.h"

//#define randseed 25634

/*randomは範囲を決められる*/
//double GetRandom(double min, double max)
//{
//	return min + ((double)rand() / ((double)RAND_MAX + 1)) * (max - min);
//}

double GetRandom(double min, double max)
{
	int subseed = randseed % RAND_MAX;
	int vR = rand() % (subseed + 5);
	double pr = (double)vR/(double)(subseed + 5);
	return min + (pr * (max - min));
}

int Size1Box(position_t *size, double *origin){
	int boxtype = 0, sort = 0;
	double ogap = 0.0;
	position_t newLL = {0.,0.,0.};

	if (pointVor == 1){
		boxtype = 1;
		ogap = 0.0;
	}
	else {
		printf("Choose how to create grid\n[0] No grid\n[1] Equal grid\n[2] Slope grid\n");
		scanf("%d", &boxtype);
		printf("%g | %g | %g \n", latticeLength.x, latticeLength.y , latticeLength.z);
		printf("Input border gap [Angstrom]: ");
		scanf("%lf", &ogap);
	}

	newLL.x = latticeLength.x - (ogap * 2 * ao);
	newLL.y = latticeLength.y - (ogap * 2 * ao);
	newLL.z = latticeLength.z - (ogap * 2 * ao);

	// printf("%g | %g | %g \n", newLL.x, newLL.y , newLL.z);
	if ((newLL.x < 0)||(newLL.y < 0)||(newLL.z < 0)) fail("Outside gap is too big");

	for (int i = 0; i < pointVor; ++i)
	{
		// if (ADDITION < 0) {
		// 	size[i].x = ((newLL.x + (2 * addR.x * latticeConst.x)) / pointVor) * (i+1);
		// 	size[i].y = ((newLL.y + (2 * addR.y * latticeConst.y)) / pointVor) * (i+1);
		// 	size[i].z = ((newLL.z + (2 * addR.z * latticeConst.z)) / pointVor) * (i+1);
		// }
		// else {
			size[i].x = (ogap * ao) + (newLL.x / pointVor) * (i+1) ;
			size[i].y = (ogap * ao) + (newLL.y / pointVor) * (i+1) ;
			size[i].z = (ogap * ao) + (newLL.z / pointVor) * (i+1) ;

			//printf("S%dBoxX %g | S%dBoxY %g | S%dBoxZ %g\n", i, size[i].x/latticeConst.x, i, size[i].y/latticeConst.y, i, size[i].z/latticeConst.z);
		// }
	}

	if (boxtype == 0){
		for (int i = 0; i < pointVor; ++i)
		{
			size[i].x = (ogap * ao) + newLL.x;
			size[i].y = (ogap * ao) + newLL.y;
			size[i].z = (ogap * ao) + newLL.z;
		}
	}
	else if (boxtype == 1);
	else if (boxtype == 2){
		char axis[10];
		char *X, *Y , *Z;
		printf("Choose axis x, y or z: ");
		scanf("%s", axis);
		printf("[1]L to H\n[2]H to L\n");
		scanf("%d", &sort);
		double ab, bc;
		X = strchr(axis, 'x');
		Y = strchr(axis, 'y');
		Z = strchr(axis, 'z');
		if (X != NULL){
			ab = newLL.x/(double)(pointVor*pointVor);
			bc = newLL.x/((double)pointVor/2.0);
			if (sort == 1){
				for (int i = 0; i < pointVor; ++i) size[i].x = (ab*(i+1)*(i+1)) ;
			}
			else {
				for (int i = 0; i < pointVor; ++i) size[i].x = -(ab*(i+1)*(i+1)) + (bc*(i+1));
			}
		}

		if (Y != NULL){
			ab = newLL.y/(double)(pointVor*pointVor);
			bc = newLL.y/((double)pointVor/2.0);
			if (sort == 1){
				for (int i = 0; i < pointVor; ++i) size[i].y = (ab*(i+1)*(i+1)) ;
			}
			else {
				for (int i = 0; i < pointVor; ++i) size[i].y = -(ab*(i+1)*(i+1)) + (bc*(i+1)); 
			}
		}

		if (Z != NULL){
			ab = newLL.z/(double)(pointVor*pointVor);
			bc = newLL.z/((double)pointVor/2.0);
			if (sort == 1){
				for (int i = 0; i < pointVor; ++i) size[i].z = (ab*(i+1)*(i+1)) ;
			}
			else {
				for (int i = 0; i < pointVor; ++i) size[i].z = -(ab*(i+1)*(i+1)) + (bc*(i+1));
			}
		}
	}
	else printf("Used default grid : Equal grid\n");

	(*origin) = (ogap*ao);
	return boxtype;
}

void SetupBox()
{
	double boxOri = 0;
	typbox = 1;
	FILE *fb;
	if((fb = fopen("box.csv", "w")) == NULL) {
		fail("Box file open error!");
	}
	int numbBox = pointVor * pointVor * pointVor;
	smallBox = (box_t *) calloc(numbBox, sizeof(box_t));
	//position_t AddRatio;
	// AddRatio.x = (latticeSizeX * (double)ADDITION) / (double)SizeRotateBox;
	// AddRatio.y = (latticeSizeY * (double)ADDITION) / (double)SizeRotateBox;
	// AddRatio.z = (latticeSizeZ * (double)ADDITION) / (double)SizeRotateBox;

	position_t *sizePerBox;
	sizePerBox = (position_t *) calloc (pointVor, sizeof(position_t));
	typbox = Size1Box(sizePerBox, &boxOri);
	for (int i = 0; i < pointVor; ++i)
	{
		fprintf(fb, "S%dBoxX,%g,S%dBoxY,%g,S%dBoxZ,%g\n", i, sizePerBox[i].x/latticeConst.x, i, sizePerBox[i].y/latticeConst.y, i, sizePerBox[i].z/latticeConst.z);
	}
	fprintf(fb, "\n");

	int nbx = 0;
	for (int i = 0; i < pointVor; ++i)
	{
		for (int j = 0; j < pointVor; ++j)
		{
			for (int k = 0; k < pointVor; ++k)
			{
				// if (ADDITION < 0) {
				// 	if (k == 0) smallBox[nbx].minBox.x =  - (AddRatio.x * latticeConst.x);
				// 	else smallBox[nbx].minBox.x = sizePerBox[k - 1].x - (AddRatio.x * latticeConst.x);
				// 	if (j == 0) smallBox[nbx].minBox.y =  - (AddRatio.y * latticeConst.y);
				// 	else smallBox[nbx].minBox.y = sizePerBox[j - 1].y - (AddRatio.y * latticeConst.y);
				// 	if (i == 0) smallBox[nbx].minBox.z =  - (AddRatio.z * latticeConst.z);
				// 	else smallBox[nbx].minBox.z = sizePerBox[i - 1].z - (AddRatio.z * latticeConst.z);
				// 	smallBox[nbx].maxBox.x =  sizePerBox[k].x - (AddRatio.x * latticeConst.x);
				// 	smallBox[nbx].maxBox.y =  sizePerBox[j].y - (AddRatio.y * latticeConst.y);
				// 	smallBox[nbx].maxBox.z =  sizePerBox[i].z - (AddRatio.z * latticeConst.z);
				// }
				// else {
					if ((k == 0) || (typbox == 0)) smallBox[nbx].minBox.x = boxOri;
					else smallBox[nbx].minBox.x = sizePerBox[k - 1].x;
					if ((j == 0) || (typbox == 0)) smallBox[nbx].minBox.y = boxOri;
					else smallBox[nbx].minBox.y = sizePerBox[j - 1].y;
					if ((i == 0) || (typbox == 0)) smallBox[nbx].minBox.z = boxOri;
					else smallBox[nbx].minBox.z = sizePerBox[i - 1].z;
					smallBox[nbx].maxBox.x = sizePerBox[k].x ;
					smallBox[nbx].maxBox.y = sizePerBox[j].y ;
					smallBox[nbx].maxBox.z = sizePerBox[i].z ;
				// }
				smallBox[nbx].used = 0;
				fprintf(fb, "MIN%d,%g,%g,%g\n", nbx, smallBox[nbx].minBox.x/latticeConst.x, smallBox[nbx].minBox.y/latticeConst.y, smallBox[nbx].minBox.z/latticeConst.z);
				fprintf(fb, "MAX%d,%g,%g,%g\n", nbx, smallBox[nbx].maxBox.x/latticeConst.x, smallBox[nbx].maxBox.y/latticeConst.y, smallBox[nbx].maxBox.z/latticeConst.z);
				nbx++;
			}
		}
	}
	fclose(fb);
}

#if 0
double Random(int magic){
	static long int previous = 4771;//4771
	int multi = magic;//16897
	static int add = 0;
	static long int divide = 2147483647;

	double cardinal = (double)previous;
	double factor = (double)multi;
	double increment = (double)add;
	double modulus = (double)divide;
	cardinal = fmod(cardinal * factor + increment, modulus);
	previous = (long int)cardinal;
	return cardinal / modulus;
}

double GetRandom(double min, double max)
{
	return ((max * Random(17643)) + min);
}
#endif
