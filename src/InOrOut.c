#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atom.h"

int NA = 0;
int room[2][2][2];

void settingParameter(){
	config();
	pointVor = 0;
	PLUS = 2;
	int rm = 0;
	avgsizeV = 0.0;
	char prd[10];
	char *X, *Y, *Z;
	simboxori.x = -0.0;
	simboxori.y = -0.0;
	simboxori.z = -0.0;
	for (int i = 0; i < 2; ++i)
		for (int j = 0; j < 2; ++j)
			for (int k = 0; k < 2; ++k) 
				room[i][j][k] = ++rm;


	if (readpointdata() == 0){
		printf("Insert lattice size [Integer Positive Number]:-\n");
		printf("X: ");
		scanf("%d", &latticeSizeX);
		printf("Y: ");
		scanf("%d", &latticeSizeY);
		printf("Z: ");
		scanf("%d", &latticeSizeZ);

		if ((latticeSizeX < 0)||(latticeSizeY < 0)||(latticeSizeZ < 0)) fail("Size less than 0");

		latticeLength.x = latticeConst.x * latticeSizeX;
		latticeLength.y = latticeConst.y * latticeSizeY;
		latticeLength.z = latticeConst.z * latticeSizeZ;

		/*ボロノイ母点*/
		while (pointVor == 0){
			printf("Insert number of voronoi point: ");
			scanf("%d", &pointVor);
		}
		/*全ボロノイ母点*/
		mN = pointVor * 8;
	}

	if (pointVor == 1) PLUS = 0;
	else PLUS = 2;
	
	N = ParticleNum * latticeSizeX * latticeSizeY * latticeSizeZ;
	alloribox.x = - (latticeLength.x/2);
	alloribox.y = - (latticeLength.y/2);
	alloribox.z = - (latticeLength.z/2);

	if(latticeLength.x >= latticeLength.y){
		if(latticeLength.x >= latticeLength.z) {SizeRotateBox = latticeSizeX; latConstBox = latticeConst.x;}
		else {SizeRotateBox = latticeSizeZ; latConstBox = latticeConst.z;}
	}
	else {
		if (latticeLength.y >= latticeLength.z) {SizeRotateBox = latticeSizeY; latConstBox = latticeConst.y;}
		else {SizeRotateBox = latticeSizeZ; latConstBox = latticeConst.z;}
	}

	if(latticeSizeZ == 1) SizeRotateBoxZ = latticeSizeZ + 1;
	else SizeRotateBoxZ = SizeRotateBox;

	Rbox = sqrt(pow(SizeRotateBox*latConstBox,2)+pow(SizeRotateBox*latConstBox,2));
	RNum = ParticleNum * (SizeRotateBox+PLUS) * (SizeRotateBox+PLUS) * (SizeRotateBoxZ+PLUS);

	/* ボロノイ設定(分割数) */
	/*ボロノイとボロノイの隙*/
	if (pointVor == 1){
		gap = 0.0;
		nearAtom = 0.0;
		//strcpy(prd, "q");
	}
	else {
		printf("Insert gap size between voronoi: ");
		scanf("%f", &gap);
		printf("Insert near atom remover range: ");
		scanf("%f", &nearAtom);
		printf("Add periodic boundaries x, y, z or q for no boundaries: ");
		scanf("%s", prd);
	}


	X = strchr(prd, 'x');
	Y = strchr(prd, 'y');
	Z = strchr(prd, 'z');
	latorigin.x = 0.0;
	latorigin.y = 0.0;
	latorigin.z = 0.0;

	if (X != NULL) latorigin.x = alloribox.x;
	if (Y != NULL) latorigin.y = alloribox.y;
	if ((Z != NULL) && (latticeSizeZ != 1)) latorigin.z = alloribox.z; //

	if (strchr(prd,'q') != NULL) printf("No periodic boundaries.\n");
	else if (pointVor == 1);
	else if ((X == NULL) && (Y == NULL) && (Z == NULL)) printf("No periodic boundaries selected.\n");

	FILE *log;
	if((log = fopen("output/Size.log", "w")) == NULL) fail("Unable to create Size.log file");
	fprintf(log, "Output for every points: \n");
	fclose(log); 

}

double distance(position_t a, position_t b){
	if (latticeSizeZ != 1) return sqrt(pow(a.x-b.x,2.0)+pow(a.y-b.y,2.0)+pow(a.z-b.z,2.0));
	else return sqrt(pow(a.x-b.x,2.0)+pow(a.y-b.y,2.0));
}

void RandomBox(mpoint_t mpoint[], dpoint_t datapoint[]) //
{
	if (remake >= 1) return;
	int cnt = 0;
	int timeout = 0;
	int inputval = 1;
	SetupBox();
	int numbBox = pointVor * pointVor * pointVor;
	position_t *posp;
	posp = (position_t *) calloc (pointVor, sizeof(position_t));

	for (size_t i = 0; i < mN; i++) {
		datapoint[i].oriPoint = 0;
		datapoint[i].psi = 0;//x30
		datapoint[i].theta = 0;//z55
		datapoint[i].phi = 0;//y45
	}

	if (pointVor == 1) {
		mpoint[cnt].pos.x = (((float)latticeSizeX / 2.0) * latticeConst.x);//latticeLength.x / 2;
		mpoint[cnt].pos.y = (((float)latticeSizeY / 2.0) * latticeConst.y);//latticeLength.y / 2;
		if (latticeSizeZ > 1) mpoint[cnt].pos.z = (((float)latticeSizeZ / 2.0) * latticeConst.z);//latticeLength.z / 2;
		else mpoint[cnt].pos.z = 0;//(latticeLength.z)/2;
		mpoint[cnt].voronoiNum = 0;
		posp[cnt] = mpoint[cnt].pos;
		//printf("Numb box for point %d : 1\n", cnt + 1);
		cnt++;
	}

	if (typbox == 0){
		printf("[0] Input coordinate or [1] Random coordinate : ");
		scanf("%d", &inputval);
	}

	if (inputval == 0){
		printf("Simulation box size [Angstrom]:\n");
		printf("x lo: %lf  x hi: %lf \n", simboxori.x/ao, latticeLength.x/ao);
		printf("y lo: %lf  y hi: %lf \n", simboxori.y/ao, latticeLength.y/ao);
		printf("z lo: %lf  z hi: %lf \n", simboxori.z/ao, latticeLength.z/ao);
		printf("\n");
		while (cnt < pointVor){
			printf("Coordinate position %d x: ", cnt+1);
			scanf("%lf", &mpoint[cnt].pos.x);
			mpoint[cnt].pos.x *= ao;
			printf("Coordinate position %d y: ", cnt+1);
			scanf("%lf", &mpoint[cnt].pos.y);
			mpoint[cnt].pos.y *= ao;
			printf("Coordinate position %d z: ", cnt+1);
			scanf("%lf", &mpoint[cnt].pos.z);
			mpoint[cnt].pos.z *= ao;
			posp[cnt] = mpoint[cnt].pos;
			mpoint[cnt].voronoiNum = 0;
			datapoint[cnt].oriPoint = -1;
			cnt++;
		}
	}
	else {
		while(cnt < pointVor){
			timeout++;
			int noBox = rand() % numbBox;
			if(smallBox[noBox].used > 0) continue;
			else {
				position_t countN;
				countN.x = noBox % pointVor;
				countN.y = floor((noBox % (pointVor*pointVor)) / pointVor);
				countN.z = floor(noBox / (pointVor*pointVor));

				for (int j = 0; j < numbBox; j++) {
					if ((j % pointVor == countN.x) ||
						(floor((j % (pointVor*pointVor)) / pointVor) == countN.y) ||
						(floor(j / (pointVor*pointVor)) == countN.z))
					{
						smallBox[j].used = 2;
					}
				}

				mpoint[cnt].pos.x = GetRandom(smallBox[noBox].minBox.x, smallBox[noBox].maxBox.x);
				mpoint[cnt].pos.y = GetRandom(smallBox[noBox].minBox.y, smallBox[noBox].maxBox.y);
				if (latticeSizeZ > 1) mpoint[cnt].pos.z = GetRandom(smallBox[noBox].minBox.z, smallBox[noBox].maxBox.z);
				else 
					mpoint[cnt].pos.z = (latticeLength.z)*0.25;
				posp[cnt] = mpoint[cnt].pos;
				mpoint[cnt].voronoiNum = 0;
				datapoint[cnt].oriPoint = -1;
				smallBox[noBox].used = 1;
				printf("Numb box for point %d : %d\t Time count : %d\n", cnt + 1, noBox, timeout);
				cnt++;
			}
		}
	}
	free(smallBox);

	double *dispt;
	dispt = (double *) calloc (pointVor, sizeof(double));
	for (int i = 0; i < pointVor; ++i)
	{
		dispt[i] = distance(simboxori, posp[i]);
	}
	swaping_small_big(posp, dispt, pointVor, "rePoint", ao);
	for (int i = 0; i < pointVor; ++i)
	{
		mpoint[i].pos = posp[i];
	}
	free(dispt);
	free(posp);
}

void RandomCheck(mpoint_t mpoint[]){
	FILE *fc;
	fc = fopen("randcheck.csv", "w");
	position_t def;
	double range = 0.0;
	int row = 0;

	for (int i = 0; i < (mN-1); i++) {
		for (int j = (i+1); j < mN; j++) {
			def.x = mpoint[i].pos.x - mpoint[j].pos.x;
			def.y = mpoint[i].pos.y - mpoint[j].pos.y;
			def.z = mpoint[i].pos.z - mpoint[j].pos.z;

			range = sqrt((def.x * def.x) + (def.y * def.y) + (def.z * def.z));
			fprintf(fc, "V%d,V%d,%g\n", i, j, range * 1.0e10);
			row++;
		}
	}
	fprintf(fc, ",Sum,=SUM(C1:C%d)\n", row);
	fprintf(fc, ",Avg,=AVERAGE(C1:C%d)\n", row);
	fclose(fc);
}

int ainOrOut(position_t atom, mpoint_t mpoint[], double *disv, int pt){

	//if(!inSide(atom, simboxori, latticeLength)) return 0;
	if (pointVor == 1) return 1;
	if(!inSide(mpoint[pt].pos, simboxori, latticeLength)) {
		if(!inSide(atom, simboxori, latticeLength)) return 0;
	}
	else {
		if(!inSide(atom, alloribox, latticeLength)) return 0;
	}

	double disCent;
	double disPoint, nearedg;
	double diss = 99999.;

	disCent = distance(atom, mpoint[pt].pos);
	for (int i = 0; i < mN; ++i)
	{
		if (i == pt) continue;
		if (!inSide(mpoint[i].pos, latorigin, latticeLength)) continue;
		if (far(mpoint[pt].pos, mpoint[i].pos)) continue;
		
		disPoint = distance(atom, mpoint[i].pos);
		if (disPoint <= (disCent + (gap*ao))) return 0;
		nearedg = (disPoint - disCent)/ao;
		if ((diss == 99999.)||(diss > nearedg)) diss = nearedg;

		if ((mpoint[i].voronoiNum != 0) && (nearedg < (nearAtom*2))) {
			for (int j = 0; j < mpoint[i].voronoiNum; ++j)
			{
				if(mpoint[i].dis[j] > (nearAtom*2)) continue;
				if(distance(atom, mpoint[i].voronoi[j]) < (nearAtom*ao)){
					if(inSide(mpoint[i].voronoi[j], simboxori,latticeLength)){
						nearAtm[NA] = mpoint[i].voronoi[j];
						NA++;
					}
					if(inSide(atom, simboxori, latticeLength)){
						nearAtm[NA] = atom;
						NA++;
					}
					if (NA >= N*2) fail("Too much atoms been thrown out");
					return 0;
				}
			}
		}
	}
	(*disv) = diss;

	return 1;
}

/*もう一度原子と原子の間の距離を確認する。この関数を使わなくてもいい。*/
void checkrange(mpoint_t mpoint[]){

	double range;
	position_t delta;
	int numb = 0;
	FILE *fr;
	fr = fopen("range.csv", "w");
	for (int s = 1; s < mN; s++) {
		if (!inSide(mpoint[s].pos, latorigin, latticeLength)) continue;
		for (int f = 0; f < mpoint[s].voronoiNum; f++) {
			for (int g = 0; g <= s-1 ; g++){
				for (int h = 0; h < mpoint[g].voronoiNum; h++){
					delta.x = (mpoint[g].voronoi[h].x - mpoint[s].voronoi[f].x) * (mpoint[g].voronoi[h].x - mpoint[s].voronoi[f].x);
					delta.y = (mpoint[g].voronoi[h].y - mpoint[s].voronoi[f].y) * (mpoint[g].voronoi[h].y - mpoint[s].voronoi[f].y);
					delta.z = (mpoint[g].voronoi[h].z - mpoint[s].voronoi[f].z) * (mpoint[g].voronoi[h].z - mpoint[s].voronoi[f].z);
					range = (sqrt(delta.x + delta.y + delta.z))/ao;
					if (range < nearAtom){
						fprintf(fr, "%d,N%d,N%d,A%d,A%d,%g,%g,%g,%g,%g,%g,%g\n", ++numb, s, g, f, h,
						mpoint[s].voronoi[f].x, mpoint[s].voronoi[f].y, mpoint[s].voronoi[f].z,
						mpoint[g].voronoi[h].x, mpoint[g].voronoi[h].x, mpoint[g].voronoi[h].x, range);
					}
				}
			}
		}
	}
	if (numb == 0) fprintf(fr, "No merge atom\n");
	fclose (fr);
}

/*バラバラの原子を小さい順からまとめるすること*/
void relocate(mpoint_t mpoint[]){
	printf("Relocating atoms...");
	int sumVnum = 0;

	for(int i=0;i<mN;i++){
		sumVnum += mpoint[i].voronoiNum;
	}

	sample = (position_t *) calloc(sumVnum, sizeof(position_t));

	int k = 0;
	for (int i = 0; i < mN; i++) {
		for (size_t j = 0; j < mpoint[i].voronoiNum; j++) {
			sample[k].x = mpoint[i].voronoi[j].x;
			sample[k].y = mpoint[i].voronoi[j].y;
			sample[k].z = mpoint[i].voronoi[j].z;
			k++;
		}
	}
	printf("Done!\n");

}

position_t addcopy(position_t orig, position_t addL, int ax, int ay, int az){
	position_t npt;
	npt.x = orig.x + (addL.x * ax);
	npt.y = orig.y + (addL.y * ay);
	npt.z = orig.z + (addL.z * az);
	return npt;
}

int positive(double n){
	if (n >= 0) return 1;
	else return 0;
}

void copyPoint(mpoint_t mpoint[], dpoint_t datapoint[]){  //
	int addto = 0;
	if (mN == pointVor) return;

	for (int i = -1; i <= 1; ++i)
	{
		for (int j = -1; j <= 1; ++j)
		{
			for (int k = -1; k <= 1; ++k)
			{
				if ((i == 0) && (j == 0) && (k == 0)) continue;
				for (int l = 0; l < pointVor; ++l)
				{
					if ((pointVor + addto) > mN) fail ("Point overflow.");
					mpoint[pointVor + addto].pos = addcopy(mpoint[l].pos, latticeLength, k, j, i);
					if (inSide(mpoint[pointVor + addto].pos, alloribox, latticeLength)) {
						datapoint[pointVor + addto].oriPoint = l;
						addto++;
					}
				}
			}
		}
	}

	FILE *fp;
	fp = fopen("Point.csv", "w");
	double R;
	position_t deltaP;
	int clss;

	fprintf(fp, ",");
	for (int i = 0; i < mN; i++) fprintf(fp, "%d,", i+1);
	fprintf(fp, "\n");
	for (int i = 0; i < mN; i++) {
		fprintf(fp, "%d,", i+1);
		if (!inSide(mpoint[i].pos, latorigin, latticeLength)) {
				fprintf(fp, "\n");
				continue;
			}
		for (int j = 0; j < mN; j++) {
			if (!inSide(mpoint[j].pos, latorigin, latticeLength)) {
				fprintf(fp, ",");
				continue;
			}
			deltaP.x = mpoint[j].pos.x - mpoint[i].pos.x;
			deltaP.y = mpoint[j].pos.y - mpoint[i].pos.y;
			deltaP.z = mpoint[j].pos.z - mpoint[i].pos.z;
			R = sqrt(pow(deltaP.x,2) + pow(deltaP.y,2) + pow(deltaP.z,2));
			clss = room[positive(deltaP.x)][positive(deltaP.y)][positive(deltaP.z)];
			if (R/ao == 0) clss = 0;
			fprintf(fp, "%g-%d,",R/ao, clss);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	if (pointVor != 1) printf("Done copy vorinoi point.\n");
}

int far(position_t place, position_t compare)
{
	int toofar = 0;
	double R = 0;
	position_t deltaP;
	deltaP.x = place.x - compare.x;
	deltaP.y = place.y - compare.y;
	deltaP.z = place.z - compare.z;
	R = sqrt(pow(deltaP.x,2) + pow(deltaP.y,2) + pow(deltaP.z,2));
	if(R >= Rbox) toofar = 1;
	else toofar = 0;

	return toofar;
}

double sizeVoronoi(position_t ato[], int total, double map[]){
	FILE *log;
	double sizeV = 0.0;
	double dis = 0.0;
	static int logpoit = 0;
	for (int i = 0; i < total-1; ++i)
	{
		if (map[i] > (nearAtom*3)) continue;
		for (int j = i; j < total; ++j)
		{
			if (map[j] > (nearAtom*3)) continue;
			dis = distance(ato[i],ato[j]);
			if(sizeV < dis) sizeV = dis;
		}
	}
	printf("Size this point is: %g\n", sizeV/ao );
	if ((log = fopen("output/Size.log", "a")) == NULL) fail("Unable to open Size.log file");
	logpoit += 1;
	fprintf(log, "[%d] Size for this point: %g\n", logpoit, sizeV/ao);
	fclose(log);
	return sizeV;
}

int inSide(position_t pt, position_t ori, position_t inR){
	int in = 0;
	if ((pt.x >= ori.x) && (pt.x < (inR.x - ori.x))) in += 1;
	if ((pt.y >= ori.y) && (pt.y < (inR.y - ori.y))) in += 1;
	if ((pt.z >= ori.z) && (pt.z < (inR.z - ori.z))) in += 1;
	
	if (in == 3) return 1;
	else return 0;
}

