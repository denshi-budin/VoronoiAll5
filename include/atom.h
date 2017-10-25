#ifndef _ATOM_H_
#define _ATOM_H_

/* 規定値 */
#define kb 8.617330350e-5
#define ao 1.0e-10
#define PI 3.14159

int latticeSizeX, latticeSizeY, latticeSizeZ, N, RNum, pointVor, mN, SizeRotateBox, SizeRotateBoxZ, randseed, PLUS, remake, typbox;
float gap, nearAtom;
double Rbox, a0, avgsizeV, latConstBox;

//#define ADDITION 0

typedef struct{
	double x;
	double y;
	double z;
}position_t;

typedef struct{
	position_t pos; /* ボロノイ母点の位置 */
	position_t *voronoi; /* ボロノイ領域内の原子の位置 */
	double *dis;
	int voronoiNum; /* ボロノイボテンのナンバリング */
}mpoint_t;

typedef struct{
    int oriPoint;
    double psi;
	double theta;
	double phi;
	double A[3][3];
}dpoint_t;

typedef struct{
  position_t minBox;
  position_t maxBox;
  int used;
}box_t;

position_t *sample;/*最後の結果の構造体*/
box_t *smallBox;
position_t *atoms;
mpoint_t *mpoint;
dpoint_t *datapoint;
position_t *nearAtm;
position_t latticeConst;
position_t latticeLength;
position_t latorigin; // origin of simulation box with boundary
position_t alloribox; // origin of all simulation
position_t simboxori; // origin of simulation box
double massa;
double *adis;

char atom1[10];
char lattice[50];
int ParticleNum;
extern int NA;

/* 構造体、関数のプロトタイプ宣言 */
extern void config();
extern void settingParameter();
extern double GetRandom(double min, double max);
extern int ainOrOut(position_t atom, mpoint_t mpoint[], double *disv, int pt);
extern void checkrange(mpoint_t mpoint[]);
extern void relocate(mpoint_t mpoint[]);
extern void swaping_small_big(position_t Ddata[], double swapBy[], int Ndata, const char *fname, double adjust);
extern void writeFile(mpoint_t mpoint[], dpoint_t datapoint[]);
extern void SetupBox();
extern void RandomBox(mpoint_t mpoint[], dpoint_t datapoint[]);
extern void RandomCheck(mpoint_t mpoint[]);
extern void copyPoint(mpoint_t mpoint[], dpoint_t datapoint[]);
extern int far(position_t place, position_t compare);
extern void StartTimer();
extern double GetTimer();
extern void PrintTimer(double timer);
extern int inSide(position_t pt, position_t ori, position_t inR);
extern double sizeVoronoi(position_t ato[], int total, double map[]);
extern int readpointdata();
extern void eulervalue(dpoint_t datapoint[]);
extern void MakeRotateEle(position_t atoms[], mpoint_t mpoint[], dpoint_t datapoint[], int pt);
extern position_t arotate2d(position_t atoms, position_t cent, dpoint_t data);
extern position_t arotate3d(position_t atoms, position_t cent, dpoint_t data);
extern position_t rotate(position_t atoms, position_t cent, dpoint_t data);
extern void Vregister(position_t atoms[], mpoint_t mpoint[], int total, int pt);

extern void fail(const char *error);

#endif
