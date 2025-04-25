/* Note:Your choice is C IDE */
/*画一条NURBS曲线*/
#include "stdafx.h"
#include "stdio.h"
#include "math.h"
#include <time.h>
#include <math.h>
#include <iostream>
#include "engine.h"
#include <stdlib.h>
#include <string.h>

using namespace std;
double mx[30000], my[30000], mz[30000];
double rmx[3000], rmy[3000], rmz[3000];
double smx[300000], smy[300000], smz[300000];
int countnum = 0;
int rcountnum = 0;
int scountnum = 0;
void deboor(double coeff[][3], double knot[], double weight[], int degree, int l, int dense, double d, int e, double *val);
void calp(double coeff[][3], double knot[], double weight[], double inputu, double *val);
void CalNurbs(double coeff[][3], double knot[], double weight[], int degree, int Section, int dense);
void CalNurbs1(double coeff[][3], double knot[], double weight[], int degree, int Section, int dense);//画交线用
void CalNurbs2(double coeff[][3], double knot[], double weight[], int degree, int Section, int dense);//画曲面用
void uccpnt(double coeff[9][9][3], double weight[9][9], double knotw[], double w, int wmici, int k, double *val);
//void wccpnt(double u, int mici);
void udbor(double coeff[9][9][3], double weight[9][9], double knotw[], double d, int e, int n, double *val);
void conpoint(int n, double px[], double py[], double pz[], double *valx, double *valy, double *valz);
void reversep(int n, double xingzhip[10][3], double *valx, double *valy, double *valz);
//void wdbor(double d, int e, int n);

int main()
{
	Engine *ep;
	mxArray *TX = NULL, *TY = NULL, *TZ = NULL, *TA = NULL, *TB = NULL, *TC = NULL;
	mxArray *TTX = NULL, *TTY = NULL, *TTZ = NULL;
	mxArray *RTX = NULL, *RTY = NULL, *RTZ = NULL;
	mxArray *STX = NULL, *STY = NULL, *STZ = NULL;
	double tx[1000], ty[1000];// tz[1000];
	double coeff[9][9][3] = { { { 30,20,20 },{ 30,40,20 },{ 30,65,20 },{ 30,85,20 },{ 30,110,20 },{ 30,130,20 },{ 30,155,20 },{ 30,180,20 },{ 30,200,20 } },
	{ { 60,20,20 },{ 60,45,30 },{ 60,75,40 },{ 60,100,50 },{ 60,130,60 },{ 60,155,50 },{ 60,185,40 },{ 60,215,30 },{ 60,240,20 } },
	{ { 100,20,20 },{ 100,50,35 },{ 100,85,50 },{ 100,115,65 },{ 100,150,80 },{ 100,180,65 },{ 100,215,50 },{ 100,250,35 },{ 100,280,20 } },
	{ { 140,20,20 },{ 140,55,40 },{ 140,95,60 },{ 140,130,80 },{ 140,170,100 },{ 140,205,80 },{ 140,245,60 },{ 140,280,40 },{ 140,320,20 } },
	{ { 200,20,20 },{ 205,60,45 },{ 200,105,70 },{ 195,145,95 },{ 200,190,120 },{ 190,230,95 },{ 200,275,70 },{ 205,310,45 },{ 210,360,20 } },
	{ { 240,20,20 },{ 240,55,40 },{ 240,95,60 },{ 240,130,80 },{ 235,170,100 },{ 240,205,80 },{ 245,245,60 },{ 240,280,40 },{ 240,320,20 } },
	{ { 300,20,20 },{ 295,50,35 },{ 300,85,50 },{ 305,115,65 },{ 300,150,80 },{ 300,180,65 },{ 310,215,50 },{ 300,250,35 },{ 300,280,20 } },
	{ { 350,20,20 },{ 360,45,30 },{ 350,75,40 },{ 340,100,50 },{ 350,130,60 },{ 350,155,50 },{ 360,185,40 },{ 350,215,30 },{ 350,240,20 } },
	{ { 400,20,20 },{ 400,40,20 },{ 400,65,20 },{ 400,85,20 },{ 400,110,20 },{ 400,130,20 },{ 400,155,20 },{ 400,180,20 },{ 400,200,20 } } }; //控制点，四边直线*/
	double knotu[13] = { 0,0,0,0,0.17,0.33,0.50,0.67,0.83,1,1,1,1 };  //节点
	double knotw[13] = { 0,0,0,0,0.17,0.33,0.50,0.67,0.83,1,1,1,1 };
	double weight[9][9] = { 1,1,1,1,1,1,1,1,1,
		1,1,1,1,1,1,1,1,1,
		1,1,1,1,1,1,1,1,1,
		1,1,1,1,1,1,1,1,1,
		1,1,1,1,1,1,1,1,1,
		1,1,1,1,1,1,1,1,1,
		1,1,1,1,1,1,1,1,1,
		1,1,1,1,1,1,1,1,1,
		1,1,1,1,1,1,1,1,1 }; //权值
	double revweight[12] = { 1,1,1,1,1,1,1,1,1,1,1,1}; //权值

	double knotc[16] = { 0,0,0,0,0.04,0.17,0.33,0.50,0.67,0.83,0.95,0.98,1,1,1,1 };//插值节点后的节点向量
	//double coeff[17][2]={40,240,100,140,170,260,230,320,300,280,360,200,420,120,480,220,540,300,600,200,700,400,600,300,500,100,400,200,300,400,200,500,100,300}; /*控制点*/
	//double knot[21]= {0,0,0,0,0.04,0.08,0.17,0.25,0.33,0.41,0.5,0.58,0.67,0.75,0.83,0.89,0.96,1,1,1,1}; /*节点*/
	//double weight[17]= {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}; /*权值 */
	//double tu;
	//int ts;
	int num = 0;
	double a;
	double b;
	double c;
	double pdata[3];
	double ma[1], mb[1], mc[1];
	double w=0.0;
	int wmici=3;
	double conval[4];
//	double convalw;
	double coeffu[10][9][3];
	double weightu[10][9];
	double revconp[10][3];
	double rcpx[20], rcpy[20], rcpz[20];
	double revcoeff[9][3];
	double scoeffu[9][3];
	double sweightu[9];
	//while (w <=1)
	//{
		//if (w > knotw[wmici + 1])
			//wmici++;
		//cout << "w=" << wmici << endl;

	//画曲面，看实际效果
	while (w <= 1)
	{
		if (w > knotw[wmici + 1])
			wmici++;
		//cout << "w=" << scountnum << endl;

		for (int k = 0; k <= 8; k++)
		{
			uccpnt(coeff, weight, knotw, w, wmici, k, conval);//固定一个参数后求另一个方向曲线的控制点，以便生成曲线
			scoeffu[k][0] = conval[0];
			scoeffu[k][1] = conval[1];
			scoeffu[k][2] = conval[2];
			sweightu[k] = conval[3];
			//cout << scoeffu[k][0] << "   " << scoeffu[k][1] << "   " << scoeffu[k][2] << "   " << sweightu[k] << endl;
		}

		CalNurbs2(scoeffu, knotu, sweightu, 3, 6, 50);

		w = w + 0.001;
	}
	//求要和平面相交的曲线的控制点
	wmici = 3;
	for (int i = 3; i <= 12; i++)
	{
		w = knotc[i];
		if (w > knotw[wmici + 1])
			wmici++;
		//wmici = i;
//		cout << w << endl;
		//if (wmici > 8)
			//wmici = 8;
		for (int k = 0; k <=8; k++)
		{
			uccpnt(coeff, weight, knotw, w, wmici, k, conval);
			coeffu[i - 3][k][0] = conval[0];
			coeffu[i - 3][k][1] = conval[1];
			coeffu[i - 3][k][2] = conval[2];
			weightu[i - 3][k] = conval[3];
		}
	}
		//w = w + 0.05;
	//}
	
	//求和平面相交的曲线
	for (int i = 0; i <= 9; i++)
	{
		CalNurbs(coeffu[i], knotu, weightu[i], 3, 6, 100);
	}
	
	if (!(ep = engOpen("\0")))
	{
		fprintf(stderr, "\nCan't start MATLAB engine\n");
		//return EXIT_FAILURE;
	}
	TX = mxCreateDoubleMatrix(1, (countnum - 1), mxREAL);
	TY = mxCreateDoubleMatrix(1, (countnum - 1), mxREAL);
	TZ = mxCreateDoubleMatrix(1, (countnum - 1), mxREAL);
	TA = mxCreateDoubleMatrix(1, 1, mxREAL);
	TB= mxCreateDoubleMatrix(1, 1, mxREAL);
	TC = mxCreateDoubleMatrix(1, 1, mxREAL);

	TTX = mxCreateDoubleMatrix(1, 1000, mxREAL);
	TTY = mxCreateDoubleMatrix(1, 1000, mxREAL);
	TTZ = mxCreateDoubleMatrix(1, 1000, mxREAL);

	STX = mxCreateDoubleMatrix(1, (scountnum - 1), mxREAL);
	STY = mxCreateDoubleMatrix(1, (scountnum - 1), mxREAL);
	STZ = mxCreateDoubleMatrix(1, (scountnum - 1), mxREAL);

	memcpy((void *)mxGetPr(TX), (void *)mx, (countnum - 1) * sizeof(double));
	memcpy((void *)mxGetPr(TY), (void *)my, (countnum - 1) * sizeof(double));
	memcpy((void *)mxGetPr(TZ), (void *)mz, (countnum - 1) * sizeof(double));
	engPutVariable(ep, "TX", TX);
	engPutVariable(ep, "TY", TY);
	engPutVariable(ep, "TZ", TZ);
	memcpy((void *)mxGetPr(STX), (void *)smx, (scountnum - 1) * sizeof(double));
	memcpy((void *)mxGetPr(STY), (void *)smy, (scountnum - 1) * sizeof(double));
	memcpy((void *)mxGetPr(STZ), (void *)smz, (scountnum - 1) * sizeof(double));
	engPutVariable(ep, "STX", STX);
	engPutVariable(ep, "STY", STY);
	engPutVariable(ep, "STZ", STZ);
	for (int i = 0; i <=400; i++)
	{
		tx[i] = 50+i;
	}
	for (int j = 0; j <= 100; j++)
	{
		ty[j] = 100+j;
	}
	
	memcpy((void *)mxGetPr(TTY), (void *)ty, 101 * sizeof(double));
	memcpy((void *)mxGetPr(TTX), (void *)tx, 401 * sizeof(double));
	
	engPutVariable(ep, "TTY", TTY);
	engPutVariable(ep, "TTX", TTX);

	engEvalString(ep, "[TMX,TMY]=meshgrid(TTX,TTY);");
	engEvalString(ep, "TTZ=0.1*TMX-0.4*TMY+35;");//0.1*TMX-0.4*TMY+35
	//engEvalString(ep, "TTZ=0.5*TMX-0.5*TMY+20;");
	//engEvalString(ep, "plot(TTY,TTX);");
	//engEvalString(ep, "TTX=0.1*TTZ+0.1*TTY;");
	engEvalString(ep, "mesh(TTX,TTY,TTZ);");//平面
	
	//start=clock();
	
	engEvalString(ep, "hold on;");
	engEvalString(ep, "plot3(TX,TY,TZ,'b.');");
	engEvalString(ep, "plot3(STX,STY,STZ,'b.');");//曲面
	//engEvalString(ep, "hold on;");
	//engEvalString(ep, "ezmesh('y-5*x+700');");
	
	//end=clock();
	//printf("the time is%f ", (end-start)/CLK_TCK);

	//求曲线和平面的交点
	for (int i = 0; i <= 9; i++)
	{
		a = 0;
		b = 1;
		while ((b - a) / 2 > 0.000000001)
		{
			c = (a + b) / 2;
			calp(coeffu[i], knotu, weightu[i], c, pdata); 
			if (fabs(0.1* pdata[0] - 0.4 * pdata[1] - pdata[2] + 35) < 0.0001)
				break;
			else if (0.1* pdata[0] - 0.4 * pdata[1] - pdata[2] + 35< 0)
				b = c;
			else
				a = c;
		}

		ma[0] = pdata[0];
		mb[0] = pdata[1];
		mc[0] = pdata[2];
		revconp[i][0] = pdata[0];//求出的交点坐标
		revconp[i][1] = pdata[1];
		revconp[i][2] = pdata[2];
		memcpy((void *)mxGetPr(TA), (void *)ma, 1 * sizeof(double));
		memcpy((void *)mxGetPr(TB), (void *)mb, 1 * sizeof(double));
		memcpy((void *)mxGetPr(TC), (void *)mc, 1 * sizeof(double));
		engPutVariable(ep, "TA", TA);
		engPutVariable(ep, "TB", TB);
		engPutVariable(ep, "TC", TC);
		engEvalString(ep, "hold on;");
		engEvalString(ep, "plot3(TA,TB,TC,'y*');");//交点
		printf("c=%f  ", c);
		printf("x=%lf ", ma[0]);
		printf("y=%lf ", mb[0]);
		printf("z=%lf\n ", mc[0]);
		//printf("x=%d ", num);
	}
	//反求过交点的曲线控制点
	reversep(10, revconp, rcpx, rcpy, rcpz);
	for (int m = 0; m <= 11; m++)
	{
		revcoeff[m][0] = rcpx[m];
		revcoeff[m][1] = rcpy[m];
		revcoeff[m][2] = rcpz[m];
		printf("x=%lf   ", revcoeff[m][0]);
		printf("y=%lf   ", revcoeff[m][1]);
		printf("z=%lf\n ", revcoeff[m][2]);
	}
	RTX = mxCreateDoubleMatrix(1, 1000, mxREAL);
	RTY = mxCreateDoubleMatrix(1, 1000, mxREAL);
	RTZ = mxCreateDoubleMatrix(1, 1000, mxREAL);
	//画过交点的交线
	CalNurbs1(revcoeff, knotc, revweight, 3, 9, 100);
	
	memcpy((void *)mxGetPr(RTX), (void *)rmx, 1000 * sizeof(double));
	memcpy((void *)mxGetPr(RTY), (void *)rmy, 1000 * sizeof(double));
	memcpy((void *)mxGetPr(RTZ), (void *)rmz, 1000 * sizeof(double));
	engPutVariable(ep, "RTX", RTX);
	engPutVariable(ep, "RTY", RTY);
	engPutVariable(ep, "RTZ", RTZ);
	engEvalString(ep, "hold on;");
	engEvalString(ep, "plot3(RTX,RTY,RTZ,'r.');");//交线
	fgetc(stdin);
	mxDestroyArray(TX);
	mxDestroyArray(TY);
	mxDestroyArray(TZ);
	mxDestroyArray(TA);
	mxDestroyArray(TB);
	mxDestroyArray(TC);
	mxDestroyArray(TTX);
	mxDestroyArray(TTY);
	mxDestroyArray(TTZ);
	mxDestroyArray(RTX);
	mxDestroyArray(RTY);
	mxDestroyArray(RTZ);
	mxDestroyArray(STX);
	mxDestroyArray(STY);
	mxDestroyArray(STZ);
	engEvalString(ep, "close;");
	engClose(ep);

	return 0;
}
void deboor(double coeff[][3], double knot[], double weight[], int degree, int l, int dense, double d, int e, double *val)
{
	double pointx, pointy, pointz, unpoint, rpointx, rpointy, rpointz;
	double t1, t2;//coeffa[9][2],weighta[9];
	//double rep[2];
	//int k, i, j, m;
	int j;
	double **coeffa, *weighta;
	coeffa = new double*[degree + l];
	
	for (int n3 = 0; n3<degree + l; n3++)
		coeffa[n3] = new double[3];

	weighta = new double[degree + l];

	for (int m = e - degree; m <= e; m++)
	{
		weighta[m] = weight[m];
		for (int i = 0; i <= 2; i++)
		{
			coeffa[m][i] = coeff[m][i] * weight[m];
		}
	}
	for (int k = 1; k <= degree; k++)
	{
		for (j = e; j >= e - degree + k; j--)
		{
			t1 = (knot[j + degree - k + 1] - d) / (knot[j + degree - k + 1] - knot[j]);
			t2 = 1.0 - t1;
			coeffa[j][0] = t1*coeffa[j - 1][0] + t2*coeffa[j][0];
			coeffa[j][1] = t1*coeffa[j - 1][1] + t2*coeffa[j][1];
			coeffa[j][2] = t1*coeffa[j - 1][2] + t2*coeffa[j][2];
			weighta[j] = t1*weighta[j - 1] + t2*weighta[j];
		}
	}

	pointx = coeffa[j + 1][0];
	pointy = coeffa[j + 1][1];
	pointz = coeffa[j + 1][2];
	unpoint = weighta[j + 1];
	rpointx = pointx / unpoint;
	rpointy = pointy / unpoint;
	rpointz = pointz / unpoint;
	val[0] = rpointx;
	val[1] = rpointy;
	val[2] = rpointz;
}
//计算一个NURBS点
void calp(double coeff[][3], double knot[], double weight[], double inputu, double *val)
{
	int ts;
	for (int tk = 3; tk < 9; tk++)
	{
		if (knot[tk] <= inputu&&inputu <= knot[tk + 1])
		{
			ts = tk;
			break;
		}
	}

	deboor(coeff, knot, weight, 3, 6, 10, inputu, ts, val);
}

void CalNurbs(double coeff[][3], double knot[], double weight[], int degree, int Section, int dense)
{
	double data[3];
	double u;
	//int countnum = 0;
	for (int kk = degree; kk<Section + degree; kk++)
	{
		if (knot[kk + 1]>knot[kk])
		{
			for (int ii = 0; ii<dense; ii++)
			{
				u = knot[kk] + ii*(knot[kk + 1] - knot[kk]) / dense;
				deboor(coeff, knot, weight, degree, Section, dense, u, kk, data);
				mx[countnum] = data[0];
				my[countnum] = data[1];
				mz[countnum] = data[2];
				countnum++;
			}
		}
	}
}

void CalNurbs1(double coeff[][3], double knot[], double weight[], int degree, int Section, int dense)
{
	double data[3];
	double u;
	//int countnum = 0;
	for (int kk = degree; kk<Section + degree; kk++)
	{
		if (knot[kk + 1]>knot[kk])
		{
			for (int ii = 0; ii<dense; ii++)
			{
				u = knot[kk] + ii*(knot[kk + 1] - knot[kk]) / dense;
				deboor(coeff, knot, weight, degree, Section, dense, u, kk, data);
				rmx[rcountnum] = data[0];
				rmy[rcountnum] = data[1];
				rmz[rcountnum] = data[2];
				rcountnum++;
			}
		}
	}
}

void CalNurbs2(double coeff[][3], double knot[], double weight[], int degree, int Section, int dense)
{
	double data[3];
	double u;
	//int countnum = 0;
	for (int kk = degree; kk<Section + degree; kk++)
	{
		if (knot[kk + 1]>knot[kk])
		{
			for (int ii = 0; ii<dense; ii++)
			{
				u = knot[kk] + ii*(knot[kk + 1] - knot[kk]) / dense;
				deboor(coeff, knot, weight, degree, Section, dense, u, kk, data);
				smx[scountnum] = data[0];
				smy[scountnum] = data[1];
				smz[scountnum] = data[2];
				scountnum++;
			}
		}
	}
}


void uccpnt(double coeff[9][9][3], double weight[9][9], double knotw[],  double w, int wmici, int k, double *val)
{
	double pval[4];
	udbor(coeff, weight, knotw, w, wmici, k, pval);
	val[0] = pval[0];
	val[1] = pval[1];
	val[2] = pval[2];
	val[3] = pval[3];
}

void udbor(double coeff[9][9][3], double weight[9][9], double knotw[], double d, int e, int n, double *val)
{
	int k, i, j, m;
	double t1, t2;
	double coeffa[9][3];
	double weighta[9];

	for (m = e - 3; m <= e; m++)
	{
		weighta[m] = weight[m][n];
		for (i = 0; i <= 2; i++)
		{
			coeffa[m][i] = coeff[m][n][i] * weight[m][n];
			//printf("%lf\n", coeff[m][n][i]);
		}
	}
	for (k = 1; k <= 3; k++)
	{
		for (j = e; j >= e - 3 + k; j--)
		{
			t1 = (knotw[j + 3 - k + 1] - d) / (knotw[j + 3 - k + 1] - knotw[j]);
			t2 = 1.0 - t1;
			coeffa[j][0] = t1*coeffa[j - 1][0] + t2*coeffa[j][0];
			coeffa[j][1] = t1*coeffa[j - 1][1] + t2*coeffa[j][1];
			coeffa[j][2] = t1*coeffa[j - 1][2] + t2*coeffa[j][2];
			weighta[j] = t1*weighta[j - 1] + t2*weighta[j];
		}

	}
	val[0] = coeffa[j + 1][0];
	val[1] = coeffa[j + 1][1];
	val[2] = coeffa[j + 1][2];
	val[3] = weighta[j + 1];
}
//反求控制点，它调用conpoint函数，参数为矩阵的阶数，型值点和控制点返回值。
void reversep(int n, double xingzhip[10][3], double *valx, double *valy, double *valz)
{
	double px[20], py[20], pz[20];
//	double valx[20], valy[20], valz[20];
	px[0] = 6 * xingzhip[0][0];
	py[0] = 6 * xingzhip[0][1];
	pz[0] = 6 * xingzhip[0][2];
	px[n - 1] = 6 * xingzhip[n - 1][0];
	py[n - 1] = 6 * xingzhip[n - 1][1];
	pz[n - 1] = 6 * xingzhip[n - 1][2];
	for (int mm = 1; mm <= n - 2; mm++)
	{
		px[mm] = xingzhip[mm][0];
		py[mm] = xingzhip[mm][1];
		pz[mm] = xingzhip[mm][2];
	}
	conpoint(n, px, py, pz, valx,valy, valz);
	valx[0] = xingzhip[0][0];
	valy[0] = xingzhip[0][1];
	valz[0] = xingzhip[0][2];
	valx[n + 1] = xingzhip[n - 1][0];
	valy[n + 1] = xingzhip[n - 1][1];
	valz[n + 1] = xingzhip[n - 1][2];
}
//追赶法反求控制点，参数分别是型值点的三个坐标以及返回的控制点
void conpoint(int n, double px[], double py[], double pz[], double *valx, double *valy, double *valz)
{
	double c[9]={-3.0, 1.0/6.0,  1.0 / 6.0, 1.0 / 6.0,1.0 / 6.0,1.0/6.0, 1.0/6.0,1.0/6.0,1.0/4.0};
	double b[10]={9.0, 7.0/12.0,  2.0 / 3.0,2.0 / 3.0,2.0 / 3.0, 2.0/3.0,2.0/3.0,2.0/3.0,7.0/12.0, 9.0};
	double d[9]={1.0/4.0, 1.0/6.0,  1.0 / 6.0,1.0 / 6.0,1.0 / 6.0,1.0/6.0, 1.0/6.0,1.0/6.0,-3.0};
	double l[20];
	double s[20];
	double vx[20];
	double vy[20];
	double vz[20];
	double ax[20];
	double ay[20];
	double az[20];
	//int i;
	//int j;

	l[0] = b[0];
	vx[0] = px[0] / l[0];
	vy[0] = py[0] / l[0];
	vz[0] = pz[0] / l[0];
	for (int i = 0; i<n - 1; i++)
	{
		s[i] = c[i] / l[i];
		l[i + 1] = b[i + 1] - d[i] * s[i];
		vx[i + 1] = (px[i + 1] - d[i] * vx[i]) / l[i + 1];
		vy[i + 1] = (py[i + 1] - d[i] * vy[i]) / l[i + 1];
		vz[i + 1] = (pz[i + 1] - d[i] * vz[i]) / l[i + 1];
	}
	ax[n - 1] = vx[n - 1];
	ay[n - 1] = vy[n - 1];
	az[n - 1] = vz[n - 1];
	valx[n] = ax[n - 1];
	valy[n] = ay[n - 1];
	valz[n] = az[n - 1];

	for (int i = n - 2; i >= 0; i--)
	{
		ax[i] = vx[i] - s[i] * ax[i + 1];
		ay[i] = vy[i] - s[i] * ay[i + 1];
		az[i] = vz[i] - s[i] * az[i + 1];
		valx[i + 1] = ax[i];
		valy[i + 1] = ay[i];
		valz[i + 1] = az[i];
		//cout << valx[i + 1] << "  " << valy[i + 1] << "   " << valz[i + 1] << endl;
	}
	//cout << valx[n] << "  " << valy[n] << "   " << valz[n] << endl;
}
/*求每条w向曲线的控制点*/
/*void wccpnt(double u, int mici, int k)
{
for (int i = k - 3; i <= k; i++)
{
wdbor(u, mici, i);
wcoeff[i][0] = wcpnx;
wcoeff[i][1] = wcpny;
wcoeff[i][2] = wcpnz;
wwei[i] = wweight;
}
}*/

/*void wdbor(double d, int e, int n)
{
	int k, i, j, m;
	double t1, t2;
	double coeffa[9][3];
	double weighta[9];

	for (m = e - degree; m <= e; m++)
	{
		weighta[m] = weight[n][m];
		for (i = 0; i <= 2; i++)
		{
			coeffa[m][i] = coeff[n][m][i] * weight[n][m];
		}
	}
	for (k = 1; k <= degree; k++)
	{
		for (j = e; j >= e - degree + k; j--)
		{
			t1 = (knotu[j + degree - k + 1] - d) / (knotu[j + degree - k + 1] - knotu[j]);
			t2 = 1.0 - t1;
			coeffa[j][0] = t1*coeffa[j - 1][0] + t2*coeffa[j][0];
			coeffa[j][1] = t1*coeffa[j - 1][1] + t2*coeffa[j][1];
			coeffa[j][2] = t1*coeffa[j - 1][2] + t2*coeffa[j][2];
			weighta[j] = t1*weighta[j - 1] + t2*weighta[j];
		}

	}
	wcpnx = coeffa[j + 1][0];
	wcpny = coeffa[j + 1][1];
	wcpnz = coeffa[j + 1][2];
	wweight = weighta[j + 1];
}*/






