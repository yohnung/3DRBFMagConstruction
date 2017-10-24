// procedure used definition source file
#include <iomanip>
#include <mkl.h>
#include "MacroAndMostUsedLibrary.h"
using namespace std;
#include "OverallVariablesDeclarification.h"
#include "ClassClarification.h"
#include "ReadingData.h"
#include "ProcedureClarification.h"

void CrdTrf2FlyDirec(int num, Point* Position, double observedValue[][Dim], double* Alpha, 
	double* max_posit, double* min_posit)
{
/******	Coordinates Transform to Flying Direction make it x-direction,			******/
/******	given num, position to find flying_direction, then specify Alpha and 	******/
/******	transform position, observed value. finally find max and min position	******/
	double flying_direction[Dim] = { 0 };
	int Num_obserPosit = M;								// Number of observed position
	int Num_Crafts = GroupNumber;						// Number of Crafts
	double(* CraftLeap)[Dim] = new double[Num_Crafts*(Num_obserPosit-1)][Dim];
	double value[Dim];									// use in coordinates transform by representing old coordinates
	int i, g;
  /******	default value, collimate with x-direction which means no coordinates transform	******/
	//flying_direction[0] = 1;
	//if (Dim > 1)
	//	flying_direction[1] = 0;
	//if (Dim > 2)
	//	flying_direction[2] = 0;
  /******	calculate flying direction by average distance vector	******/
	for (g = 0; g < Num_Crafts; g++)
	{
		for (i = 0; i < Num_obserPosit - 1; i++)
		{
			
			CraftLeap[g*(Num_obserPosit-1) + i][0] = Position[g*(Num_obserPosit) + i + 1].getx() 
				- Position[g*(Num_obserPosit) + i].getx();
			flying_direction[0] += CraftLeap[g*(Num_obserPosit-1) + i][0];
			if (Dim > 1)
			{
				CraftLeap[g*(Num_obserPosit-1) + i][1] = Position[g*(Num_obserPosit) + i + 1].gety() 
					- Position[g*(Num_obserPosit) + i].gety();
				flying_direction[1] += CraftLeap[g*(Num_obserPosit-1) + i][1];
			}
			if (Dim > 2)
			{
				CraftLeap[g*(Num_obserPosit-1) + i][2] = Position[g*(Num_obserPosit) + i + 1].getz() 
					- Position[g*(Num_obserPosit) + i].getz();
				flying_direction[2] += CraftLeap[g*(Num_obserPosit-1) + i][2];
			}
		}
	}	
  /******	find transform degree	******/
	specifyDegree(flying_direction, Alpha);		
  /******	according given Alpha[0] = alpha, Alpha[1] = gamma to transform coordinates	******/
	if (Dim == 1)
		cout << "This is a one-dimensional problem, need not coordinates transform!" << endl;
	if (Dim > 1)
	{
		for (i = 0; i< num; i++)
		{
			value[0] = Position[i].getx();	value[1] = Position[i].gety();	// old coordinates
    /****** rotate oldx-oldy to new x'-y coordinates by z-axis ******/
			ZAxisCrdTrans(-Alpha[0], value);								// need couterclock
			Position[i].assignx(value[0]); Position[i].assigny(value[1]);	// assign new coordinates
			ZAxisCrdTrans(-Alpha[0], observedValue[i]);						// need couterclock
    /****** rotate x'-oldz to new x-z coordinates by y-axis	******/
			if (Dim > 2)
			{
				value[2] = Position[i].getz();
				YAxisCrdTrans(Alpha[1], value);
				Position[i].assignx(value[0]); Position[i].assignz(value[2]);
				YAxisCrdTrans(Alpha[1], observedValue[i]);
			}
		}
	}
  /******	find maximum and minimum position value	******/
	max_posit[0] = Position[0].getx();
	min_posit[0] = Position[0].getx();
	if (Dim > 1)
	{
		max_posit[1] = Position[0].gety();
		min_posit[1] = Position[0].gety();
	}
	if (Dim > 2)
	{
		max_posit[2] = Position[0].getz();
		min_posit[2] = Position[0].getz();
	}
	for (i = 1; i < num; i++)
	{
		if (Position[i].getx() < min_posit[0])
			min_posit[0] = Position[i].getx();
		if (Position[i].getx() > max_posit[0])
			max_posit[0] = Position[i].getx();
		if (Dim > 1)
		{
			if (Position[i].gety() < min_posit[1])
				min_posit[1] = Position[i].gety();
			if (Position[i].gety() > max_posit[1])
				max_posit[1] = Position[i].gety();
		}
		if (Dim > 2)
		{
			if (Position[i].getz() < min_posit[2])
				min_posit[2] = Position[i].getz();
			if (Position[i].getz() > max_posit[2])
				max_posit[2] = Position[i].getz();
		}
	}
  /******	recycle dynamic memory	******/
	delete[] CraftLeap;
}
void specifyDegree(double* direction, double* Alpha)
{
/******	according to direction, find alpha and gamma to write in Alpha[0] and Alpha[1]	******/
	if (Dim == 1)								// don't need transform
		Alpha[0] = 0;
	if (Dim == 2)
	{
		Alpha[0] = acos(direction[0] / sqrt(direction[0] * direction[0]
			+ direction[1] * direction[1]));
		Alpha[1] = 0;							// useless
	}
	if (Dim == 3)
	{
		Alpha[2] = sqrt(direction[0] * direction[0]		// use as a temp value
			+ direction[1] * direction[1]);				// sqrt(x*x + y*y)
		Alpha[0] = acos(direction[0] / Alpha[2]);
		Alpha[1] = acos(Alpha[2] / sqrt(Alpha[2] * Alpha[2] + direction[2] * direction[2]));
		Alpha[2] = 0;							// useless
	}
}
void ZAxisCrdTrans(double Alpha, double value[Dim])
{
/******	according given Alpha transform coordinates around z axis	******/
	double x, y;			// new coordinates
	double oldx, oldy;		// old coordinates
	double valuex, valuey;	// new value
/*if (Dim == 1)*/
	//cout << "This is a one-dimensional problem, need not coordinates transform!" << endl;
	oldx = value[0]; oldy = value[1];
/****** rotate oldx-oldy to new x-y coordinates ******/
	x = oldx * cos(Alpha) - oldy * sin(Alpha);
	y = oldx * sin(Alpha) + oldy * cos(Alpha);
	value[0] = x; value[1] = y;
}
void YAxisCrdTrans(double Alpha, double value[Dim])
{
/******	according given Alpha transform coordinates around y axis	******/
	double x, z;			// new coordinates
	double oldx, oldz;		// old coordinates
	double valuex, valuez;	// new value
	oldx = value[0]; oldz = value[2];
/****** rotate oldx-oldz to new x-z coordinates ******/
	x =  oldx * cos(Alpha) + oldz * sin(Alpha);
	z = -oldx * sin(Alpha) + oldz * cos(Alpha);
	value[0] = x; value[2] = z;
}
void XAxisCrdTrans(double Alpha, double value[Dim])
{
/******	according given Alpha transform coordinates around x axis	******/
	double y, z;			// new coordinates
	double oldy, oldz;		// old coordinates
	double valuey, valuez;	// new value
/****** rotate oldy-oldz to new y-z coordinates by x-axis	******/
	oldy = value[1]; oldz = value[2];
	y = oldy * cos(Alpha) - oldz * sin(Alpha);
	z = oldy * sin(Alpha) + oldz * cos(Alpha);
	value[1] = y; value[2] = z;
}
void constructNodesWeb(double* max_posit, double* min_posit, Point* Node, Point* DistParam, double* TypicalLength)
{
/******	construct the web of nodes and specify the typical length of RBF according to distance between adjacent nodes	******/
	double stretchParam = 3.;							// DistanceParameter = stretchParameter * meanDelta;
	Point max_Position, min_Position;
	max_Position.specify(max_posit); min_Position.specify(min_posit);
	double* meandist = new double[Dim]();				// mean distance between Nodes
	Point* meanDelta = new Point[Dim]();				// to instruct the web of nodes
	int i, j, k;
/******	specify meandist and meanDelta by 1.4*(max_posit-min_posit)/(Num-1)	**********************************/
	if (Nx == 1)
	{
		cout << "Caution!   There is only 1 Node" << endl;
		meandist[0] = 0;
	}
	else
	{
		meandist[0] = 1.4 * (max_posit[0] - min_posit[0]) / (Nx - 1);
		TypicalLength[0] = meandist[0];								// Typical length of Radial Basic Function
	}
	meanDelta[0].assignx(meandist[0]);
	DistParam[0].assignx(meandist[0]);
	if (Dim > 1)
	{
		if (Ny == 1)
		{
			cout << " This is a two dimensional problem, but there is only 1-dimension ";
			cout << " distributed nodes along x-direction" << endl;
			meandist[1] = 0;
		}
		else
		{
			meandist[1] = 1.4 * (max_posit[1] - min_posit[1]) / (Ny - 1);
			TypicalLength[1] = meandist[1];							// Typical length of Radial Basic Function
		}
		meanDelta[1].assigny(meandist[1]);
		DistParam[0].assigny(meandist[1]);
	}
	if (Dim > 2)
	{
		if (Nz == 1)
		{
			cout << "This is a three dimensional problem, but there is only 2-dimension";
			cout << "distributed nodes in x-y plane" << endl;
			meandist[2] = 0;
		}
		else
		{
			meandist[2] = 1.4 * (max_posit[2] - min_posit[2]) / (Nz - 1);
			TypicalLength[2] = meandist[2];							// Typical length of Radial Basic Function
		}
		meanDelta[2].assignz(meandist[2]);
		DistParam[0].assignz(meandist[2]);
	}
	Node[0] = min_Position - 0.2 * (max_Position - min_Position);
/******* Specify DistParam, usually 3 times Node's distance which means 3 * 1.4 in dimensionless case	******/
	double* temptemp = new double[2]();				// DistParam[0] = stretchParam * DistParam[0];
	temptemp[0] = 4.2;								// DistParam[0] = stretchParam * DistParam[0];
	DistParam[0].specify(temptemp);					// DistParam[0] = stretchParam * DistParam[0];
	delete[] temptemp;
/******* construct the web of Nodes		**********************************************************************/
	for (i = 1; i < Nx; i++)
	{
		if (Dim > 1)
		{
			for (j = 1; j < Ny; j++)
			{
				if (Dim > 2)
					for (k = 1; k < Nz; k++)
					{																
						Node[k*Nx*Ny + (j - 1)*Nx + (i - 1)] = Node[(k - 1)*Nx*Ny
							+ (j - 1)*Nx + (i - 1)] + meanDelta[2];
						DistParam[k*Nx*Ny + (j - 1)*Nx + (i - 1)] = DistParam[0];	// every time z plus a meanDelta[2]
					}
				Node[j*Nx + (i - 1)] = Node[(j - 1)*Nx + (i - 1)] + meanDelta[1];	// every time y plus a meanDelta[1]
				DistParam[j*Nx + (i - 1)] = DistParam[0];
			}
			if (Dim > 2)
				for (k = 1; k < Nz; k++)
				{
					Node[k*Nx*Ny + (Ny - 1)*Nx + (i - 1)] = Node[(k - 1)*Nx*Ny
						+ (Ny - 1)*Nx + (i - 1)] + meanDelta[2];
					DistParam[k*Nx*Ny + (Ny - 1)*Nx + (i - 1)] = DistParam[0];
				}
		}
		Node[i] = Node[i - 1] + meanDelta[0];										// every time x plus a meanDelta[0]
		DistParam[i] = DistParam[0];
	}
	if (Dim > 1)
	{
		for (j = 1; j < Ny; j++)
		{
			if (Dim > 2)
				for (k = 1; k < Nz; k++)
				{
					Node[k*Nx*Ny + (j - 1)*Nx + Nx - 1] =
						Node[(k-1)*Nx*Ny + (j - 1)*Nx + Nx - 1] + meanDelta[2];
					DistParam[k*Nx*Ny + (j - 1)*Nx + Nx - 1] = DistParam[0];
				}
			Node[j*Nx + Nx - 1] = Node[(j - 1)*Nx + Nx - 1] + meanDelta[1];
			DistParam[j*Nx + Nx - 1] = DistParam[0];
		}
		if (Dim > 2)
			for (k = 1; k < Nz; k++)
			{
				Node[k*Nx*Ny + (Ny - 1)*Nx + Nx - 1] =
					Node[(k - 1)*Nx*Ny + (Ny - 1)*Nx + Nx - 1] + meanDelta[2];
				DistParam[k*Nx*Ny + (Ny - 1)*Nx + Nx - 1] = DistParam[0];
			}
	}
	delete[] meandist;
	delete[] meanDelta;
	ofstream dimout("DimensionInformation.dat");	// write out dimension infromation
	dimout << Dim << endl << Nx << endl << Ny << endl << Nz << endl;
	dimout.close();
}
void CrdTrfFromFlyDirec(int num, double* Alpha, Point* Position, double observedValue[][Dim], 
	Point* Node, Point* DistParam, double* TypicalLenght)
{
/****** coordinates transform form x-flying direction to original one according to degree Alpha ******/
	double value[Dim];
	int i;
	if (Dim == 1)
		cout << "This is a one-dimensional problem, need not coordinates transform!" << endl;
	if (Dim > 1)
	{
		for (i = 0; i< num; i++)		// transform position and observed value
		{
			value[0] = Position[i].getx(); value[1] = Position[i].gety();		// old coordinates
  /****** rotate new y-z back to y'-oldz coordinates by x-axis	******/
			if (Dim > 2)
			{
				value[2] = Position[i].getz();
				YAxisCrdTrans(-Alpha[1], value);								// rotate back
				Position[i].assignx(value[0]); Position[i].assignz(value[2]);
				YAxisCrdTrans(-Alpha[1], observedValue[i]);
			}
  /****** rotate new x-y' to oldx-oldy coordinates by z-axis ******/
			ZAxisCrdTrans(Alpha[0], value);						// rotate back, minus compared with 'CrdTrf2FlyDirec'
			Position[i].assignx(value[0]); Position[i].assigny(value[1]);	// assign new coordinates
			ZAxisCrdTrans(Alpha[0], observedValue[i]);
		}
		for (i = 0; i< N; i++)			// transform Node
		{

			value[0] = Node[i].getx(); value[1] = Node[i].gety();				// old coordinates
  /****** rotate new y-z to y'-oldz coordinates by x-axis	******/
			if (Dim > 2)
			{
				value[2] = Node[i].getz();
				YAxisCrdTrans(-Alpha[1], value);
				Node[i].assignx(value[0]); Node[i].assignz(value[2]);
			}
  /****** rotate new x-y' to oldx-oldy coordinates by z-axis ******/
			ZAxisCrdTrans(Alpha[0], value);
			Node[i].assignx(value[0]); Node[i].assigny(value[1]);	// assign new coordinates  
		}
		for (i = 0; i< N; i++)			// transform DistParam
		{
			value[0] = DistParam[i].getx(); value[1] = DistParam[i].gety();		// old coordinates
  /****** rotate y'-oldz to new y-z coordinates by x-axis	******/
			if (Dim > 2)
			{
				value[2] = DistParam[i].getz();
				YAxisCrdTrans(-Alpha[1], value);
				DistParam[i].assignx(value[0]); DistParam[i].assignz(value[2]);
			}
  /****** rotate oldx-oldy to new x-y' coordinates by z-axis ******/
			ZAxisCrdTrans(Alpha[0], value);
			DistParam[i].assignx(value[0]); DistParam[i].assigny(value[1]);	// assign new coordinates
		}
	}
  /******	redefine typicla length used in Radial Basic Function	******/
	if (Dim > 2)
		YAxisCrdTrans(-Alpha[1], TypicalLenght);
	ZAxisCrdTrans(Alpha[0], TypicalLenght);
}
void CrdTrf2MinVarDir(int num, Point* Position, double observedValue[][Dim], double* Alpha, 
	Point* Node, Point* DistParam, double* TypicalLength, double* abs_max_value, double* abs_min_value)
{
/******	According to position and observe value to find a min varing direction,			*******/
/******	then specify degree Alpha, to make it z-directionand. Thus transform 'position',*******/
/******	'observed value', 'node', 'distparam', and find absolute maximum observed value.*******/
/******	Finally transform 'typical length' to new coordinates									*******/
	double min_vary_direc[Dim];		// minimum varing direction
	double value[Dim];				// use in coordinates transform by representing old coordinates
	int i, d;
	MVAfunction(Position, observedValue, min_vary_direc);		// find minimum varying direction
  /******	find transform degree	******/
	specifyDegree(min_vary_direc, Alpha);
  /******	according given Alpha[0] = alpha, Alpha[1] = gamma to transform coordinates	******/
	if (Dim == 1)
		cout << "This is a one-dimensional problem, need not coordinates transform!" << endl;
	if (Dim > 1)
	{
		for (i = 0; i< num; i++)		// transform position and observed value
		{
    /****** rotate oldx-oldy to new x-y' coordinates by z-axis ******/
			value[0] = Position[i].getx();	value[1] = Position[i].gety();	// old coordinates
			ZAxisCrdTrans(Pi/2.-Alpha[0], value);
			Position[i].assignx(value[0]); Position[i].assigny(value[1]);	// assign new coordinates
			ZAxisCrdTrans(Pi/2.-Alpha[0], observedValue[i]);
    /****** rotate y'-oldz to new y-z coordinates by x-axis	******/
			if (Dim > 2)
			{
				value[2] = Position[i].getz();
				XAxisCrdTrans(-Alpha[1], value);
				Position[i].assigny(value[1]); Position[i].assignz(value[2]);
				XAxisCrdTrans(-Alpha[1], observedValue[i]);
			}
		}
		for (i = 0; i< N; i++)			// transform Node
		{
    /****** rotate oldx-oldy to new x-y' coordinates by z-axis ******/
			value[0] = Node[i].getx();	value[1] = Node[i].gety();	// old coordinates
			ZAxisCrdTrans(Pi / 2. - Alpha[0], value);
			Node[i].assignx(value[0]); Node[i].assigny(value[1]);	// assign new coordinates
    /****** rotate y'-oldz to new y-z coordinates by x-axis	******/
			if (Dim > 2)
			{
				value[2] = Node[i].getz();
				XAxisCrdTrans(-Alpha[1], value);
				Node[i].assigny(value[1]); Node[i].assignz(value[2]);
			}
		}
		for (i = 0; i< N; i++)			// transform DistParam
		{
    /****** rotate oldx-oldy to new x-y' coordinates by z-axis ******/
			value[0] = DistParam[i].getx();	value[1] = DistParam[i].gety();	// old coordinates
			ZAxisCrdTrans(Pi / 2. - Alpha[0], value);
			DistParam[i].assignx(value[0]); DistParam[i].assigny(value[1]);	// assign new coordinates
    /****** rotate y'-oldz to new y-z coordinates by x-axis	******/
			if (Dim > 2)
			{
				value[2] = DistParam[i].getz();
				XAxisCrdTrans(-Alpha[1], value);
				DistParam[i].assigny(value[1]); DistParam[i].assignz(value[2]);
			}
		}
	}
  /******	find the maximun or minum observed value	***************************************/
	for (d = 0; d < Dim; d++)
	{
		abs_max_value[d] = abs(observedValue[0][d]);
		abs_min_value[d] = abs(observedValue[0][d]);
	}
	for (i = 1; i < num; i++)
	{
		for (d = 0; d < Dim; d++)
		{
			if (abs(observedValue[i][d]) < abs_min_value[d])
				abs_min_value[d] = abs(observedValue[i][d]);
			if (abs(observedValue[i][d]) > abs_max_value[d])
				abs_max_value[d] = abs(observedValue[i][d]);
		}
	}
  /******	redefine typicla length used in Radial Basic Function	******/
	ZAxisCrdTrans(Pi / 2. - Alpha[0], TypicalLength);
	if (Dim > 2)
		XAxisCrdTrans(-Alpha[1], TypicalLength);
}
void MVAfunction(Point* Position, double observedValue[][Dim], double* min_vary_direc)
{
/******	find minimum varying direction ******/
	min_vary_direc[0] = 0;
	if (Dim > 1)
	{
		min_vary_direc[1] = 1;
	}
	if (Dim > 2)
	{
		min_vary_direc[2] = 0;
	}
}
void specifyAY(double A[][N_Alpha], double* Y,	double observedValue[][Dim],
	Point* Position, int Number_obserPosit, BasicFunction* Chi)
{
	double* FirstDeriv = new double[Dim]();
	double(*SecondDeriv)[Dim] = new double[Dim][Dim]();
	double* ScalarDeriv = new double();
	double x, y, z, r;
	int i, j, k;
	if (Dim == 1)
	{
		for (i = 0; i < Number_obserPosit; i++)
		{
			for (j = 0; j < N; j++)
				A[i][j] = Chi[j].getValue(Position[i]);				// scalar B
			Y[i] = observedValue[i][0];
		}
	}
	if (Dim == 2)
	{
		for (i = 0; i < Number_obserPosit; i++)
		{
			for (j = 0; j < N; j++)
				A[2 * i][j] = Chi[j].get_yDeriv(Position[i]);		// Bx
			for (j = 0; j < N; j++)
				A[2 * i + 1][j] = -Chi[j].get_xDeriv(Position[i]);	// By
			Y[2 * i] = observedValue[i][0];							// Bx
			Y[2 * i + 1] = observedValue[i][1];						// By
		}
	}
	if (Dim == 3)
	{
		for (i = 0; i < Number_obserPosit; i++)
		{
			x = Position[i].getx(); y = Position[i].gety(); z = Position[i].getz();
			r = Position[i].getradialLength();
			for (j = 0; j < N; j++)
			{
				Chi[j].getDerivative(Position[i], FirstDeriv, SecondDeriv, ScalarDeriv);
				A[3 * i][2 * j] = 1 / r * (z*FirstDeriv[1] - y*FirstDeriv[2]);	// 1/r * (z * Partial_y - y * Partial_z)
				A[3 * i][2*j+1] = 2 * FirstDeriv[0] 
					+ (x*SecondDeriv[0][0] + y*SecondDeriv[0][1] + z*SecondDeriv[0][2]) 
					- x*(*ScalarDeriv);
				A[3*i+1][2 * j] = 1 / r*(x*FirstDeriv[2] - z*FirstDeriv[0]);	// 1/r * (x * Partial_z - z * Partial_x)
				A[3*i+1][2*j+1] = 2 * FirstDeriv[1] 
					+ (x*SecondDeriv[0][1] + y*SecondDeriv[1][1] + z*SecondDeriv[1][2]) 
					- y*(*ScalarDeriv);
				A[3*i+2][2 * j] = 1 / r*(y*FirstDeriv[0] - x*FirstDeriv[1]);	// 1/r * (y * Partial_x - x * Partial_y)
				A[3*i+2][2*j+1] = 2 * FirstDeriv[2] 
					+ (x*SecondDeriv[0][2] + y*SecondDeriv[1][2] + z*SecondDeriv[2][2]) 
					- z*(*ScalarDeriv);
			}
			Y[3 * i] = observedValue[i][0];							// Bx
			Y[3 * i + 1] = observedValue[i][1];						// By
			Y[3 * i + 2] = observedValue[i][2];						// Bz
		}
	}
	delete ScalarDeriv;
	delete[] FirstDeriv;
	delete[] SecondDeriv;
}
void generalSolver(double* Alpha, Point* Position, double observedValue[][Dim],
	int Number_obserPosit, BasicFunction* Chi)
{
/******	generally solve alpha	******************************************************/
	oncefit fitting;
	double(*fittedValue)[Dim] = new double[max_length][Dim]();
	double *tempAlpha = new double[N]();
	double SquareSum = 0;
	double tempSquareSum = 0;
	fitting.assignParameter(tempAlpha);
	int i, j, k;
/******	construct square sum of difference between observed and fitted value	******/
	for (i = 0; i < Number_obserPosit; i++)
	{
		fitting.getValue(fittedValue[i], Position[i], Chi);
		SquareSum += (fittedValue[i][0] - observedValue[i][0])*(fittedValue[i][0] - observedValue[i][0]);
		if (Dim > 1)
			SquareSum += (fittedValue[i][1] - observedValue[i][1])*(fittedValue[i][1] - observedValue[i][1]);
		if (Dim > 2)
			SquareSum += (fittedValue[i][2] - observedValue[i][2])*(fittedValue[i][2] - observedValue[i][2]);
	}
/****** using tempAlpha[N] and tempSquareSum to optimize fitting::Alpha[N]...
		... to be continued! to be continued! to be continued! to be continued!	******/
	for (i = 0; i < N; i++)
		Alpha[i] = tempAlpha[i];
	delete[] tempAlpha;
	delete[] fittedValue;
}
double LinearLUSolver(double* Alpha, double A[][N_Alpha], double* Y, int Number_obserPosit)
{
/******	redefining A \equiv A^T * A, and B \equiv A^T * y to LU solve A * alpha = B	******/
	double *tempAlpha = new double[N_Alpha]();
	double(*reA)[N_Alpha] = new double[N_Alpha][N_Alpha]();			// defined as A^T * A
	double* re_a = new double[N_Alpha*N_Alpha]();					// array a of reA
	double* B = new double[N_Alpha]();
	int* ipiv = new int[N_Alpha]();
	double anorm_infty = 0;											// reA's infinity-norm, which is used to compute condition number of reA
	double temp_anorm;												// temp reA's norm
	double acond;													// reA's condition number
	int i, j, k, d;
	for (i = 0; i < N_Alpha; i++)
	{
		temp_anorm = 0;
		for (j = 0; j < N_Alpha; j++)
		{
			for (k = 0; k < Number_obserPosit; k++)
			{
				for (d = 0; d < Dim; d++)
					reA[i][j] += A[Dim*k + d][i] * A[Dim*k + d][j];	// A^T * A
			}
			re_a[i * N_Alpha + j] = reA[i][j];
			temp_anorm += abs(reA[i][j]);							// temp norm of every row in reA
		}
		if (anorm_infty < temp_anorm)
			anorm_infty = temp_anorm;								// maximum array-norm of reA, which is infinity-anorm of reA
/******	define Matrix B	with B = A^T * Y	*********************************************/
		for (k = 0; k < Number_obserPosit; k++)
		{
			for (d = 0; d < Dim; d++)
				B[i] += A[Dim*k + d][i] * Y[Dim*k + d];
		}
		tempAlpha[i] = B[i];
	}
/****** solving the problem by resorting to LAPACK(Linear Algebric Package) library	******/
	LAPACKE_dgetrf(LAPACK_ROW_MAJOR, N_Alpha, N_Alpha, re_a, N_Alpha, ipiv);					// This and the following sentence usese LAPACK linear algebra, specificly LU factorization
	LAPACKE_dgecon(LAPACK_ROW_MAJOR, 'I', N_Alpha, re_a, N_Alpha, anorm_infty, &acond);			// estimate the condition number of A
	LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', N_Alpha, 1, re_a, N_Alpha, ipiv, tempAlpha, 1);		// Using LU factored Matrix to solve A * alpha = B, overwrite B
	ofstream MatrixOut("Condition Number and SqaureMatrix.dat");
	write(&acond, reA, MatrixOut);
	MatrixOut.close();
	ofstream LUout("LU_Alpha.dat");
	write(tempAlpha, N_Alpha, LUout);
	LUout.close();
/****** pass temp-Alpha to real Alpha	**************************************************/
	for (i = 0; i < N_Alpha; i++)
		Alpha[i] = tempAlpha[i];
/****** recycle dynamic memory and stop program	******************************************/
	delete[] ipiv;
	delete[] B;
	delete[] re_a;
	delete[] reA;
	delete[] tempAlpha;
	return acond;
}
void LinearQRSolver(double* Alpha, double A[][N_Alpha], double* Y, int Number_obserPosit, double* abs_max_value)
{
/******	directly solve A * alpha = y using QR Factorization	******/
	double *tempAlpha = new double[N_Alpha]();
	double* a = new double[Dim * Number_obserPosit * N_Alpha]();			// restore A matrix
	double* c = new double[Dim * Number_obserPosit]();						// restore observed value
	double* tau = new double[N_Alpha]();									// temp matrix
	int i, j, k, d;															// define a, c matrix or array using advanced approach
	for (i = 0; i < Number_obserPosit; i++)
	{
		for (d = 0; d < Dim; d++)
		{
			if (abs(Y[Dim * i + d]) > tolerance * abs_max_value[d])			// d=0, for B or Bx; d=1 for By; d=2 for Bz
			{
				for (j = 0; j < N_Alpha; j++)
					a[(Dim * i + d) * N_Alpha + j] = A[Dim * i + d][j];// / abs(Y[Dim*i+d]);
				c[Dim * i + d] = Y[Dim * i + d];// / abs(Y[Dim*i+d]);
			}
			else
			{
				for (j = 0; j < N_Alpha; j++)
					a[(Dim * i + d) * N_Alpha + j] = A[Dim * i + d][j];// / (tolerance * abs_max_value[d]);
				c[Dim * i + d] = Y[Dim * i + d];// / (tolerance * abs_max_value[d]);
			}
		}
	}
/****** solving the problem by resorting to LAPACK(Linear Algebric Package) library	******/
	LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, Dim * Number_obserPosit, N_Alpha, a, N_Alpha, tau);							// QR factorization, which gives that A = Q * R which Q is 2*Number_obserPosit-by-N matrix and R is N-by-N matrix
	LAPACKE_dormqr(LAPACK_ROW_MAJOR, 'L', 'T', Dim * Number_obserPosit, 1, N_Alpha, a, N_Alpha, tau, c, 1);			// compute C = Q^T * C, Q is restored in a
	cblas_dtrsm(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, N_Alpha, 1, 1, a, N_Alpha, c, 1);	// Comute R * x = C and C is overwittedn by the solution matrix x
	for (i = 0; i < N_Alpha; i++)
		tempAlpha[i] = c[i];
	ofstream QRout("QR_Alpha.dat");
	write(tempAlpha, N_Alpha, QRout);
	QRout.close();
/****** pass temp-Alpha to real Alpha	**************************************************/
	for (i = 0; i < N_Alpha; i++)
		Alpha[i] = tempAlpha[i];
/****** recycle dynamic memory and stop program	******************************************/
	delete[] tau;
	delete[] c;
	delete[] a;
	delete[] tempAlpha;
}
void LinearSVDSolver(double* Alpha, double A[][N_Alpha], double* Y, int Number_obserPosit, double* abs_max_value, double acond)
{
	double *tempAlpha = new double[N_Alpha]();
	double* a = new double[Dim * Number_obserPosit * N_Alpha]();		// restore A matrix
	double* c = new double[Dim * Number_obserPosit]();					// restore observed value
	double* singularValue = new double[N_Alpha]();
	int* rank = new int();
	double rcond;
	if (acond > 0)
		rcond = cutoff * sqrt(acond);								// cond(A^T*A) = cond(A)^2 
	else
		rcond = -1;
	int i, j, k, d;														// define a, c matrix or array again! Remember again!
	for (i = 0; i < Number_obserPosit; i++)
	{
		for (d = 0; d < Dim; d++)
		{
			if (abs(Y[Dim * i + d]) > tolerance * abs_max_value[d])		// d=0, for B or Bx; d=1 for By; d=2 for Bz
			{
				for (j = 0; j < N_Alpha; j++)
					a[(Dim * i + d) * N_Alpha + j] = A[Dim * i + d][j];// / abs(Y[Dim*i+d]);
				c[Dim * i + d] = Y[Dim * i + d];// / abs(Y[Dim*i+d]);
			}
			else
			{
				for (j = 0; j < N_Alpha; j++)
					a[(Dim * i + d) * N_Alpha + j] = A[Dim * i + d][j];// / (tolerance * abs_max_value[d]);
				c[Dim * i + d] = Y[Dim * i + d];// / (tolerance * abs_max_value[d]);
			}
		}
	}
	LAPACKE_dgelss(LAPACK_ROW_MAJOR, Dim * Number_obserPosit, N_Alpha, 1, a, N_Alpha, c, 1, singularValue, rcond, rank);
	for (i = 0; i < N_Alpha; i++)
		tempAlpha[i] = c[i];
	ofstream SVDout("SVD_Rank_Alpha_SingularValue.dat");
	SVDout << *rank << endl;
	write(tempAlpha, N_Alpha, SVDout);
	write(singularValue, N_Alpha, SVDout);
	SVDout.close();
/****** pass temp-Alpha to real Alpha	**************************************************/
	for (i = 0; i < N_Alpha; i++)
		Alpha[i] = tempAlpha[i];
/****** recycle dynamic memory and stop program	******************************************/	
	delete rank;
	delete[] singularValue;
	delete[] c;
	delete[] a;
	delete[] tempAlpha;
}
void write(Point* posit, int num, ofstream& filename)
{
/******	write out observation positions	******/
	int i;
	filename << setiosflags(ios::scientific) << setprecision(5);
	for (i = 0; i < num; i++)
		filename << setw(14) << posit[i].getx() << "  ";
	filename << endl;
	if (Dim > 1)
	{
		for (i = 0; i < num; i++)
			filename << setw(14) << posit[i].gety() << "  ";
		filename << endl;
	}
	if (Dim > 2)
	{
		for (i = 0; i < num; i++)
			filename << setw(14) << posit[i].getz() << "  ";
		filename << endl;
	}
}
void write(double value[][Dim], int num, ofstream& filename)
{
/******	write out observed magnetic field value	******/
	int i;
	filename << setiosflags(ios::scientific) << setprecision(5);
	for (i = 0; i < num; i++)
		filename << setw(14) << value[i][0] << "  ";
	filename << endl;
	if (Dim > 1)
	{
		for (i = 0; i < num; i++)
			filename  << setw(14) << value[i][1] << "  ";
		filename << endl;
	}
	if (Dim > 2)
	{
		for (i = 0; i < num; i++)
			filename << setw(14) << value[i][2] << "  ";
		filename << endl;
	}
	filename << endl;
}
void write(double* alpha_singularvalue, int num, ofstream& filename)
{
/******	write out optimized 'alpha's or singular value of Matrix	******/
	int n = N_Alpha / (N);		// if N is not enclosed by parenthesis, N_Alpha/N = Nx*Ny*Nz*2 / Nx*Ny*Nz is not 2 
	int i;
	filename << setiosflags(ios::scientific) << setprecision(15);
	for (i = 0; i < num / n; i++)
		filename << setw(25) << alpha_singularvalue[i*n] << "  ";
	filename << endl;
	if (n > 1)
	{
		for (i = 0; i < num / n; i++)
			filename << setw(25) << alpha_singularvalue[i * n + 1] << " ";
		filename << endl;
	}
}
void write(double* rcond, double Matrix[][N_Alpha], ofstream& filename)
{
/******	write out Matrix and its condition number	******/
	int i, j;
	filename << setiosflags(ios::scientific) << setprecision(15);
	filename << *rcond << endl;
	for (i = 0; i < N_Alpha; i++)
	{
		for (j = 0; j < N_Alpha; j++)
		{
			filename << setw(25) << Matrix[i][j] << " ";
		}
		filename << endl;
	}

}