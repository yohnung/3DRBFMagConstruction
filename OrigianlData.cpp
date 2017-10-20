#include<random>
#include<time.h>
#include "MacroAndMostUsedLibrary.h"
using namespace std;
#include "ClassClarification.h"
#include "ProcedureClarification.h"
#include "ReadingData.h"


using namespace std;
int read(Point* Position, double observedValue[][Dim], double* max_posit,
	double* min_posit, double* abs_max_value, double* abs_min_value)
{		
/******	to do a coordinates transformation and specify the craft relative position	*******************************/
	int Num_obserPosit;												// Number of observed position
	int Num_Crafts;													// Number of Crafts
	double* starting_posit = new double[Dim]();						// starting position (x,y,z)
	double* flying_posit = new double[Dim]();						// flying position
	double(*CraftsInterval)[Dim] = new double[10][Dim]();			// distance vector between maximum 10 different space-craft
	double(*CraftLeap)[Dim] = new double[max_length][Dim]();		// space-craft's interval between two adjacent time
	ifstream datain("OriginalObsevedData.dat");
/******	do a coordinates transformation and specify magnetic field value in new coordinates	***********************/
	CoordinatesTransform(observedValue, CraftsInterval, &Num_Crafts,
		CraftLeap, &Num_obserPosit, datain);
	datain.close();
	int i, j, k, d, g;
	for (i = 0; i < Num_obserPosit; i++)
	{
		CraftLeap[i][0] = Leapx;
		if (Dim > 1)
			CraftLeap[i][1] = Leapy;
		if (Dim > 2)
			CraftLeap[i][2] = Leapz;
	}
	starting_posit[0] = startx;
	flying_posit[0] = startx;
	max_posit[0] = startx; min_posit[0] = startx;
	CraftsInterval[0][0] = Interval12x;								// Interval[0] means distance vector between 1st and 2nd
	CraftsInterval[1][0] = Interval13x;								// Interval[1] for 1st and 3rd
	CraftsInterval[2][0] = Interval14x;								// Interval[2] for 1st and 4th
	if (Dim > 1)
	{
		starting_posit[1] = starty;
		flying_posit[1] = starty;
		max_posit[1] = starty; min_posit[1] = starty;
		CraftsInterval[0][1] = Interval12y;	
		CraftsInterval[1][1] = Interval13y;
		CraftsInterval[2][1] = Interval14y;
	}
	if (Dim > 2)
	{
		starting_posit[2] = startz;
		flying_posit[2] = startz;
		max_posit[2] = startz; min_posit[2] = startz;
		CraftsInterval[0][2] = Interval12z;	
		CraftsInterval[1][2] = Interval13z;
		CraftsInterval[2][2] = Interval14z;
	}
	Position[0].specify(flying_posit);
	default_random_engine generator(time(NULL));
	uniform_real_distribution<double> distribution(0.1, 0.3);		// mean value is 0.25; can be modified
/******	flying a space-craft, generating positionas		***********************************/
	for (g = 0; g < Num_Crafts; g++)								// g for different craft
	{
		if (g > 0)													// set for 2nd and 3rd craft
		{
			for (d = 0; d < Dim; d++)
			{
				flying_posit[d] = starting_posit[d] + CraftsInterval[g - 1][d];
				if (flying_posit[d] < min_posit[d])
					min_posit[d] = flying_posit[d];
				if (flying_posit[d] > max_posit[d])
					max_posit[d] = flying_posit[d];
			}
			Position[g*Num_obserPosit].specify(flying_posit);
		}
		for (i = 1; i < Num_obserPosit; i++)						// for a specific craft
		{
			for (d = 0; d < Dim; d++)
			{
				flying_posit[d] += CraftLeap[i][d];	// + distribution(generator);		
				if (flying_posit[d] < min_posit[d])
					min_posit[d] = flying_posit[d];
				if (flying_posit[d] > max_posit[d])
					max_posit[d] = flying_posit[d];
			}
			Position[g*Num_obserPosit + i].specify(flying_posit);
		}		
	}
/******	model the observed value, should be done by `CoordinatesTransform` subprogram	***/
	if (Dim == 1)
		for (i = 0; i < Num_Crafts*Num_obserPosit; i++)
			observedValue[i][0] = tanhtofit(Position[i]);			// f(x,y,z) to be fitted, specified manually
	if (Dim > 1)
	{
		Point* Node = new Point[N]();
		Point* DistParam = new Point[N]();
		constructNodesWeb(Node, DistParam, max_posit, min_posit);	// instruct web of node
		for (i = 0; i < Num_Crafts*Num_obserPosit; i++)
		{			
			RBF_initial(observedValue[i], Position[i], Node, DistParam);
//			magisland(observedValue[i], Position[i]);
		}
		delete[] Node;
		delete[] DistParam;
	}
/******	find the maximun or minum observed value	***************************************/
	for (d = 0; d < Dim; d++)
	{
		abs_max_value[d] = abs(observedValue[0][d]);
		abs_min_value[d] = abs(observedValue[0][d]);
	}
	for (i = 1; i < Num_Crafts*Num_obserPosit; i++)
	{
		for (d = 0; d < Dim; d++)
		{
			if (abs(observedValue[i][d]) < abs_min_value[d])
				abs_min_value[d] = abs(observedValue[i][d]);
			if (abs(observedValue[i][d]) > abs_max_value[d])
				abs_max_value[d] = abs(observedValue[i][d]);
		}
	}
	ofstream dataout("observeddata.dat");							// write position and magnetic filed value
//	dataout << "data generating from c++ program not really from observation" << endl;
//	dataout << "Position,     observerdValue" << endl;
	write(Position, Num_Crafts*Num_obserPosit, dataout);
	write(observedValue, Num_Crafts*Num_obserPosit, dataout);
	dataout.close();
/****** recycle dynacim memory and stop the program	****************************************/
	delete[] CraftLeap;
	delete[] CraftsInterval;
	delete[] flying_posit;
	delete[] starting_posit;
	return Num_Crafts*Num_obserPosit;
}
void CoordinatesTransform(double observedValue[][Dim], double CraftsInterval[][Dim],
	int *Num_Crafts, double CraftLeap[][Dim], int *Num_obserPosit, ifstream& filename)
{
/******	a lot to do ******/
	*Num_Crafts = GroupNumber;
	*Num_obserPosit = M;
}
double tanhtofit(Point& var)
{
	double value;
	value = tanh(var.getx() / Length_Scale);			// can be modified
	return value;
}
void magisland(double* value, Point& var)
{
	double x, y, z;
	double oldx, oldy, oldz;							// old coordinates
	double Bx, By, Bz;
	double delta_psi = island_Magnit;					// delta_psi = island_Magnit*cos(kx*x)*cos(ky*y)
	double kx = 2 * 3.141592653589793 / (3 * Length_Scale);
	double ky = 3.141592653589793 / (3 * Length_Scale);
	if (Dim > 1)
	{
		x = var.getx(); y = var.gety();
/****** rotate x-y to older coordinates ******/
		oldx = x * cos(degreez) - y * sin(degreez);
		oldy = x * sin(degreez) + y * cos(degreez);
/****** rotate oldy-z to older coordinates in 3D	******/	
		if (Dim > 2)
		{
			z = var.getz();
			oldz = oldy * sin(degreex) + z * cos(degreex);
			oldy = oldy * cos(degreex) - z * sin(degreex);
		}
/****** compute magnetic field value in oldx-lody-oldz	******/			
		Bx = tanh(oldy / Length_Scale);
		Bx += -ky*delta_psi*cos(kx*oldx)*sin(ky*oldy);
		By = kx*delta_psi*sin(kx*oldx)*cos(ky*oldy);
		if (Dim > 2)
			Bz = 0.5;			
/******* rotate vector from oldy-oldz to oldy-z	******/
		if (Dim > 2)
		{
			value[1] = By * cos(degreex) + Bz * sin(degreex);
			value[2] = -By * sin(degreex) + Bz * cos(degreex);
		}
		else
			value[1] = By;
/******* rotate vector from x-oldy to x-y	******/
		value[0] = Bx * cos(degreez) + value[1] * sin(degreez);
		value[1] = Bx * (-sin(degreez)) + value[1] * cos(degreez);
	}	
}
void RBF_initial(double* value, Point& var, Point* Node, Point* DistParam)
{
	BasicFunction *Chi = new BasicFunction[N]();
	oncefit fitting;
	double* Alpha = new double[N_Alpha]();
	int i;
	for (i = 0; i < N; i++)
		Chi[i].specify(Node[i], DistParam[i]);
	for (i = 0; i < N_Alpha; i++)
		Alpha[i] = i+10;
	fitting.assignParameter(Alpha);
	fitting.getValue(value, var, Chi);
	delete[] Alpha;
	delete[] Chi;
}


