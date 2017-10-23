#include<random>
#include<time.h>
#include "MacroAndMostUsedLibrary.h"
#include "OverallVariablesDeclarification.h"
using namespace std;
#include "ClassClarification.h"
#include "ProcedureClarification.h"
#include "ReadingData.h"

int read(Point* Position, double observedValue[][Dim])
{
/******	generate the magnetic field value	*******************************/
	int Num_obserPosit = M;										// Number of observed position
	int Num_Crafts = GroupNumber;								// Number of Crafts
	double* starting_posit = new double[Dim]();						// starting position (x,y,z)
	double* flying_posit = new double[Dim]();						// flying position
	double(*CraftsInterval)[Dim] = new double[10][Dim]();			// distance vector between maximum 10 different space-craft
	double(*CraftLeap)[Dim] = new double[max_length][Dim]();		// space-craft's interval between two adjacent time

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
	CraftsInterval[0][0] = Interval12x;								// Interval[0] means distance vector between 1st and 2nd
	CraftsInterval[1][0] = Interval13x;								// Interval[1] for 1st and 3rd
	CraftsInterval[2][0] = Interval14x;								// Interval[2] for 1st and 4th
	if (Dim > 1)
	{
		starting_posit[1] = starty;
		flying_posit[1] = starty;
		CraftsInterval[0][1] = Interval12y;
		CraftsInterval[1][1] = Interval13y;
		CraftsInterval[2][1] = Interval14y;
	}
	if (Dim > 2)
	{
		starting_posit[2] = startz;
		flying_posit[2] = startz;
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
				flying_posit[d] = starting_posit[d] + CraftsInterval[g - 1][d];
			Position[g*Num_obserPosit].specify(flying_posit);
		}
		for (i = 1; i < Num_obserPosit; i++)						// for a specific craft
		{
			for (d = 0; d < Dim; d++)
				flying_posit[d] += CraftLeap[i][d];	// + distribution(generator);
			Position[g*Num_obserPosit + i].specify(flying_posit);
		}
	}
/******	model the observed value	******/
	if (Dim == 1)
		for (i = 0; i < Num_Crafts*Num_obserPosit; i++)
			observedValue[i][0] = tanhtofit(Position[i]);			// f(x,y,z) to be fitted, specified manually
	if (Dim > 1)
		for (i = 0; i < Num_Crafts*Num_obserPosit; i++)
			magisland(observedValue[i], Position[i]);				// magnetic island configuration
	
	
/****** recycle dynacim memory and stop the program	****************************************/
	delete[] CraftLeap;
	delete[] CraftsInterval;
	delete[] flying_posit;
	delete[] starting_posit;
	return Num_Crafts*Num_obserPosit;
}
int read(Point* Position, double observedValue[][Dim], ifstream& datain)
{
/******	simply read from datain the position of craft and observed magnetic field value	******/
	double x, y, z;
	double valuex, valuey, valuez;
	char ch;
	int i, num;
	for (i = 0; i < max_length; i++)
	{
		datain.get(ch);
		if (!ch)
			break;
		else
		{
			datain >> x >> y >> z >> valuex >> valuey >> valuez;
			num = i + 1;
		}
	}
	M = num;							// every craft's observed position's number
	GroupNumber = 4;					// craft's number
	return GroupNumber*M;
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
	double kx = 3.141592653589793 / (3 * Length_Scale);
	double ky = 2 * 3.141592653589793 / (3 * Length_Scale);
	double kz = 2 * 3.141592653589793 / (3 * Length_Scale);
	if (Dim > 1)	// By is guide field
	{
/****** rotate x-y to older coordinates ******/
		x = var.getx(); y = var.gety();
		oldx = x * cos(degreez) + y * sin(degreez);
		oldy = -x * sin(degreez) + y * cos(degreez);
		if (Dim > 2)
		{
/****** rotate oldx-z to older coordinates in 3D	******/
			z = var.getz();
			oldz = oldx * sin(degreey) + z * cos(degreey);
			oldx = oldx * cos(degreey) - z * sin(degreey);
		}
/****** compute magnetic field value in oldx-lody-oldz	******/	
		By = tanh(oldx / Length_Scale);
		Bx = ky*delta_psi*cos(kx*oldx)*sin(ky*oldy);
		By += -kx*delta_psi*sin(kx*oldx)*cos(ky*oldy);
		value[0] = Bx;
		if (Dim > 2)			//  y component is guide field
		{
			Bz = tanh(oldx / Length_Scale);
			Bx = -kz*delta_psi*cos(kx*oldx)*sin(kz*oldz);
			Bz += kx*delta_psi*sin(kx*oldx)*cos(kz*oldz);
			By = 0.5;
/******* rotate vector from oldx-oldz to oldx-z	******/
			value[0] = Bx * cos(degreey) + Bz * sin(degreey);
			value[2] = -Bx * sin(degreey) + Bz * cos(degreey);
		}
/******* rotate vector from oldx-oldy to x-y	******/
		value[1] = value[0] * sin(degreez) + By * cos(degreez);
		value[0] = value[0] * cos(degreez) - By * sin(degreez);
	}	
}
void RBFModelField(int num, Point* Node, Point* DistParam, Point* Position, double observedValue[][Dim])
{
	BasicFunction *Chi = new BasicFunction[N]();
	oncefit fitting;
	double* Alpha = new double[N_Alpha]();
	int i, j;
	for (i = 0; i < N; i++)
		Chi[i].specify(Node[i], DistParam[i]);
//	for (i = 0; i < N_Alpha; i++)
//		Alpha[i] = i + 10;
	Alpha[0] = 15;
	Alpha[3] = -7;
	Alpha[5] = 23;
	fitting.assignParameter(Alpha);
	for (i = 0; i < num; i++)
		fitting.getValue(observedValue[i], Position[i], Chi);
	delete[] Alpha;
	delete[] Chi;
}


