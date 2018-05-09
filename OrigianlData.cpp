#include<math.h>
#include<iomanip>
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
	CraftsInterval[3][0] = Interval15x;								// Interval[3] for 1st and 5th
	CraftsInterval[4][0] = Interval16x;								// Interval[4] for 1st and 6th
	CraftsInterval[5][0] = Interval17x;								// Interval[5] for 1st and 7th
	CraftsInterval[6][0] = Interval18x;								// Interval[6] for 1st and 8th
	CraftsInterval[7][0] = Interval19x;								// Interval[7] for 1st and 9th

	if (Dim > 1)
	{
		starting_posit[1] = starty;
		flying_posit[1] = starty;
		CraftsInterval[0][1] = Interval12y;
		CraftsInterval[1][1] = Interval13y;
		CraftsInterval[2][1] = Interval14y;
		CraftsInterval[3][1] = Interval15y;
		CraftsInterval[4][1] = Interval16y;
		CraftsInterval[5][1] = Interval17y;
		CraftsInterval[6][1] = Interval18y;
		CraftsInterval[7][1] = Interval19y;
	}
	if (Dim > 2)
	{
		starting_posit[2] = startz;
		flying_posit[2] = startz;
		CraftsInterval[0][2] = Interval12z;
		CraftsInterval[1][2] = Interval13z;
		CraftsInterval[2][2] = Interval14z;
		CraftsInterval[3][2] = Interval15z;
		CraftsInterval[4][2] = Interval16z;
		CraftsInterval[5][2] = Interval17z;
		CraftsInterval[6][2] = Interval18z;
		CraftsInterval[7][2] = Interval19z;

	}
	Position[0].specify(flying_posit);
	default_random_engine generator(time(NULL));
	uniform_real_distribution<double> distribution(0.1, 0.3);		// mean value is 0.25; can be modified
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
	int Num_obserPosit = M;
	double x[3];
	double valuex, valuey, valuez;
	double Hour, Minute, Second, tmp;
	int num;
	num = 0;
	while (datain.good())
	{
		datain >> Hour >> Minute >> Second;
		datain >> valuex >> valuey >> valuez >> tmp >> x[0] >> x[1] >> x[2];
		Position[num].specify(x); 
		observedValue[num][0] = valuex; 
		observedValue[num][1] = valuey; 
		observedValue[num][2] = valuez;
		datain >> valuex >> valuey >> valuez >> tmp >> x[0] >> x[1] >> x[2];
		Position[Num_obserPosit + num].specify(x); 
		observedValue[Num_obserPosit + num][0] = valuex; 
		observedValue[Num_obserPosit + num][1] = valuey; 
		observedValue[Num_obserPosit + num][2] = valuez;
		datain >> valuex >> valuey >> valuez >> tmp >> x[0] >> x[1] >> x[2];
		Position[2 * Num_obserPosit + num].specify(x);
		observedValue[2 * Num_obserPosit + num][0] = valuex;
		observedValue[2 * Num_obserPosit + num][1] = valuey;
		observedValue[2 * Num_obserPosit + num][2] = valuez;
		datain >> valuex >> valuey >> valuez >> tmp >> x[0] >> x[1] >> x[2];
		Position[3 * Num_obserPosit + num].specify(x);
		observedValue[3 * Num_obserPosit + num][0] = valuex;
		observedValue[3 * Num_obserPosit + num][1] = valuey;
		observedValue[3 * Num_obserPosit + num][2] = valuez;
		num = num + 1;
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
	double kx = Pi / (3 * Length_Scale);
	double ky = 2 * Pi / (3 * Length_Scale);
	double kz = 2 * Pi / (3 * Length_Scale);
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
		By = tanh(oldx / Length_Scale);			// B\vec = e\hat_z \times \nabla\psi
		Bx = ky*delta_psi*cos(kx*oldx)*sin(ky*oldy);
		By += -kx*delta_psi*sin(kx*oldx)*cos(ky*oldy);
		value[0] = Bx;
		if (Dim > 2)			//  y component is guide field
		{
			Bz = tanh(oldx / Length_Scale);
			Bx =   kz*delta_psi*cos(kx*oldx)*sin(kz*oldz);
			Bz += -kx*delta_psi*sin(kx*oldx)*cos(kz*oldz);
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
void ModelField(int num, Point* Node, Point* DistParam, Point* Position, double Value[][Dim])
{
	int i;
	if (model == SEPARATORMODEL)
		SeparatorField(num, Position, Value);	// use Radial Basic Function to re-model magnetic field
	else if (model == RBFMODEL)
		RBFModelField(num, Node, DistParam, Position, Value);	// use Separator Model to model a magnetic field
	else if (model == CURRENTSHEET)
	{
		if (Dim == 1)
			for (i = 0; i < num; i++)
				Value[i][0] = tanhtofit(Position[i]);			// f(x,y,z) to be fitted, specified manually
		if (Dim > 1)
			for (i = 0; i < num; i++)
				magisland(Value[i], Position[i]);				// magnetic island configuration
	}
	else
		cout << "You should read data from real observation file, do you really do that?" << endl;
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
void SeparatorField(int num, Point* Position, double Value[][Dim])
{
/******	use Separator Model to model a magnetic field, please refer to Pontin[2011], Advances in Space Research	******/
	 double Bx, By, Bz;
	 double x, y, z;
	 int i, j;
	 
	 for (i = 0; i < num; i++)
	 {
		 x = Position[i].getx();
		 y = Position[i].gety();
		 z = Position[i].getz();

		 y = y - Null_OffsetY;			// This is means that magneti-null is located at (Null_OffsetY +/- Null_Posit)
		 Bx = Bmagnit * (x * (y - 3 * Null_Posit)
			 + 0.5 * Current_along_Separator * (-z * exp(-(x*x + z*z) / ( 2 * Cur_Width * Cur_Width)))); 
		 By = Bmagnit * (Null_Posit * Null_Posit - y * y + 0.5 * (x * x + z * z));
		 Bz = Bmagnit * (z * (y + 3 * Null_Posit)
			 + 0.5 * Current_along_Separator * ( x * exp(-(x*x + z*z) / ( 2 * Cur_Width * Cur_Width))));
		 Value[i][0] = Bx;
		 Value[i][1] = By;
		 Value[i][2] = Bz;
	 }
	 write_script_tecplot();
}
void write_script_tecplot()
{
/******	create 5 circle geometries which are the origin of streamtrace to be plotted in tecplot	******/
	int nn = 7;							// 7 points
	double rr = 0.1 * Null_Posit;
	double* theta = new double[nn];
	double tmptheta;
	double x, y, z;
	int i, j;
	for (i = 0; i < nn; i++)
		theta[i] = 2. * Pi / nn * (0.5 + i);
	ofstream fileID("streamtrace_circle_geometry_for3D.mcr");
	fileID << "#!MC 1410" << endl;
	fileID << "$!VarSet |MFBD| = 'C:\\Program Files\\Tecplot\\Tecplot 360 EX 2016 R1' " << endl;
	default_random_engine generator(time(NULL));
	uniform_real_distribution<double> distribution(0., 0.5*(theta[1]-theta[0]));		// mean value is 0.25; can be modified
	tmptheta = distribution(generator);
	for (i = 0; i < nn; i++)			// C1 circle
	{
		x = rr * cos(theta[i] + tmptheta);
		z = rr * sin(theta[i] + tmptheta);
		y = Null_OffsetY;
		fileID << " $!STREAMTRACE ADD\n  STREAMTYPE = VOLUMELINE\n  STREAMDIRECTION = BOTH" << endl;
		fileID << "  STARTPOS{X=" << x << " Y=" << y << " Z=" << z << "}" << endl;
	}
	tmptheta = distribution(generator);
	for (i = 0; i < nn; i++)			// C2 circle
	{
		x = -.5;
		z = 0 + rr * cos(theta[i] + tmptheta);
		y = Null_OffsetY-Null_Posit + rr * sin(theta[i] + tmptheta);
		fileID << " $!STREAMTRACE ADD\n  STREAMTYPE = VOLUMELINE\n  STREAMDIRECTION = BOTH" << endl;
		fileID << "  STARTPOS{X=" << x << " Y=" << y << " Z=" << z << "}" << endl;;
	}
	tmptheta = distribution(generator);
	for (i = 0; i < nn; i++)			// C3 circle
	{
		x = .5;
		z = 0 + rr * cos(theta[i] + tmptheta);
		y = Null_OffsetY-Null_Posit + rr * sin(theta[i] + tmptheta);
		fileID << " $!STREAMTRACE ADD\n  STREAMTYPE = VOLUMELINE\n  STREAMDIRECTION = BOTH" << endl;
		fileID << "  STARTPOS{X=" << x << " Y=" << y << " Z=" << z << "}" << endl;
	}
	tmptheta = distribution(generator);
	for (i = 0; i < nn; i++)			// C4 circle
	{
		z = -.5;
		x = 0 + rr * cos(theta[i] + tmptheta);
		y = Null_OffsetY+Null_Posit + rr * sin(theta[i] + tmptheta);
		fileID << " $!STREAMTRACE ADD\n  STREAMTYPE = VOLUMELINE\n  STREAMDIRECTION = BOTH" << endl;
		fileID << "  STARTPOS{X=" << x << " Y=" << y << " Z=" << z << "}" << endl;
	}
	tmptheta = distribution(generator);
	for (i = 0; i < nn; i++)			// C5 circle
	{
		z = .5;
		x = 0 + rr * cos(theta[i] + tmptheta);
		y = Null_OffsetY+Null_Posit + rr * sin(theta[i] + tmptheta);
		fileID << " $!STREAMTRACE ADD\n  STREAMTYPE = VOLUMELINE\n  STREAMDIRECTION = BOTH" << endl;
		fileID << "  STARTPOS{X=" << x << " Y=" << y << " Z=" << z << "}" << endl;
	}
	fileID << "$!RemoveVar |MFBD|";
}
