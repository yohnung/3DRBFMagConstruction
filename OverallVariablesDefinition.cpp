#include "MacroAndMostUsedLibrary.h"
/******	related to space craft, when read observed value form real craft data file, these are useless	******/
double M = 21;				// number of flying position, should be read from other file
/******	flying 5 space-crafts, generating position as	******/
/******			        **3 					******/
/******			    **8		**7					******/
/******			**4     **5     **2				******/
/******			    **9		**6					******/
/******			        **1						******/
double GroupNumber = 4;		// number of lines or spacecrafts, maximum 4
double startx = 3.2;			// starting position's x component
double starty = -27;			// starting position's y component
double startz = -2;		// starting position's z component
double Interval12x = 2;	// x distance between 1st and 2nd spacecraft
double Interval12y = 0;	// y distance between 1st and 2nd spacecraft
double Interval12z = 2;	// z distance between 1st and 2nd spacecraft
double Interval13x = 0;		// x distance between 1st and 3rd spacecraft
double Interval13y = 0;		// y distance between 1st and 3rd spacecraft
double Interval13z = 4;	// z distance between 1st and 3rd spacecraft
double Interval14x = -2;	// x distance between 1st and 4th spacecraft
double Interval14y = 0;	// y distance between 1st and 4th spacecraft
double Interval14z = 2;	// z distance between 1st and 4th spacecraft

double Interval15x = 0.;	// x distance between 1st and 5th spacecraft
double Interval15y = 0;		// y distance between 1st and 5th spacecraft
double Interval15z = .5;	// z distance between 1st and 5th spacecraft

double Interval16x = 0.25;	// x distance between 1st and 6th spacecraft
double Interval16y = 0;		// y distance between 1st and 6th spacecraft
double Interval16z = 0.25;	// z distance between 1st and 6th spacecraft
double Interval17x = 0.25;	// x distance between 1st and 7th spacecraft
double Interval17y = 0;		// y distance between 1st and 7th spacecraft
double Interval17z = .75;	// z distance between 1st and 7th spacecraft
double Interval18x = -0.25;	// x distance between 1st and 8th spacecraft
double Interval18y = 0;		// y distance between 1st and 8th spacecraft
double Interval18z = 0.75;	// z distance between 1st and 8th spacecraft
double Interval19x = -0.25;	// x distance between 1st and 9th spacecraft
double Interval19y = 0;		// y distance between 1st and 9th spacecraft
double Interval19z = 0.25;	// z distance between 1st and 9th spacecraft

double Leapx = 0;			// x distance of a flying craft between adjacent time
double Leapy = 1;			// y distance of a flying craft between adjacent time 
double Leapz = 0;			// z distance of a flying craft between adjacent time
/******	model the magnetic field value, useless when read from file	******/
double Length_Scale = 3.;				// preset Model 'island Magnetic Field's typical length used in tanh(x/Length_Scale)
double island_Magnit = 0.5;				// related to 'island Magnetic Field's island width = sqrt(4*Length_Scale*island_Magnit)
double degreez = Pi * 10/ 180;			// 10 degree rotate coordinates of z-axis
double degreey = Pi / 4;				// 45 degree rotate corrdinate of x-axis
/******	model the magnetic field value using 'Separator model'	******/
double Bmagnit = 1;						// magnitude of magnetic field
double Current_along_Separator = 2*sqrt(8);	// current flow along separator line
double Null_Posit = 1;					// null's positions, corresponding to +/- value
double Null_OffsetY = -4;
double Cur_Width = 10.;
/****** add some random error to the observed value		******/
int	   random_switch = 0;				// '1' for on; '0' for off;
double Random_Mean = 0.;
double Random_Sigma = 0.5 * Bmagnit;
double Random_Bound = 0.2 * Bmagnit;

int pos_rand_switch = 0;				// random position on when '1'
double xPos_Rand_Mean = 0;
double xPos_Rand_Sigma = 0.5 * Leapy;
double xPos_Rand_Bound = 0.2 * Leapy;
double yPos_Rand_Mean = 0;
double yPos_Rand_Sigma = 0.75 * Leapy;	// almost uniform
double yPos_Rand_Bound = 0.5 * Leapy;	// y-direction is the most random position
double zPos_Rand_Mean = 0;
double zPos_Rand_Sigma = 0.5 * Leapy;
double zPos_Rand_Bound = 0.2 * Leapy;
/******	choose what model you want use	******/
char   model = CURRENTSHEET;			// SEPARATORMODEL for separator model; RBFMODEL for RBF-model; DEFAULT for no model
/******	Parameter used in solver	******************************/
double cutoff = 3e2;			// When we solve the problem using SVD, omit singular value using this parameter
double tolerance = 0.5;		// when B is lower than tolerance*max_value, the related equation will be wholely divided by B 

