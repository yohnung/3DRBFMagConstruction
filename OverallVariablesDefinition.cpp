#include "MacroAndMostUsedLibrary.h"
/******	related to space craft, when read observed value form real craft data file, these are useless	******/
double M = 21;				// number of flying position, useless when read from file
double GroupNumber = 4;		// number of lines or spacecrafts, maximum 4
double startx = 0;			// starting position's x component
double starty = -10;		// starting position's y component
double startz = -2;			// starting position's z component
double Interval12x = 2;	// x distance between 1st and 2nd spacecraft
double Interval12y = 0;	// y distance between 1st and 2nd spacecraft
double Interval12z = 2;	// z distance between 1st and 2nd spacecraft
double Interval13x = 0;	// x distance between 2nd and 3rd spacecraft
double Interval13y = 0;	// y distance between 2nd and 3rd spacecraft
double Interval13z = 4;	// z distance between 2nd and 3rd spacecraft
double Interval14x = -2;	// x distance between 1st and 2nd spacecraft
double Interval14y = 0;	// y distance between 1st and 2nd spacecraft
double Interval14z = 2;	// z distance between 1st and 2nd spacecraft
double Leapx = 0;			// x distance of a flying craft between adjacent time
double Leapy = 1;		// y distance of a flying craft between adjacent time
double Leapz = 0;			// z distance of a flying craft between adjacent time
/******	model the magnetic field value, useless when read from file	******/
double Length_Scale = 3.;				// preset Model 'island Magnetic Field's typical length used in tanh(x/Length_Scale)
double island_Magnit = 0.5;			// related to 'island Magnetic Field's island width = sqrt(4*Length_Scale*island_Magnit)
double degreez = Pi * 10/ 180;	// 10 degree rotate coordinates of z-axis
double degreey = Pi / 4;				// 45 degree rotate corrdinate of x-axis
/******	choose what model you want use	******/
char   model = SEPARATORMODEL;			// SEPARATORMODEL for separator model; RBFMODEL for RBF-model; DEFAULT for no model
/******	Parameter used in solver	******************************/
double cutoff = -1e5;		// When we solve the problem using SVD, omit singular value using this parameter
double tolerance = 0.5;		// when B is lower than tolerance*max_value, the related equation will be wholely divided by B 

