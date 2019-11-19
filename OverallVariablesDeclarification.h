/******	useless when the code can read real craft's data	******/
extern double M;				// number of flying position, useless when read from file
extern double GroupNumber;		// number of lines or spacecrafts
extern double startx; 	// starting position's x component
extern double starty;	// starting position's y component
extern double startz;	// starting position's z component
extern double Interval12x;		// x distance between 1st and 2nd spacecraft
extern double Interval12y;		// y distance between 1st and 2nd spacecraft
extern double Interval12z;		// z distance between 1st and 2nd spacecraft
extern double Interval13x;		// x distance between 1st and 3rd spacecraft
extern double Interval13y;		// y distance between 1st and 3rd spacecraft
extern double Interval13z;		// z distance between 1st and 3rd spacecraft
extern double Interval14x;		// x distance between 1st and 4th spacecraft
extern double Interval14y;		// y distance between 1st and 4th spacecraft
extern double Interval14z;		// z distance between 1st and 4th spacecraft
extern double Interval15x;		// x distance between 1st and 5th spacecraft
extern double Interval15y;		// y distance between 1st and 5th spacecraft
extern double Interval15z;		// z distance between 1st and 5th spacecraft
extern double Interval16x;		// x distance between 1st and 6th spacecraft
extern double Interval16y;		// y distance between 1st and 6th spacecraft
extern double Interval16z;		// z distance between 1st and 6th spacecraft
extern double Interval17x;		// x distance between 1st and 7th spacecraft
extern double Interval17y;		// y distance between 1st and 7th spacecraft
extern double Interval17z;		// z distance between 1st and 7th spacecraft
extern double Interval18x;		// x distance between 1st and 8th spacecraft
extern double Interval18y;		// y distance between 1st and 8th spacecraft
extern double Interval18z;		// z distance between 1st and 8th spacecraft
extern double Interval19x;		// x distance between 1st and 9th spacecraft
extern double Interval19y;		// y distance between 1st and 9th spacecraft
extern double Interval19z;		// z distance between 1st and 9th spacecraft
extern double Leapx;			// x distance of a flying craft between adjacent time
extern double Leapy; 			// y distance of a flying craft between adjacent time
extern double Leapz;			// z distance of a flying craft between adjacent time
extern double Length_Scale;		// preset Model 'island Magnetic Field's typical length used in tanh(x/Length_Scale)
extern double island_Magnit;	// related to 'island Magnetic Field's island width = sqrt(4*Length_Scale*island_Magnit)
extern double degreez;			// 45 degree rotate coordinates of z-axis
extern double degreey;			// 10 degree rotate corrdinate of x-axis#pragma once
extern double Bmagnit;					// magnitude of magnetic field
extern double Current_along_Separator;	// current flow along separator line
extern double Null_Posit;				// null's positions, corresponding to +/- value
extern double Null_OffsetY;
extern double Cur_Width;
extern int	  random_switch;
extern double Random_Mean;
extern double Random_Sigma;
extern double Random_Bound;
extern int	  pos_rand_switch;			// random position on when '1'
extern double xPos_Rand_Mean;
extern double xPos_Rand_Sigma;
extern double xPos_Rand_Bound;
extern double yPos_Rand_Mean;
extern double yPos_Rand_Sigma;
extern double yPos_Rand_Bound;
extern double zPos_Rand_Mean;
extern double zPos_Rand_Sigma;
extern double zPos_Rand_Bound;
/******	choose what model you want use	******/
extern char   model;			// SEPARATORMODEL for separator model; RBFMODEL for RBF-model; DEFAULT for no model
/******	Parameter used in solver	******************************/
extern double cutoff;		// When we solve the problem using SVD, omit singular value using this parameter
extern double tolerance;	// when B is lower than tolerance*max_value, the related equation will be wholely divided by B 
