/******	useless when the code can read real craft's data	******/
extern double M;				// number of flying position, useless when read from file
extern double GroupNumber;		// number of lines or spacecrafts
extern double startx; 	// starting position's x component
extern double starty;	// starting position's y component
extern double startz;	// starting position's z component
extern double Interval12x;		// x distance between 1st and 2nd spacecraft
extern double Interval12y;		// y distance between 1st and 2nd spacecraft
extern double Interval12z;		// z distance between 1st and 2nd spacecraft
extern double Interval13x;		// x distance between 2nd and 3rd spacecraft
extern double Interval13y;		// y distance between 2nd and 3rd spacecraft
extern double Interval13z;		// z distance between 2nd and 3rd spacecraft
extern double Interval14x;		// x distance between 1st and 2nd spacecraft
extern double Interval14y;		// y distance between 1st and 2nd spacecraft
extern double Interval14z;		// z distance between 1st and 2nd spacecraft
extern double Leapx;			// x distance of a flying craft between adjacent time
extern double Leapy; 			// y distance of a flying craft between adjacent time
extern double Leapz;			// z distance of a flying craft between adjacent time
extern double Length_Scale;		// preset Model 'island Magnetic Field's typical length used in tanh(x/Length_Scale)
extern double island_Magnit;	// related to 'island Magnetic Field's island width = sqrt(4*Length_Scale*island_Magnit)
extern double degreez;			// 45 degree rotate coordinates of z-axis
extern double degreey;			// 10 degree rotate corrdinate of x-axis#pragma once
/******	Parameter used in solver	******************************/
extern double cutoff;		// When we solve the problem using SVD, omit singular value using this parameter
extern double tolerance;	// when B is lower than tolerance*max_value, the related equation will be wholely divided by B 
