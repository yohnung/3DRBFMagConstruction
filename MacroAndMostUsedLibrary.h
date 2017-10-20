// Macro definition and most used library header file
#include <iostream>
#include <fstream>
/******	Configurating problem	*********************************/
#define Dim 3				// dimension
#define Nx 10				// number of Nodes in x direction, must be at least 1
#define Ny 3				// number of Nodes in x direction, must be at least 1
#define Nz 3				// number of Nodes in x direction, must be at least 1
#define N  Nx*Ny*Nz			// number of Nodes in total
#define N_Alpha N*2			// Number of 'alpha's; for 1D/2D, N_Alpha=N; for 3D N_Alpha=2N
#define max_length 1000		// maximum number of observed positions
#define startx -2		// starting position's x component
#define starty -2		// starting position's y component
#define startz -5			// starting position's z component
#define Lx 0.4				// Radial Basic Function's scale length in x direction
#define Ly 2				// Radial Basic Function's scale length in y direction
#define Lz 1				// Radial Basic Function's scale length in z direction
/******	Parameter used in solver	******************************/
#define cutoff  1e+33		// When we solve the problem using SVD, omit singular value using this parameter
#define tolerance 0.5		// when B is lower than tolerance*max_value, the related equation will be wholely divided by B 
/******	useless when the code can read real craft's data	******/
#define M 21					// number of flying position, useless when read from file
#define GroupNumber 3		// number of lines or spacecrafts
#define Interval12x 0		// x distance between 1st and 2nd spacecraft
#define Interval12y 2		// y distance between 1st and 2nd spacecraft
#define Interval12z 0		// z distance between 1st and 2nd spacecraft
#define Interval13x	0		// x distance between 2nd and 3rd spacecraft
#define Interval13y	4		// y distance between 2nd and 3rd spacecraft
#define Interval13z	0		// z distance between 2nd and 3rd spacecraft
#define Interval14x 0		// x distance between 1st and 2nd spacecraft
#define Interval14y 6		// y distance between 1st and 2nd spacecraft
#define Interval14z 0		// z distance between 1st and 2nd spacecraft
#define Leapx 0.2				// x distance of a flying craft between adjacent time
#define Leapy 0					// y distance of a flying craft between adjacent time
#define Leapz 0.2				// z distance of a flying craft between adjacent time
#define Length_Scale	3.		// preset Model 'island Magnetic Field's typical length used in tanh(x/Length_Scale)
#define island_Magnit	.5		// related to 'island Magnetic Field's island width = sqrt(4*Length_Scale*island_Magnit)
#define degreez		0//3.141592653589793/4						// 45 degree rotate coordinates of z-axis
#define degreex		0//3.141592653589793*10/180				// 10 degree rotate corrdinate of x-axis