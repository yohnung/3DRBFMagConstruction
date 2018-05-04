// Macro definition and most used library header file
#include <iostream>
#include <fstream>
/******	Configurating problem	*********************************/
#define Dim 3				// dimension
#define Nx 10				// number of Nodes in x direction, must be at least 1
#define Ny 2				// number of Nodes in y direction, must be at least 1
#define Nz 2				// number of Nodes in z direction, must be at least 1
#define N  Nx*Ny*Nz			// number of Nodes in total
#define GNx 40				// number of Grid in x direction
#define GNy 100				// number of Grid in y direction
#define GNz 100				// number of Grid in z direction
#define GN GNx*GNy*GNz		// number of Grids in total
#define N_Alpha N*2			// Number of 'alpha's; for 1D/2D, N_Alpha=N; for 3D N_Alpha=2N
#define max_length 1000		// maximum number of observed positions
/****** model magnetic field type	******/
#define RBFMODEL 'r'		// use RBF-mode to generate a model field
#define SEPARATORMODEL 's'	// use separator-model
#define CURRENTSHEET  'c'	// use current sheet in 1d or magisland in 2d
#define DEFAULT 'd'			// do not use model field, instead, use reading data from observation
#define Pi	3.141592653589793
