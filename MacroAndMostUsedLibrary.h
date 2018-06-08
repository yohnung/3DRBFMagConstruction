// Macro definition and most used library header file
#include <iostream>
#include <fstream>
/******	Configurating problem	*********************************/
#define Dim 3				// dimension
#define Nx 20				// number of Nodes in x direction, must be at least 1
#define Ny 2				// number of Nodes in y direction, must be at least 1
#define Nz 2				// number of Nodes in z direction, must be at least 1
#define N  Nx*Ny*Nz			// number of Nodes in total
#define GNx 101				// number of Grid in x direction
#define GNy 41				// number of Grid in y direction
#define GNz 41				// number of Grid in z direction
#define GN GNx*GNy*GNz		// number of Grids in total
#define N_Alpha N*2			// Number of 'alpha's; for 1D/2D, N_Alpha=N; for 3D N_Alpha=2N
#define max_length 1500		// maximum number of observed positions
/****** model magnetic field type	******/
#define RBFMODEL 'r'			// use RBF-mode to generate a model field
#define SEPARATORMODEL1 's'		// use separator-model-1
#define SEPARATORMODEL2 't'		// use separator-model-2
#define CURRENTSHEET  'c'		// use current sheet in 1d or magisland in 2d
#define DEFAULT 'd'				// do not use model field, instead, use reading data from observation
#define Pi	3.141592653589793
