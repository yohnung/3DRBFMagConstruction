// Macro definition and most used library header file
#include <iostream>
#include <fstream>
/******	Configurating problem	*********************************/
#define Dim 2				// dimension
#define Nx 20				// number of Nodes in x direction, must be at least 1
#define Ny 3				// number of Nodes in x direction, must be at least 1
#define Nz 1				// number of Nodes in x direction, must be at least 1
#define N  Nx*Ny*Nz			// number of Nodes in total
#define N_Alpha N*1			// Number of 'alpha's; for 1D/2D, N_Alpha=N; for 3D N_Alpha=2N
#define max_length 1000		// maximum number of observed positions
#define Pi	3.141592653589793
