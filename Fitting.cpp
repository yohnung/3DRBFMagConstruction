/*
	This is a Code Source File of 'Reconstrction of Magnetic Field'.
	Basis Functions are Radial Basic Function \chi denoted `Chi` in code.
	Specifically for 1D:
				B(x) = \sum_j \alpha_j * \chi(x, R_j ,D_j)
				for 2D:
				B\vec = \nabla \psi ]times e\hat_z
						with \psi(x\vec) = \sum_j \alpha_j * \chi(x\vec, R\vec_j, D_j)
				for 3D:
				B\vec = \nabla \times (\psi_1 * r\vec)
						+\nabla \times \nabla \times (\psi_2 * r\vec)
						with \psi_1,2(x\vec) = \sum_j \alpha_1,j * \chi(x\vec, R\vec_j, D_j).
				R\vec_j above means jth Node's position.
	Here in code, we construct a web of Nodes containing N Nodes, and correspondingly
	there are
		N alphas in 1D or 2D; 2*N alphas in 3D.

	Meaning of constants or variables' name:
		Dim			: Dimension of the problem
		Nx,Ny,Nz	: Number of Nodes in x,y,z direction
		N			: Number of Nodes in total
		N_Alpha		: Number of Alphas in total. = N in 1D/2D; = 2N in 3D
		max_length	: maximum number of observation point
		tolerance	: Using in specification of Matrix A and B when solving A*alpha=B
*/
#include <mkl.h>
#include "MacroAndMostUsedLibrary.h"
using namespace std;
#include "ClassClarification.h"
#include "ProcedureClarification.h"
#include "ReadingData.h"

int main()
{
	Point* Position = new Point[max_length]();		// observed position
	double(*observedValue)[Dim] = new double[max_length][Dim]();				// observed magnetic field value
	int Number_obserPosit;							// Number of observed position
	Point Node[N];									// Nodal position, x-coordinate change first, then y, finnaly z; index=k*Ny*Nx+j*Ny+i
	Point DistParam[N];								// Distance Parameter used in Basic Function
	BasicFunction Chi[N];							// Corresponding to different Nodal Function
	double Alpha[N_Alpha];							// for 1D/2D N_alpha = N; 
	ofstream Resultout("RBF_Length_Node_Dist_Alpha.dat");		// write out Node, Distance Parameter, Alpha
/******	specify craft's position relate to whole Structure, magnetic field value and number of observation data	******/
	double* max_posit = new double[Dim]();							// max_position_x, max_position_y, max_position_z;
	double* min_posit = new double[Dim]();							// min position
	double* abs_max_value = new double[Dim]();						// max observed value
	double* abs_min_value = new double[Dim]();						// minimum observed value
	Number_obserPosit = read(Position, observedValue, max_posit, 
		min_posit, abs_max_value, abs_min_value);					// reading by generating data
//	Number_obserPosit = read(Position, observedValue, max_posit,
//		min_posit, abs_max_value, abs_min_value, filename);			// reading from real observed value
/******	construct a web of Nodes and specify Basic Functions	******************************************************/
	constructNodesWeb(Node, DistParam, max_posit, min_posit);		// instruct web of node
	Resultout << Lx << " " << Ly << " " << Lz << endl;
	write(Node, N, Resultout);
	write(DistParam, N, Resultout);
	for (int i = 0; i < N; i++)
		Chi[i].specify(Node[i], DistParam[i]);						// specify Basic Function
/******	define Matrix A and Y to solve Alpha of A*Alpha=Y, using LU, QR and SVD method	******************************/
	double (*A)[N_Alpha] = new double[Dim * Number_obserPosit][N_Alpha]();	// define Matrix A, later in used in solving A* alpha = Y
	double* Y = new double[Dim * Number_obserPosit]();						// observed value
	specifyAY(A, Y, observedValue, Position, Number_obserPosit, Chi);		// specify A and Y using RBF and observed value
	double acond = 1;														// condition numver of A^T * A
	acond = LinearLUSolver(Alpha, A, Y, Number_obserPosit);					// LU solving A * alpha = Y
	LinearQRSolver(Alpha, A, Y, Number_obserPosit, abs_max_value);
	LinearSVDSolver(Alpha, A, Y, Number_obserPosit, abs_max_value, acond);	// SVD to solve the problem
	write(Alpha, N_Alpha, Resultout);										// output fittiing parameter
/******	using optimized 'Alpha's to reconstruct magnetic filed		**************************************************/
	double(*fittedValue)[Dim] = new double[Number_obserPosit][Dim]();
	oncefit fitting;
	fitting.assignParameter( Alpha );										
	for (int i = 0; i < Number_obserPosit; i++)
		fitting.getValue(fittedValue[i], Position[i], Chi);
/****** close file, recycle dynamic memory and terminate the code	**************************************************/
	Resultout.close();
	delete[] Y;
	delete[] A;
	delete[] Position;
	delete[] observedValue;
	delete[] max_posit;
	delete[] min_posit;
	delete[] fittedValue;
	delete[] abs_max_value;
	delete[] abs_min_value;
	return 0;
}