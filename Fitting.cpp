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
#include "OverallVariablesDeclarification.h"
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
	double TypicalLength[Dim];						// typical length used in radial basic function
	double Alpha[N_Alpha];							// for 1D/2D N_alpha = N; 
	ifstream datain("craft_observed_value");
	ofstream dataout("modified_obsered_data.dat");	// write position and magnetic filed value
	ofstream Resultout("RBF_Length_Node_Dist_Alpha.dat");		// write out Node, Distance Parameter, Alpha
/******	specify craft's position relate to whole Structure, magnetic field value and number of observation data	******/
	double* max_posit = new double[Dim]();									// max_position_x, max_position_y, max_position_z;
	double* min_posit = new double[Dim]();									// min position
	Number_obserPosit = read(Position, observedValue);						// reading by generating data whose configuration is a magnetic island
//	Number_obserPosit = read(Position, observedValue, datain);				// reading from real observed value
	double* Alpha1 = new double[Dim]();										// the first transformation degree
	CrdTrf2FlyDirec(Number_obserPosit, Position, observedValue, Alpha1,   	// Coordinates Transform to Flying Direction, make it x-direction
		max_posit, min_posit);
/******	construct a web of Nodes and specify Basic Functions	******************************************************/
	constructNodesWeb(max_posit, min_posit, Node, DistParam, TypicalLength);// construct web of node and specify Typical Length and max,min position
	CrdTrfFromFlyDirec(Number_obserPosit, Alpha1, Position, observedValue, 
		Node, DistParam, TypicalLength);									// coordinates transform form x-flying direction to original one
	double* abs_max_value = new double[Dim]();								// max observed value
	double* abs_min_value = new double[Dim]();								// minimum observed value
	double* Alpha2 = new double[Dim]();										// second transformation degree
	CrdTrf2MinVarDir(Number_obserPosit, Position, observedValue, Alpha2, 	// according to position and observe value to
		Node, DistParam, TypicalLength, abs_max_value, abs_min_value);		// specify a min varing direction, make it z-direction and redifine Lx, Ly, Lz
//	TypicalLength[0] = 0.4; TypicalLength[1] = 2;
	BasicFunction::assign_TypicalLength(TypicalLength);
/******	use RBF to model magnetic field and write position, value, node, distance parameter and Typical Length	*******/
	RBFModelField(Number_obserPosit, Node, DistParam, Position, observedValue);	// use Radial Basic Function to re-model magnetic field
//	dataout << "data modified by c++ program " << endl;	dataout << "Position,     observerdValue" << endl;
	write(Position, Number_obserPosit, dataout); write(observedValue, Number_obserPosit, dataout);
	Resultout << BasicFunction::Lx << " " << BasicFunction::Ly << " " << BasicFunction::Lz << " " << endl;
	write(Node, N, Resultout);	write(DistParam, N, Resultout);
/******	specify Radial Basic Function	*******/
	for (int i = 0; i < N; i++)
		Chi[i].specify(Node[i], DistParam[i]);
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
	delete[] fittedValue;
	delete[] Y;
	delete[] A;
	delete[] Alpha2;
	delete[] abs_min_value;
	delete[] abs_max_value;
	delete[] Alpha1;
	delete[] min_posit;
	delete[] max_posit;
	delete[] observedValue;
	delete[] Position;
	Resultout.close();
	dataout.close();
	datain.close();
	return 0;
}