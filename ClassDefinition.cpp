#include "MacroAndMostUsedLibrary.h"
using namespace std;
#include "ClassClarification.h"

Point::Point()
{
	Coordinate = new double[Dim]();
	Coordinate[0] = 0;
	if (Dim > 1)
		Coordinate[1] = 0;
	if (Dim > 2)
		Coordinate[2] = 0;
}
Point::Point(double* x)
{
	Coordinate = new double[Dim]();
	Coordinate[0] = x[0];
	if (Dim > 1)
		Coordinate[1] = x[1];
	if (Dim > 2)
		Coordinate[2] = x[2];
}
Point::Point(Point& x)
{
	Coordinate = new double[Dim]();
	Coordinate[0] = x.getx();
	if (Dim > 1)
		Coordinate[1] = x.gety();
	if (Dim > 2)
		Coordinate[2] = x.getz();
}
Point::~Point()
{
	delete[] Coordinate;
}
void Point::specify(double *x)
{
	Coordinate[0] = x[0];
	if (Dim > 1)
		Coordinate[1] = x[1];
	if (Dim > 2)
		Coordinate[2] = x[2];
}
double Point::getradialLength()
{
	radialLength = sqrt(Coordinate[0] * Coordinate[0]);
	if (Dim > 1)
		radialLength = sqrt(Coordinate[0] * Coordinate[0]
			+ Coordinate[1] * Coordinate[1]);
	if (Dim > 2)
		radialLength = sqrt(Coordinate[0] * Coordinate[0]
			+ Coordinate[1] * Coordinate[1]
			+ Coordinate[2] * Coordinate[2]);
	return radialLength;
}
double Point::getx()
{
	return Coordinate[0];
}
double Point::gety()
{
	if (Dim > 1)
		return Coordinate[1];
	else
	{
		cout << "Wrong, there is no y-component!" << endl;
		return 0;
	}
}
double Point::getz()
{
	if (Dim > 2)
		return Coordinate[2];
	else
	{
		cout << "Wrong, there is no z-componnet" << endl;
		return 0;
	}
}
void Point::assignx(double x)
{
	Coordinate[0] = x;
}
void Point::assigny(double y)
{
	if (Dim > 1)
		Coordinate[1] = y;
	else
		cout << "Wrong, there is no y-component" << endl;
}
void Point::assignz(double z)
{
	if (Dim > 2)
		Coordinate[2] = z;
	else
		cout << "Wrong, there is no z-component" << endl;
}
Point& Point::operator = (Point& b)
{
	Coordinate[0] = b.getx();
	if (Dim > 1)
		Coordinate[1] = b.gety();
	if (Dim > 2)
		Coordinate[2] = b.getz();
	return *this;
}
Point operator + (Point& a, Point& b)
{
	Point c;
	c.Coordinate[0] = a.Coordinate[0] + b.Coordinate[0];
	if (Dim > 1)
		c.Coordinate[1] = a.Coordinate[1] + b.Coordinate[1];
	if (Dim > 2)
		c.Coordinate[2] = a.Coordinate[2] + b.Coordinate[2];
	return c;
}
Point operator - (Point& a, Point& b)
{
	Point c;
	c.Coordinate[0] = a.Coordinate[0] - b.Coordinate[0];
	if (Dim > 1)
		c.Coordinate[1] = a.Coordinate[1] - b.Coordinate[1];
	if (Dim > 2)
		c.Coordinate[2] = a.Coordinate[2] - b.Coordinate[2];
	return c;
}
double operator * (Point& a, Point& b)
{
	double c;
	c = a.Coordinate[0] * b.Coordinate[0];
	if (Dim > 1)
		c = c + a.Coordinate[1] * b.Coordinate[1];
	if (Dim > 2)
		c = c + a.Coordinate[2] * b.Coordinate[2];
	return c;
}
Point operator * (double a, Point& b)
{
	Point c;
	c.Coordinate[0] = a * b.Coordinate[0];
	if (Dim > 1)
		c.Coordinate[1] = a * b.Coordinate[1];
	if (Dim > 2)
		c.Coordinate[2] = a * b.Coordinate[2];
	return c;
}
Point operator / (Point& a, double b)
{
	Point c;
	c.Coordinate[0] = 1 / b * a.Coordinate[0];
	if (Dim > 1)
		c.Coordinate[1] = 1 / b *a.Coordinate[1];
	if (Dim > 2)
		c.Coordinate[2] = 1 / b * a.Coordinate[2];
	return c;
}

BasicFunction::BasicFunction() {}
BasicFunction::BasicFunction(Point& R, Point& D)
{
	Node = R;
	DistParameter = D;
}
void BasicFunction::specify(Point& R, Point& D)
{
	Node = R;
	DistParameter = D;
}
double BasicFunction::getValue(Point& x)
{
	if (Dim == 1)
		Value = sqrt((x - Node)*(x - Node) + DistParameter*DistParameter);
	else if (Dim == 2)
	{
		Value = sqrt((x.getx() - Node.getx())*(x.getx() - Node.getx()) / (Lx*Lx)
			+ (x.gety() - Node.gety())*(x.gety() - Node.gety()) / (Ly*Ly)
			+ DistParameter*DistParameter);
	}
	else if (Dim == 3)
	{
		Value = sqrt((x.getx() - Node.getx())*(x.getx() - Node.getx()) / (Lx*Lx)
			+ (x.gety() - Node.gety())*(x.gety() - Node.gety()) / (Ly*Ly)
			+ (x.getz() - Node.getz())*(x.getz() - Node.getz()) / (Lz*Lz)
			+ DistParameter*DistParameter);
	}
	return Value;
}
void BasicFunction::getDerivative(Point& x, double* FirstDerivative,
						double SecondDerivative[][Dim], double *ScalarDerivative)
{
	Value = getValue(x);
	double xx, xy, xz, Nodex, Nodey, Nodez, Value3;
	xx = x.getx();	Nodex = Node.getx(); Value3 = Value*Value*Value;
	FirstDerivative[0] = (xx - Nodex) / Value;
	SecondDerivative[0][0] = 1 / Value - (xx - Nodex)*(xx - Nodex) / Value3;
	*ScalarDerivative = SecondDerivative[0][0];
	if (Dim > 1)
	{
		xy = x.gety();	Nodey = Node.gety();
		FirstDerivative[0] = FirstDerivative[0] / (Lx*Lx);
		FirstDerivative[1] = (xy - Nodey) / (Ly*Ly*Value);
		SecondDerivative[0][0] = 1 / (Lx*Lx*Value)
			- (xx - Nodex)*(xx - Nodex) / (Lx*Lx*Lx*Lx*Value3);
		SecondDerivative[0][1] = -(xx - Nodex)*(xy - Nodey) / (Lx*Lx*Ly*Ly*Value3);
		SecondDerivative[1][0] = SecondDerivative[0][1];
		SecondDerivative[1][1] = 1 / (Ly*Ly*Value) 
			- (xy - Nodey)*(xy - Nodey) / (Ly*Ly*Ly*Ly*Value3);
		*ScalarDerivative = SecondDerivative[0][0] + SecondDerivative[1][1];
	}
	if (Dim > 2)
	{
		xz = x.getz(); Nodez = Node.getz();
		FirstDerivative[2] = (xz - Nodez) / (Lz*Lz*Value);
		SecondDerivative[0][2] = -(xx - Nodex)*(xz - Nodez) / (Lx*Lx*Lz*Lz*Value3);
		SecondDerivative[1][2] = -(xy - Nodey)*(xz - Nodez) / (Ly*Ly*Lz*Lz*Value3);
		SecondDerivative[2][1] = SecondDerivative[1][2];
		SecondDerivative[2][0] = SecondDerivative[0][2];
		SecondDerivative[2][2] = 1 / (Lz*Lz*Value) 
			- (xz - Nodez)*(xz - Nodez) / (Lz*Lz*Lz*Lz*Value3);
		*ScalarDerivative += SecondDerivative[2][2];
	}
}
double BasicFunction::get_xDeriv(Point& x)
{
	Value = getValue(x);
	if (Dim == 1)		
		Derivative = (x.getx() - Node.getx()) / Value;
	else	// including 2D or higher condition
		Derivative = (x.getx() - Node.getx()) / (Lx * Lx * Value);
	return Derivative;
}
double BasicFunction::get_yDeriv(Point& x)
{
	if (Dim > 1)
	{
		Value = getValue(x);
		Derivative = (x.gety() - Node.gety()) / (Ly * Ly * Value);
		return Derivative;
	}
	else
		return 0;
}
double BasicFunction::get_zDeriv(Point& x)
{
	if (Dim > 2)
	{
		Value = getValue(x);
		Derivative = (x.getz() - Node.getz()) / (Lz * Lz * Value);
		return Derivative;
	}
	else
		return 0;
}
double BasicFunction::get_xxDeriv(Point& x)
{
	Value = getValue(x);
	if (Dim == 1)
		Derivative = 1 / Value 
		- (x.getx() - Node.getx())*(x.getx() - Node.getx()) / (Value*Value*Value);
	if (Dim > 1)
		Derivative = 1 / (Lx*Lx*Value) 
		- (x.getx() - Node.getx())*(x.getx() - Node.getx()) / (Lx*Lx*Lx*Lx*Value*Value*Value);
	return Derivative;
}
double BasicFunction::get_xyDeriv(Point& x)
{
	if (Dim > 1)
	{
		Value = getValue(x);
		Derivative = -(x.getx() - Node.getx())
			*(x.gety() - Node.gety()) / (Lx*Lx*Ly*Ly*Value*Value*Value);
		return Derivative;
	}
	else
		return 0;
}
double BasicFunction::get_yyDeriv(Point& x)
{
	if (Dim > 1)
	{
		Value = getValue(x);
		Derivative = 1 / (Ly*Ly*Value) 
			- (x.gety() - Node.gety())*(x.gety() - Node.gety()) / (Ly*Ly*Ly*Ly*Value*Value*Value);
		return Derivative;
	}
	else
		return 0;
}
double BasicFunction::get_xzDeriv(Point& x)
{
	if (Dim > 2)
	{
		Value = getValue(x);
		Derivative = -(x.getx() - Node.getx())
			*(x.getz() - Node.getz()) / (Lx*Lx*Lz*Lz*Value*Value*Value);
		return Derivative;
	}
	else
		return 0;
}
double BasicFunction::get_yzDeriv(Point& x)
{
	if (Dim > 2)
	{
		Value = getValue(x);
		Derivative = -(x.gety() - Node.gety())
			*(x.getz() - Node.getz()) / (Ly*Ly*Lz*Lz*Value*Value*Value);
		return Derivative;
	}
	else
		return 0;
}
double BasicFunction::get_zzDeriv(Point& x)
{
	if (Dim > 2)
	{
		Value = getValue(x);
		Derivative = 1 / (Lz*Lz*Value) 
			- (x.getz() - Node.getz())*(x.getz() - Node.getz()) / (Lz*Lz*Lz*Lz*Value*Value*Value);
		return Derivative;
	}
	else
		return 0;
}
double BasicFunction::get_ScalarDeriv(Point& x)
{
	Derivative = get_xxDeriv(x);
	if (Dim > 1)
		Derivative += get_yyDeriv(x);
	if (Dim > 2)
		Derivative += get_zzDeriv(x);
	return Derivative;
}

double oncefit::Alpha[N_Alpha] = { 1 };
void oncefit::assignParameter(double* Array)
{
	for (int i = 0; i < N_Alpha; i++)
		Alpha[i] = Array[i];
}
void oncefit::getValue(double* fittingValue, Point& x, BasicFunction* bf)
{
	int i = 0;
	if (Dim == 1)
	{
		for (i = 0; i < N; i++)
			fittingValue[0] += Alpha[i] * bf[i].getValue(x);
	}
	if (Dim == 2)
	{
		for (i = 0; i < N; i++)
		{
			fittingValue[0] += Alpha[i] * bf[i].get_yDeriv(x);
			fittingValue[1] += Alpha[i] * (-bf[i].get_xDeriv(x));
		}
	}
	if (Dim == 3)	// B = alpha1 * \nabla \times (\psi1 * r) + alpha2 * \nabla \time \nabla \times (\psi2 * r)
	{
		double Px, Py, Pz, r;
		double* FirstDeriv = new double[Dim]();
		double(*SecondDeriv)[Dim] = new double[Dim][Dim]();
		double* ScalarDeriv = new double();
		Px = x.getx(); Py = x.gety(); Pz = x.getz(); r = x.getradialLength();
		for (i = 0; i < N; i++)
		{
			bf[i].getDerivative(x, FirstDeriv, SecondDeriv, ScalarDeriv);
			fittingValue[0] += Alpha[2 * i] / r*(Pz*FirstDeriv[1] - Py*FirstDeriv[2]);	// 1/r * (z * Partial_y - y * Partial_z)
			fittingValue[0] += Alpha[2 * i + 1] * (2 * FirstDeriv[0] 
				+ (Px*SecondDeriv[0][0] + Py*SecondDeriv[0][1] + Pz*SecondDeriv[0][2]) 
				- Px*(*ScalarDeriv));
			fittingValue[1] += Alpha[2 * i] / r*(Px*FirstDeriv[2] - Pz*FirstDeriv[0]);	// 1/r * (x * Partial_z - z * Partial_x)
			fittingValue[1] += Alpha[2 * i + 1] * (2 * FirstDeriv[1] 
				+ (Px*SecondDeriv[0][1] + Py*SecondDeriv[1][1] + Pz*SecondDeriv[1][2]) 
				- Py*(*ScalarDeriv));
			fittingValue[2] += Alpha[2 * i] / r*(Py*FirstDeriv[0] - Px*FirstDeriv[1]);	// 1/r * (y * Partial_x - x * Partial_y)
			fittingValue[2] += Alpha[2 * i + 1] * (2 * FirstDeriv[2] 
				+ (Px*SecondDeriv[0][2] + Py*SecondDeriv[1][2] + Pz*SecondDeriv[2][2]) 
				- Pz*(*ScalarDeriv));
		}
		delete ScalarDeriv;
		delete[] FirstDeriv;
		delete[] SecondDeriv;
	}
}
