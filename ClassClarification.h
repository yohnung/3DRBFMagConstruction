// Class Clarification header file
class Point
{
public:
	Point();
	Point(double* x);
	Point(Point& x);
	~Point();
	void specify(double *x);
	double getradialLength();
	double getx();
	double gety();
	double getz();
	void assignx(double x);
	void assigny(double y);
	void assignz(double z);
	Point& operator = (Point& b);
	friend Point operator + (Point& a, Point& b);
	friend Point operator - (Point& a, Point& b);
	friend double operator * (Point& a, Point& b);
	friend Point operator * (double a, Point& b);
	friend Point operator / (Point& a, double b);
private:
	double* Coordinate;				// coordinates
	double radialLength;
};

class BasicFunction
{
public:
	BasicFunction();
	BasicFunction(Point& R, Point& D);
	void specify(Point& R, Point& D);
	double getValue(Point& x);
	void getDerivative(	Point& x, double*, 
						double [][Dim], double* );
	double get_xDeriv(Point& x);
	double get_yDeriv(Point& x);	
	double get_zDeriv(Point& x);	
	double get_xxDeriv(Point& x);
	double get_xyDeriv(Point& x);
	double get_xzDeriv(Point& x);
	double get_yyDeriv(Point& x);
	double get_yzDeriv(Point& x);
	double get_zzDeriv(Point& x);
	double get_ScalarDeriv(Point& x);
private:
	Point Node;
	Point DistParameter;
	double Value;
	double Derivative;
};

class oncefit
{
public:
	static void assignParameter(double* Array);
	void getValue(	double* fittingValue, Point& x, 
					BasicFunction* bf);
private:
	static double Alpha[N_Alpha];
};
