// reading data header file
int read(Point*, double[][Dim], double*, double*, double*, double*);
void CoordinatesTransform(	double[][Dim], double[][Dim], 
							int*, double[][Dim], int*, ifstream& );
double tanhtofit(Point& var);
void magisland(double* value, Point& var);
void RBF_initial(double* value, Point& var, Point* Node, Point* DistParam);