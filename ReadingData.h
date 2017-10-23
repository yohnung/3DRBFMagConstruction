// reading data header file
int read(Point*, double[][Dim]);
int read(Point*, double[][Dim], ifstream&);
double tanhtofit(Point& var);
void magisland(double* value, Point& var);
void RBFModelField(int, Point*, Point*, Point*, double[][Dim]);