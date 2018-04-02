// reading data header file
int read(Point*, double[][Dim]);
int read(Point*, double[][Dim], ifstream&);
double tanhtofit(Point& var);
void magisland(double* value, Point& var);
void ModelField(int, Point*, Point*, Point*, double[][Dim]);
void RBFModelField(int, Point*, Point*, Point*, double[][Dim]);
void SeparatorField(int, Point*, double[][Dim]);
void write_script_tecplot();