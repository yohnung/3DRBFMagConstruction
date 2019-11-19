// reading data header file
int read(Point*, double[][Dim]);
int read(Point*, double[][Dim], ifstream&);
double tanhtofit(Point& var);
void magisland(double* value, Point& var);
void addRandomPos(int, Point*, Point*);
void ModelField(int, Point*, Point*, Point*, double[][Dim]);
void addRandomError(int, double[][Dim]);
void RBFModelField(int, Point*, Point*, Point*, double[][Dim]);
void SeparatorField1(int, Point*, double[][Dim]);
void SeparatorField2(int, Point*, double[][Dim]);
void write_script_tecplot();