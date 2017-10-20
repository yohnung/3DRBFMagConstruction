// commom procedure header file
void constructNodesWeb(Point*, Point*, double*, double*);
void specifyAY(double[][N_Alpha], double*, double[][Dim], Point*, int, BasicFunction*);
void generalSolver(double*, Point*, double[][Dim], int, BasicFunction*);
double LinearLUSolver(double*, double[][N_Alpha], double*, int);
void LinearQRSolver(double*, double[][N_Alpha], double*, int, double*);
void LinearSVDSolver(double*, double[][N_Alpha], double*, int, double*, double);
void write(double[][Dim], int, ofstream&);
void write(double *, int, ofstream&);
void write(Point*, int, ofstream&);
void write(double*, double[][N_Alpha], ofstream&);
