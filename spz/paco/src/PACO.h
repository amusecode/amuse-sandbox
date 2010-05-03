
void PatternConstruct(double x[],double y[], double xd[], double yd[],int N,char * pattern,int *size);
void PartialPatternConstruct(double x0, double y0, 
			     double vx0, double vy0,
			     double x1, double y1, 
			     int *sign0, int *sign1);

float AutoCorrelate(char * pattern,float * data,int size,float Y=-1,int method = 1);
float * Delta(float * delta,int Dsize,char * pattern,int psize,int lpsize);
float PSratio(float * delta, int size); 

void Savedata(float * data,int size,char * name);
void SavePattern(char * pattern,int size,char * name);
