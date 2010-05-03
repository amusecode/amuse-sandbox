
#include <fstream> 
#include <iostream>
#include <math.h> 
#include "PACO.h" 

using namespace std; 

// PACO auxiliary functions 
void PatternConstruct(double x[],double y[], 
		      double xd[], double yd[],int N,char * pattern, int *size)
{
  int n=0;
  int sign0, sign1;
  for(int i=0; i<N-1; i++) {

    PartialPatternConstruct(x[i], y[i], xd[i], yd[i],
			    x[i+1],  y[i+1], 
			    &sign0, &sign1);
    if (sign0>=0) {
      pattern[n++] = sign0 ? 'Y': 'X';
    }
    //    cout <<"i=" << i<< "pattern[i]="<< pattern[i] << endl;
  }

  if (sign1>=0)
    pattern[n++] = sign0 ? 'Y': 'X';
  *size = n; 
}

void PartialPatternConstruct(double x0, double y0, 
			     double vx0, double vy0,
			     double x1, double y1, 
			     int *sign0, int *sign1) {

  if(x0*x1<0 && y0*y1<0) {
    double tx=fabs(x0/vx0);
    double ty=fabs(y0/vy0);
    
    if(tx<=ty) { 
      *sign0 = 0; 
      *sign1 = 1; 
    }
    else { 
      *sign0 = 1; 
      *sign1 = 0; 
    }
  }    
  else if(x0*x1<0)  {
    *sign0 = 0; 
    *sign1 = -1; 
  }
  else if(y0*y1<0) {
    *sign0 = 1;
    *sign1 = -1;
  }
  else {
    *sign0 = -1;
    *sign1 = -1;
  }
}

// Autocorrelation fonction for the pattern P constructed from the signal 
// 
float AutoCorrelate(char * P, float * data, int size, float Upsilon, int method)
{
  int i, lp; for( i=0;i<size;i++) data[i] = 0;

	// method = 1 by default, similar to autocorrelation product for binary (0,1) data
	//
	if(method)
		for(int i=0;i<size;i++)
		{
			for(int j=0;j<size;j++) { data[i] += (P[j] == P[(i+j)%size]);}
			data[i] = data[i] / size;
		}
	else for(int i=0;i<size;i++)
		 {
		 	for(int j=0;j<size-i;j++)
			     {data[i] += (P[j] == P[j+i]); }
	           	 data[i] = data[i] / (size-i);
		 }		

	// Renormalise autocorrelation data to 1 or 0 (binary) using threshold Upsilon 
	//
	if(Upsilon >= 0){
		for(int i=0;i<size;i++)
			data[i] = (data[i] >= Upsilon);
	          }

	// Identify length of the first harmonic giving non-zero corelation for this threshold
	// Skip i = 0 (strict autocorrelation point always evaluates to 1) 
	//
	i = 1 ; while( data[i] == 0 && i < size ) { i++ ; } 
	if( i == 1 || i == size ) { lp = 0 ; } else { lp = i ; } 

return lp;
}

// Cumbing fonction delta which identifies sections of the signal that contain the repeating 
// pattern p of length lp. 
//
float * Delta(float * delta,int Dsize,char * pattern,int N,int lp)
{	for(int i = 0;i<N-2*lp;i++)
		delta[i] = (pattern[i]==pattern[lp+i]);
return delta;
}

// Pattern-to-Signal ratio which is an estimate of the fraction of the total signal during which 
// the pattern repeats itself. 
// 
float PSratio(float * delta,int size) 
{
	float sum = 0;
	for(int i=0;i<size;i++) { sum += delta[i]; } 
	return sum/size;
}

// Two functions to save data to disc: applied to autocorrelation function and cumbing function delta.
//
void Savedata(float * data,int size,char * name)
{
	ofstream fout(name,ios::trunc);
	for(int i = 0;i<size;i++)
		fout << i << "    " << data[i] << endl;
	fout.close();
}

void SavePattern(char * pattern,int size,char * name)
{
	ofstream fout(name,ios::trunc);
	for(int i = 0;i<size;i++)
		fout << pattern[i];
	fout.close();
}

#ifndef TOOLBOX

/* PACO - a C++ code to identify regular and semi-regular orbits using
 * pattern recognition.  See the README file for pointers / sample
 * data.
 *
 * Code written in April 2010 by Sotiris Chatzopoulos & Christian Boily, 
 * Academy of Athens, Greece, and Observatoire astronomique, Strasbourg
 *
 * Main reference: N. Faber, F. Flitti, C.M. Boily et al. 2010, MNRAS,
 * submitted
 *
 * Contact: christian.boily@astro.unistra.fr  
 */ 

#include <fstream> 
#include <iostream>
#include <math.h> 

#include <stdlib.h>


using namespace std;

///////////////// Main source code ////////////////
//
int main(  int argc, char *argv[])
{ int N = 0; 

  // PART 1: Data input
  // Check length of signal: import data from file signal.dat
  //
  ifstream inputf ;  inputf.open("signal.dat") ; 

  if( inputf.is_open() ) {   float t, x,y, vx, vy; 
    while ( !inputf.eof() ) {  inputf >> t >>  x >> y >> vx >> vy ; N++; }
    inputf.close() ; 
    }
 else { cout << "Error on opening signal.dat\n" ; exit(-1) ; } 
    cout << " Total number of input data = " << --N << endl; 

  // Re-start input data but this time with properly sized arrays - 
  //
  double t[N], x[N], y[N], vx[N], vy[N] ; char P[N] ; 

  inputf.open("signal.dat") ; int count = 0; 
   while ( !inputf.eof() ) { inputf >> t[count] >> x[count] >> y[count] >> vx[count] >> vy[count] ; count++; }
  inputf.close() ; 

  // PART 2: Orbit crossings 
  // Construct binary (XY or 01) pattern P from input signal - hereafter N is the length of P 
  //
  PatternConstruct( x,y, vx, vy, N, P, &N) ; 

  for (int i=0;i<N; i++)
    cout << P[i];
  cout << endl;


  // PART 3 - Pattern recognition & PS ratio
  // 
        // Define threshold for detection of pattern - optional command-line parameter
	//
        float Upsilon = 0.95 ; if( argc > 1 ){ Upsilon = fmin( 0.99, atof( argv[1] ) ) ; } 

	//autocorraletion data will be stored in array autocdata of length N 
	//
	float * autocdata = new float[N];

	//Define lp, the size of the unit pattern p
        //
	int lp;

	// Compute the auto-correletion fnc and the lenght lp of the unit pattern p
	//
        lp = AutoCorrelate(P, autocdata, N, Upsilon, 1) ; cout << " Length lp = " << lp << endl ;

	// Save autocorraletion data to  file ACdata.dat
	//
	Savedata( autocdata, N, "ACdata.dat" );

	// Array delta contains all harmonics of time-phase "tau", i.e. index "i" (cf. Eq 5 & 9 of Faber at al.) 
	//
	float * delta = new float[N-lp];

	// Calculate delta function and save it to a file for diagnostics
	//
	delta =  Delta(delta,N -lp,P,N,lp); Savedata( delta, N -lp,"Delta.dat");

	//calculate PS ratio, the fraction of time a signal is in the pattern identified (Eq. 10 of Faber et al.)
	//
	float ps = PSratio(delta,N-2*lp); cout << " PS ratio is " <<ps<< " for threshold Upsilon = " <<Upsilon<< endl; 

	// Finally, identify the unit pattern p of length lp : recall that delta(i) = 0 or 1 (integer values)
	//
	int check = 0, i = 0 ; while ( check < lp && i < N - lp ) { check++; check *= delta[i]; i++; } 
	cout << " Pattern found: " ; 
	if( i < N-lp && i >= lp ) { check = i;  for ( i=i-lp; i< check; i++ ){ cout << P[i] ; } } 
	else { cout << " Wrong pattern length .. exit \n" ; exit(-1) ; }  cout << "\n" ;
	SavePattern(P+check-lp,lp,"PACO.dat") ;

return 0;
}


#endif //TOOLBOX
