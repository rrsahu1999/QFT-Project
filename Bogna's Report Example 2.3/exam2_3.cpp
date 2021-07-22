#include <iostream>
#include <fstream>
#include <complex>
#include <cmath>
using namespace std;

const double pi  = acos(-1.0);
const double d2r = pi/180;        //defining constant for converting degrees to radians
const double c   = 3e8;			  //defining constant for speed of light
complex<double> I(0,1);

struct  X{
	complex<double> wfn[4];                      //defining a structure X with a complex 1-d array with 4 elements for the wave function
};

struct  Y{					     //defining a structure Y to receive complex numbers from functions
	complex<double> z;
};

class Spinor {

	public:
	complex<double> squareket[2],squarebra[2],angleket[2],anglebra[2]; //defining the different types of Weyl spinors
	double theta,phi;
	//defining the gamma matrices in Weyl representation;
	//         row 1	      //        row 2 		//        row 3 		 //       row 4 
	complex<double> g[4][4][4] =
        {{{{0,0},{0,0},{1,0},{0,0}},{{0,0},{0,0},{0,0},{1,0}},{{1,0},{0,0},{0,0},{0,0}}, {{0,0},{1,0},{0,0},{0,0}}},   /*gamma 0*/
        {{{0,0},{0,0},{0,0},{1,0}},{{0,0},{0,0},{1,0},{0,0}},{{0,0},{-1,0},{0,0},{0,0}},{{-1,0},{0,0},{0,0},{0,0}}},   /*gamma 1*/
        {{{0,0},{0,0},{0,0},{0,-1}},{{0,0},{0,0},{0,1},{0,0}},{{0,0},{0,1},{0,0},{0,0}},{{0,-1},{0,0},{0,0},{0,0}}},   /*gamma 2*/
        {{{0,0},{0,0},{1,0},{0,0}},{{0,0},{0,0},{0,0},{-1,0}},{{-1,0},{0,0},{0,0},{0,0}},{{0,0},{1,0},{0,0},{0,0}}}};   /*gamma 3*/
	double eta[4][4]= {{-1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};

	X IXXXXX(double p[],double FMASS,int NHEL,int NSF) {
		if ( NSF == 1) {		//if inflowing particle is a fermion, then use crossing symmetry for massless particles
			NHEL = -NHEL;		//u+/-(p) = v-/+(p)
		}		
		//computing the angle parameters theta & phi from given four momentum
		theta = acos(-p[3]/p[0]);		//gives theta from 0 to 180
		phi = acos(-p[1]/(p[0]*sin(theta)));    //gives phi from 0 to 180
		
		if ( (-p[2]/(p[0]*sin(theta))) < 0 )
			phi = 2*pi - phi;		//gives phi from 0 to 360 as required
					
		//computing the square ket
		squareket[0] = sqrt(-2*p[0])*-1*sin(theta/2)*exp(-I*phi);
		squareket[1] = sqrt(-2*p[0])*cos(theta/2);

		//computing the angle ket
		angleket[0]  = squareket[1];
		angleket[1]  = -1.0*conj(squareket[0]);

		X FI;		////Creating the inflowing fermion wavefunction
		if (NHEL == 1) {
			FI.wfn[0] = squareket[0];
			FI.wfn[1] = squareket[1];
			FI.wfn[2] = 0; FI.wfn[3] = 0;
		}
		else	{
			FI.wfn[0] = 0; FI.wfn[1] = 0;			
			FI.wfn[2] = angleket[0];
			FI.wfn[3] = angleket[1];
		}
		return FI;
					
	}

	X OXXXXX(double p[],double FMASS,int NHEL,int NSF) {
		if ( NSF == -1) {		//if outflowing particle is an anti-fermion, then use crossing symmetry for massless particles
			NHEL = -NHEL;		//vbar+/-(p) = ubar-/+(p)
		}		
		//computing the angle parameters theta & phi from given four momentum
		theta = acos(-p[3]/p[0]);		//gives theta from 0 to 180
		phi = acos(-p[1]/(p[0]*sin(theta)));    //gives phi from 0 to 180
		
		if ( (-p[2]/(p[0]*sin(theta))) < 0 )
			phi = 2*pi - phi;		//gives phi from 0 to 360 as required
		
		//computing the square bra
		squarebra[0] = sqrt(-2*p[0])*cos(theta/2);
		squarebra[1] = sqrt(-2*p[0])*sin(theta/2)*exp(-I*phi);

		//computing the angle bra
		anglebra[0]  = -1.0*conj(squarebra[1]);
		anglebra[1]  = squarebra[0];

		X FO;		////Creating the outflowing fermion wavefunction
		if (NHEL == 1) {
			FO.wfn[0] = squarebra[0];
			FO.wfn[1] = squarebra[1];
			FO.wfn[2] = 0; FO.wfn[3] = 0;
		}
		else	{
			FO.wfn[0] = 0; FO.wfn[1] = 0;			
			FO.wfn[2] = anglebra[0];
			FO.wfn[3] = anglebra[1];
		}
		return FO;
					
	}

	Y FSIXXX(X FI,X FO,double SC, double GC,double sbp[],double ifp[]) {
		// sbp[] defines the scalar boson four momentum
		// ifp[] defines the inflowing fermion four momentum

		complex<double> fp[4][4] = {0};		//defining the fermion propagator
		
		double fpp[4] = {0};			//defining the fermion propagator four momentum
		double q = 0;				//variable for square of fermion propagator four momentum
		for (int i=0;i<4;i++)
			fpp[i] = sbp[i] + ifp[i];

		//calculating the square of fermion propagator four momentum to be used later
		for (int i=0;i<4;i++) {
			for (int j=0;j<4;j++)
				q = q + eta[i][j]*fpp[i]*fpp[j];
		}

		//computing the entries of the fermion propagator
		for (int i=0;i<4;i++) {
			for (int j=0;j<4;j++) {
				for (int k=0;k<4;k++)
					fp[i][j] = fp[i][j] + fpp[k]*g[k][i][j];
				fp[i][j] = -1.0*fp[i][j]/q;
			}
		}

		Y M;		//defining complex variable to store the Feynman amplitude	
		M.z = 0; 	//initializing to zero				
		for (int i=0;i<4;i++) {
			for (int j=0;j<4;j++)
				M.z = M.z + FO.wfn[i]*fp[i][j]*FI.wfn[j];	
		}
		M.z = M.z*(I*GC)*(I*GC);
		
		return M;
	}
	
};

int main() {

	double p[5][4]; //defining arrays for four momentum of the scatterers
	double FMASS; int NHEL,NSF;	          //defining variables for fermion mass, helicity(+1/-1) and index for particle(1)/anti-particle(-1)
	FMASS = 0;				  //this particular example considers massless fermions	
	double SC = 1;				  //defining the wavefunction for the scalar boson which is unity
	double GC = 1;		  		  //defining the coupling constant of fermion and scalar boson
	Y M_a,M_b;		  		  //defining variables for Feynman amplitude of each Feynman diagram
	complex<double> M;			  //final answer of total Feynman amplitude

	ifstream B;				  //Object of ifstream class to read the four momentum data of the scatterers
	B.open("four_momentum_data_.dat", ios::in);	
	for (int i=1;i<5;i++) {
		for (int j=0;j<4;j++)
			B >> p[i][j];
	}
	B.close();

//	Computing Feynman amplitude for h2=-h4=+1/2 case

	//	Computing 1st part of amplitude with 1&2 paired

	Spinor S;	//Creating an object of the Spinor Class to use its functions

			//Obtaining the wave function for the inflowing particle 2, which is a anti-fermion with +1 helicity
	X FI;		//Creating the inflowing wavefunction
	NHEL = 1; 
	NSF  = -1;	
	FI = S.IXXXXX(p[2],FMASS,NHEL,NSF);

			//Obtaining the wave function for the outflowing particle 4, which is a fermion with -1 helicity

	X FO;		//Creating the outflowing wavefunction
	NHEL = -1; 
	NSF  = 1;	
	FO = S.OXXXXX(p[4],FMASS,NHEL,NSF);

	X FP;		//defining the fermion propagator times the inflowing fermion wavefunction
	

	M_a = S.FSIXXX(FI,FO,SC,GC,p[1],p[2]);
	M_b = S.FSIXXX(FI,FO,SC,GC,p[3],p[2]);
	M = M_a.z + M_b.z;
	
	cout << "Feynman Amplitude of A4(12+34-) process is " << M << endl;
	cout << "Norm of Feynman amplitude of this process is " << norm(M) << endl;

	return 0;
}
