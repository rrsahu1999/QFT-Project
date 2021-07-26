#include <iostream>
#include <fstream>
#include <complex>
#include <cmath>
#include<time.h>
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
	double eta[4][4]= {{1,0,0,0},{0,-1,0,0},{0,0,-1,0},{0,0,0,-1}};
	double e = 0.001;		//precision parameter

	double getrand(int a) {

		double dec,rvalue;
		dec = (rand() % 100)*0.01;
		rvalue = (rand() % a) + dec;
		return rvalue;
	}
	
	X IXXXXX(double p[],double FMASS,int NHEL,int NSF) {
		if ( NSF == -1) {		//if inflowing particle is an anti-fermion, then use crossing symmetry for massless particles
			NHEL = -NHEL;		//v-/+(p) = u+/-(p)
		}		
		//computing the angle parameters theta & phi from given four momentum
		theta = acos(p[3]/p[0]);		//gives theta from 0 to 180
		
		if ( (abs(theta) < e) or (abs(theta-pi) < e))
			phi = 0;
		else {
			phi = acos(p[1]/(p[0]*sin(theta)));    //gives phi from 0 to 180
			
			if ( (p[2]/(p[0]*sin(theta))) < 0 )
				phi = 2*pi - phi;		//gives phi from 0 to 360 as required
		}
					
		//computing the square ket
		squareket[0] = sqrt(2*p[0])*-1*sin(theta/2)*exp(-I*phi);
		squareket[1] = sqrt(2*p[0])*cos(theta/2);

		//computing the angle ket
		angleket[0]  = squareket[1];
		angleket[1]  = -1.0*conj(squareket[0]);

		X FI;		////Creating the inflowing fermion wavefunction
		if (NHEL == 1) {
			FI.wfn[0] = 0; FI.wfn[1] = 0;			
			FI.wfn[2] = angleket[0];
			FI.wfn[3] = angleket[1];

		}
		else	{
			FI.wfn[0] = squareket[0];
			FI.wfn[1] = squareket[1];
			FI.wfn[2] = 0; FI.wfn[3] = 0;
		}
		return FI;
					
	}

	X OXXXXX(double p[],double FMASS,int NHEL,int NSF) {
		
		if ( NSF == -1) {		//if outflowing particle is an anti-fermion, then use crossing symmetry for massless particles
			NHEL = -NHEL;		//vbar+/-(p) = ubar-/+(p)
		}		
		//computing the angle parameters theta & phi from given four momentum
		theta = acos(p[3]/p[0]);		//gives theta from 0 to 180
		
		if ( (abs(theta) < e) or (abs(theta-pi) < e))
			phi = 0;
		else {
			phi = acos(p[1]/(p[0]*sin(theta)));    //gives phi from 0 to 180
			
			if ( (p[2]/(p[0]*sin(theta))) < 0 )
				phi = 2*pi - phi;		//gives phi from 0 to 360 as required
		}

		//computing the square bra
		squarebra[0] = sqrt(2*p[0])*cos(theta/2);
		squarebra[1] = sqrt(2*p[0])*sin(theta/2)*exp(-I*phi);

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

	Y FSIXXX(X FI1, X FO2,X FO3,X FI4,double p1[],double p2[]) {

		double ppp[4] = {0}, s12 = 0;				//variable for photon propagator four momenta and its square

		for (int i=0;i<4;i++)
			ppp[i] = p1[i] + p2[i];

		//calculating the square of photon propagator four momentum to be used later
		for (int i=0;i<4;i++)
			s12 = s12 + eta[i][i]*ppp[i]*ppp[i];

		Y M1,M2;		//defining complex variable to store the matrix elements	
		Y M;			//defining complex variable to store helicity amplitude			
		M.z = 0;		//initializing to zero
			
		for (int k=0;k<4;k++) {
			M1.z = 0; M2.z = 0;	//initializing to zero	
			
			for (int i=0;i<4;i++) {
				for (int j=0;j<4;j++) {
					M1.z = M1.z + FO2.wfn[i]*g[k][i][j]*FI1.wfn[j];
					M2.z = M2.z + FO3.wfn[i]*eta[k][k]*g[k][i][j]*FI4.wfn[j];
				}	
			}
			M.z = M.z + M1.z*M2.z;
		}

		M.z = M.z*I/s12;
		
		return M;
	}
	
};

int main() {

	double p[5][4]; //defining arrays for four momentum of the scatterers
	double FMASS; int NHEL,NSF;	          //defining variables for fermion mass, helicity(+1/-1) and index for particle(1)/anti-particle(-1)
	FMASS = 0;				  //for ultra relativistic fermions we can consider them massless	
	double GC = 1;		  		  //defining the coupling constant of electromagnetic interaction
	Y M;			  	  //defining variables for Feynman amplitude
	double energy,theta,phi;		  //defining com energy and scattering angles
	int hc[4][5] = {{0,-1,1,1,-1},
			 {0,1,-1,-1,1},
			 {0,1,-1,1,-1},
			 {0,-1,1,-1,1}};	  //initializing the helicity configurations
	double ans = 0;		 
	Spinor S;	//Creating an object of the Spinor Class to use its functions
	
	//generating the four momentum data of the scatterers
	energy = 1.0e-18;  			  //setting com energy = 2 GeV in units of c	
	srand(time(0));				  //initializing randome sequence by current value of time
	theta  = S.getrand(180)*d2r;
	phi    = S.getrand(360)*d2r;

		//initializing the initial momenta of the electron and positron
	p[1][0] = energy;	p[2][0] = energy;
	p[1][1] = 0;		p[2][1] = 0;
	p[1][2] = 0;		p[2][2] = 0;
	p[1][3] = energy;	p[2][3] = -energy;
	
		//initializing the final momenta of the muon and anti-muon			
	p[3][0] = energy;			p[4][0] = energy;
	p[3][1] = energy*sin(theta)*cos(phi);	p[4][1] = energy*sin(pi-theta)*cos(phi + pi); 
	p[3][2] = energy*sin(theta)*sin(phi);	p[4][2] = energy*sin(pi-theta)*sin(phi + pi);
	p[3][3] = energy*cos(theta);		p[4][3] = energy*cos(pi-theta);

	ofstream A;				  //Object of ifstream class to write the four momentum data of the scatterers to verify later
	A.open("four_momentum_data.dat", ios::out);	
	for (int i=1;i<5;i++) {
		for (int j=0;j<4;j++)
			A << p[i][j] << " ";
		A << endl;
	}
	A.close();

	for (int i=0;i<4;i++) {
				//Obtaining the wave function for the inflowing particles 1 & 4

		// Helicity configuration (i)
		//	Electron	Positron	Muon		Anti-muon
		//	p1,hc[i][1]	p2,hc[i][2]	p3,hc[i][3]	p4,hc[i][4]
				
		//Generating the wavefunction for the incoming electron of momentum p1
		X FI1;
		NHEL = hc[i][1]; 
		NSF  = 1;	
		FI1 = S.IXXXXX(p[1],FMASS,NHEL,NSF);

		//Generating the wavefunction for the incoming positron of momentum p2
		X FO2;
		NHEL = hc[i][2]; 
		NSF  = -1;	
		FO2 = S.OXXXXX(p[2],FMASS,NHEL,NSF);
		
		//Generating the wavefunction for the outgoing muon of momentum p3
		X FO3;
		NHEL = hc[i][3]; 
		NSF  = 1;	
		FO3 = S.OXXXXX(p[3],FMASS,NHEL,NSF);

		//Generating the wavefunction for the outgoing anti-muon of momentum p2
		X FI4;
		NHEL = hc[i][4]; 
		NSF  = -1;	
		FI4 = S.IXXXXX(p[4],FMASS,NHEL,NSF);

		//Calculating the helicity amplitude for this particular configuration

		M = S.FSIXXX(FI1,FO2,FO3,FI4,p[1],p[2]);
		
		cout << "Feynman Amplitude of the process" << "e-["<<hc[i][1]<<"]" << " + e+["<<hc[i][2]<<"] --> " <<"mu-["<<hc[i][3]<<"]" << " + mu+["<<hc[i][4]<<"] is" << M.z << endl;
		cout << "Norm of Feynman amplitude of this process is " << norm(M.z) << endl;
		ans = ans + norm(M.z);
	}
		
	cout << "Total amplitude is : " << ans << endl;
	return 0;
}