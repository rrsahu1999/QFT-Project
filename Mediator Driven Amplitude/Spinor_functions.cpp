#include <iostream>
#include <complex>
#include <cmath>
#include "constants.h"
#include "structures_list.h"
#include "Spinor.h"
using namespace std;

//function to compute inflowing fermion wavefunction
X Spinor::IFWFN(double p[],double FMASS,int HEL,int NSF) {

	HEL = 2*HEL - 1;        // mapping 0,1 ---> -1,+1
	if ( NSF == -1) {		//if inflowing particle is an anti-fermion, then use crossing symmetry for massless particles
		HEL = -HEL;		//v-/+(p) = u+/-(p)
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

	X IF;		////Creating the inflowing fermion wavefunction
	if (HEL == 1) {
		IF.wfn[0] = 0; IF.wfn[1] = 0;			
		IF.wfn[2] = angleket[0];
		IF.wfn[3] = angleket[1];

	}
	else	{
		IF.wfn[0] = squareket[0];
		IF.wfn[1] = squareket[1];
		IF.wfn[2] = 0; IF.wfn[3] = 0;
	}
	return IF;
				
}

//function to compute outflowing fermion wavefunction
X Spinor::OFWFN(double p[],double FMASS,int HEL,int NSF) {
		
	HEL = 2*HEL - 1;        // mapping 0,1 ---> -1,+1        
	if ( NSF == -1) {		//if outflowing particle is an anti-fermion, then use crossing symmetry for massless particles
		HEL = -HEL;		//vbar+/-(p) = ubar-/+(p)
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

	X OF;		////Creating the outflowing fermion wavefunction
	if (HEL == 1) {
		OF.wfn[0] = squarebra[0];
		OF.wfn[1] = squarebra[1];
		OF.wfn[2] = 0; OF.wfn[3] = 0;
	}
	else	{
		OF.wfn[0] = 0; OF.wfn[1] = 0;			
		OF.wfn[2] = anglebra[0];
		OF.wfn[3] = anglebra[1];
	}
	return OF;
				
}

//function to compute amplitude for Photon Propagator
Y Spinor::P_AMP(X IF[2],X OF[2],Y spp) {

	double GC_ph = 1;		  		  //defining the coupling constant of electromagnetic interaction

	Y M1,M2;		//defining complex variable to store the matrix elements	
	Y M;			//defining complex variable to store helicity amplitude

	for (int k=0;k<4;k++) {			
		for (int i=0;i<4;i++) {
			for (int j=0;j<4;j++) {
				M1.z = M1.z + OF[0].wfn[i]*gamma[k][i][j]*IF[0].wfn[j];
				M2.z = M2.z + OF[1].wfn[i]*diag_eta[k]*gamma[k][i][j]*IF[1].wfn[j];
			}	
		}
		M.z = M.z + M1.z*M2.z;
	}
	M.z = -1.0*GC_ph*M.z*I/spp.z;
	
	return M;
}

//function to compute amplitude for Z-Boson Propagator
Y Spinor::Z_AMP(X IF[2],X OF[2],Y spp) {

	
    double GC_zb = 1;		  		  //defining the coupling constant of neutral weak interaction
	double gz = 2.4952e9*eV/c;		  //decay period of z boson (in units of c)
	double mz = 91.9e9*eV/c;	      //mass of z boson (in units of c)
	double ga=0.2,gb=0.8;			  //weak interaction parameters
	complex<double> gamma_z[4][4][4]; //matrix for z boson vertex

	for (int i=0;i<4;i++) {
		for (int j=0;j<4;j++) {
			for (int k=0;k<4;k++)
				gamma_z[k][i][j] = gamma[k][i][j]*(ga*diag_id[j] + gb*diag_g5[j]);
		}
	}
	spp.z = spp.z - mz*mz + I*gz*mz; 

	Y M1,M2;		//defining complex variable to store the matrix elements	
	Y M;			//defining complex variable to store helicity amplitude

	for (int k=0;k<4;k++) {			
		for (int i=0;i<4;i++) {
			for (int j=0;j<4;j++) {
				M1.z = M1.z + OF[0].wfn[i]*gamma_z[k][i][j]*IF[0].wfn[j];
				M2.z = M2.z + OF[1].wfn[i]*diag_eta[k]*gamma_z[k][i][j]*IF[1].wfn[j];
			}	
		}
		M.z = M.z + M1.z*M2.z;
	}

	M.z = -1.0*GC_zb*M.z*I/spp.z;
	
	return M;
}

//function to compute amplitude for Scalar Boson Propagator
Y Spinor::S_AMP(X IF[2],X OF[2],Y spp) {

    double GC_sb = 1;		//defining the coupling constant of interaction via scalar boson

	Y M1,M2;		//defining complex variable to store the matrix elements	
	Y M;			//defining complex variable to store helicity amplitude
	
	for (int i=0;i<4;i++) {
				M1.z = M1.z + OF[0].wfn[i]*IF[0].wfn[i];
				M2.z = M2.z + OF[1].wfn[i]*IF[1].wfn[i];
	}
	M.z = M1.z*M2.z;
	M.z = GC_sb*M.z*I/spp.z;
	
	return M;
}