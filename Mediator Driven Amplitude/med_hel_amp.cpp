#include <iostream>
#include <complex>
#include <cmath>
#include <time.h>
#include "constants.h"
#include "structures_list.h"
#include "functions_list.h"
#include "Spinor.h"
using namespace std;

int main() {

	double p[5][4];					  //defining arrays for four momentum of the scatterers
	double FMASS[5] = {0}; 	          //defining variables for fermion mass, for ultra relativistic fermions we can consider them massless	
	int NSF[5]={0};				  	  //index for particle(1)/anti-particle(-1)
	int HEL[5]={0};					  //index for helicity of ith particle
	int mtr;						  //choice of mediator particle
	Y M;			  	  			  //defining variables for Feynman amplitude
	double fin_amp=0;				  //final amplitude with contribution from all helicity configurations
	double energy,theta,phi;		  //defining com energy and scattering angles
	double diag_eta[4]= {1,-1,-1,-1};
	Spinor S;	//Creating an object of the Spinor Class to use its functions


	//taking the mediator and fermion details from user
    do {
        cout << "Enter Details for the Lowest Order Tree Amplitude:\n\n";
        cout << "Select Mediator Boson:\n1. Photon\n2. Z Boson\n3. Scalar Boson\nChoice: ";
        cin  >> mtr;
        cout << "\nEnter COM Energy (in exponential form & in units of c): ";
        cin  >> energy;
        for (int i=1;i<5;i++) {
            if (i == 1)
                cout << "\nEnter Details of Incoming Fermions:\n";
            else if(i == 3)
                cout << "Enter Details of Outgoing Fermions:\n";
            cout << "#" << i <<":\n";
            //cout << "Mass: ";
            //cin  >> FMASS[i];
            cout << "Particle(1)/ Anti-Particle(-1): ";
            cin  >> NSF[i];
            cout << "\n";
        }    
    } while(false);

	//generating the four momentum data of the scatterers	
	srand(time(0));				  //initializing random sequence by current value of time
	theta  = getrand(180)*d2r;
	phi    = getrand(360)*d2r;

	//initializing the momenta of the incoming particles
	p[1][0] = energy;	p[2][0] = energy;
	p[1][1] = 0;		p[2][1] = 0;
	p[1][2] = 0;		p[2][2] = 0;
	p[1][3] = energy;	p[2][3] = -energy;
	
	//initializing the momenta of the outgoing particles			
	p[3][0] = energy;						p[4][0] = energy;
	p[3][1] = energy*sin(theta)*cos(phi);	p[4][1] = energy*sin(pi-theta)*cos(phi + pi); 
	p[3][2] = energy*sin(theta)*sin(phi);	p[4][2] = energy*sin(pi-theta)*sin(phi + pi);
	p[3][3] = energy*cos(theta);			p[4][3] = energy*cos(pi-theta);

	for (int i=0;i<16;i++) {

		//generating the binary form of 0-15 for helicity configs
		int dec = i;  
		int hctr = 4; //helicity index for each particle

		while(dec > 0) {
			HEL[hctr] = dec%2;
			dec = dec/2;
			hctr--;
		}
		double pp[4]={0};				  //defining the propagator four momenta
		Y spp;							  //square of propagator four momenta
		int ictr = 0; //inflowing fermion counter
		int octr = 0; //outflowing fermion counter
		X IF[2];	//declaring the inflowing fermion wavefunction
		X OF[2];	//declaring the outflowing fermion wavefunction
		for (int j=1;j<5;j++) {
			//Generating the wavefunction for the fermion #i

			if (((j <= 2) && (NSF[j] == 1)) || ((j > 2) && (NSF[j] == -1))) {
				IF[ictr] = S.IFWFN(p[j],FMASS[j],HEL[j],NSF[j]);
				ictr++;
				if (ictr == 1) {
					for (int k=0;k<4;k++)
						pp[k] = pp[k] + p[j][k];
				}

			}
			else {
				OF[octr] = S.OFWFN(p[j],FMASS[j],HEL[j],NSF[j]);
				octr++;
				if (octr == 1) {
					for (int k=0;k<4;k++)
						pp[k] = pp[k] + p[j][k];
				}
			}
		}

		//initialising the square of propagator four momenta
		for (int j=0;j<4;j++)
			spp.z = spp.z + diag_eta[j]*pp[j]*pp[j];

		//Calculating the helicity amplitude for this particular configuration
		switch (mtr)
		{
		case 1: M = S.P_AMP(IF,OF,spp);
			break;
		case 2: M = S.Z_AMP(IF,OF,spp);
			break;
		case 3: M = S.S_AMP(IF,OF,spp);
			break;
		}

		cout << "Feynman Amplitude of the process" << "#1["<<2*HEL[1] - 1<<"]" << " + #2["<<2*HEL[2] - 1<<"] ---> " <<"#3["<<2*HEL[3] - 1<<"]" << " + #4["<<2*HEL[4] - 1<<"] is" << M.z << "\n";
		cout << "Norm of Feynman amplitude of this process is " << norm(M.z) << "\n\n";
		fin_amp = fin_amp + norm(M.z);
	}
		
	cout << "Total amplitude is : " << fin_amp << endl;
	return 0;
}