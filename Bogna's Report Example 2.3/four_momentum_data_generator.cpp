#include<time.h>
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const double pi  = acos(-1.0);
const double d2r = pi/180;        	  //defining constant for converting degrees to radians
const double c   = 3e8;			  //defining constant for speed of light


double getrand(int a,int b) {

	double dec,rvalue;
	dec = (rand() % 100)*0.01;
	rvalue = (rand() % b) + 1 + dec;
	return rvalue;
}

int main() {
	double p[5][4] = {0},energy,theta,phi;    //defining arrays for four momentum, energy, angle data of scatterers
	srand(time(0));				  //initializing randome sequence by current value of time
	ofstream A;				  //Object of ofstream class to write a sample four momentum data of scatterers
	A.open("four_momentum_data_.dat", ios::out|ios::trunc);

	energy = getrand(1,100)*1.6e-13;  //any energy from 1 to 100MeV
	for (int i=0;i<2;i++) {		
		theta  = getrand(1,180)*d2r;
		phi    = getrand(1,180)*d2r;	
	
			p[i+1][0] = -1;			  p[i+3][0] = -1;
			p[i+1][1] = sin(theta)*cos(phi);  p[i+3][1] = sin(pi-theta)*cos(phi + pi); 
			p[i+1][2] = sin(theta)*sin(phi);  p[i+3][2] = sin(pi-theta)*sin(phi + pi);
			p[i+1][3] = cos(theta);		  p[i+3][3] = cos(pi-theta);
			
			for (int j=0;j<4;j++) {
				p[i+1][j] = p[i+1][j]*energy/c;
				p[i+3][j] = p[i+3][j]*energy/c;
			}
	}
	for (int i=1;i<5;i++) {	
		for (int j=0;j<4;j++)
			A << p[i][j] << " ";	
		A << endl;	
	}
	A.close();
	return 0;
}
