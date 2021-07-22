#include <iostream>
#include <fstream>
#include <complex>
#include <cmath>
using namespace std;

const double pi  = acos(-1.0);
const double d2r = pi/180;        //defining constant for converting degrees to radians
const double c   = 3e8;			  //defining constant for speed of light
complex<double> I(0,1);

int main() {

	double p[5][4];
	double p12[4] = {0}, p13[4] = {0}, p14[4] = {0};
	double s12 = 0, s13 = 0, s14 = 0;
	double eta[4][4]= {{1,0,0,0},{0,-1,0,0},{0,0,-1,0},{0,0,0,-1}};
	int hc[4][5] = {{0,-1,1,1,-1},
		 {0,1,-1,-1,1},
		 {0,1,-1,1,-1},
		 {0,-1,1,-1,1}};	  //initializing the helicity configurations
	double ans[4] = {0}, tans = 0;

	ifstream B;				  //Object of ofstream class to write covariant four momentum of the scatterers to a file
	B.open("four_momentum_data.dat", ios::in);	
	for (int i=1;i<5;i++) {
		for (int j=0;j<4;j++)
			B >> p[i][j];
	}
	B.close();
	
	for (int i=0;i<4;i++) {
		p12[i] = p[1][i] + p[2][i];
		p13[i] = p[1][i] + p[3][i];
		p14[i] = p[1][i] + p[4][i];
	}
		
	for (int i=0;i<4;i++) {
			s12 = s12 + eta[i][i]*p12[i]*p12[i];
			s13 = s13 + eta[i][i]*p13[i]*p13[i];
			s14 = s14 + eta[i][i]*p14[i]*p14[i];
	}

	ans[0] = 4*s13*s13/(s12*s12);	ans[1] = ans[0];
	ans[2] = 4*s14*s14/(s12*s12);	ans[3] = ans[2];

	cout << "Norm of Feynman Amplitude for the different helicity configurations are:" << endl;
	for (int i=0;i<4;i++) {

		cout << "e-["<<hc[i][1]<<"]" << " + e+["<<hc[i][2]<<"] --> " <<"mu-["<<hc[i][3]<<"]" << " + mu+["<<hc[i][4]<<"] : " << ans[i] << endl;
		tans = tans + ans[i];
	}
	
	cout << "Total amplitude is : " << tans << endl;	
	return 0;
}
