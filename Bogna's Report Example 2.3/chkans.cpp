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
	double p_A_p_B[4] = {0}, p_A_p_D[4] = {0};
	double s = 0, u = 0, ans = 0;
	double eta[4][4]= {{-1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};


	ifstream B;				  //Object of ofstream class to write covariant four momentum of the scatterers to a file
	B.open("four_momentum_data_.dat", ios::in);	
	for (int i=1;i<5;i++) {
		for (int j=0;j<4;j++)
			B >> p[i][j];
	}
	B.close();
	
	for (int i=0;i<4;i++) {
		p_A_p_B[i] = p[1][i] - p[2][i];
		p_A_p_D[i] = p[1][i] - p[4][i];
	}
		
	for (int i=0;i<4;i++) {
		for (int j=0;j<4;j++)
			s = s + eta[i][j]*p_A_p_B[i]*p_A_p_B[j];
	}
	for (int i=0;i<4;i++) {
		for (int j=0;j<4;j++)
			u = u + eta[i][j]*p_A_p_D[i]*p_A_p_D[j];
	}
	ans = (s - u)*(s - u)/(s*u);
	cout << "Norm of the Feynman Amplitude of A4(12+34-) process using Mandelstam variables is " << ans << endl;
	return 0;
}
