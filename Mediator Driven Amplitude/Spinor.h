class Spinor {

	public:
	complex<double> squareket[2],squarebra[2],angleket[2],anglebra[2]; //defining the different types of Weyl spinors
	double theta,phi;

	//defining the gamma matrices in Weyl representation;
	//         row 1	      //        row 2 		//        row 3 		 //       row 4 
	complex<double> gamma[4][4][4] =
        {{{{0,0},{0,0},{1,0},{0,0}},{{0,0},{0,0},{0,0},{1,0}},{{1,0},{0,0},{0,0},{0,0}}, {{0,0},{1,0},{0,0},{0,0}}},   /*gamma 0*/
        {{{0,0},{0,0},{0,0},{1,0}},{{0,0},{0,0},{1,0},{0,0}},{{0,0},{-1,0},{0,0},{0,0}},{{-1,0},{0,0},{0,0},{0,0}}},   /*gamma 1*/
        {{{0,0},{0,0},{0,0},{0,-1}},{{0,0},{0,0},{0,1},{0,0}},{{0,0},{0,1},{0,0},{0,0}},{{0,-1},{0,0},{0,0},{0,0}}},   /*gamma 2*/
        {{{0,0},{0,0},{1,0},{0,0}},{{0,0},{0,0},{0,0},{-1,0}},{{-1,0},{0,0},{0,0},{0,0}},{{0,0},{1,0},{0,0},{0,0}}}};   /*gamma 3*/
	double diag_eta[4]= {1,-1,-1,-1};
    double diag_g5[4] = {-1,-1,1,1};
    double diag_id[4] = {1,1,1,1};

	double e = 0.001;		//precision parameter

	X IFWFN(double p[],double FMASS,int NHEL,int NSF);
    X OFWFN(double p[],double FMASS,int NHEL,int NSF);
    Y P_AMP(X IF[2],X OF[2],Y spp); //function to compute the amplitude for photon propagator
    Y Z_AMP(X IF[2],X OF[2],Y spp); //function to compute the amplitude for z boson propagator
    Y S_AMP(X IF[2],X OF[2],Y spp); //function to compute the amplitude for scalar boson propagator
};