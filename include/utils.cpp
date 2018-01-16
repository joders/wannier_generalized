#ifndef _UTILS_CPP_
#define _UTILS_CPP_

#include<string>
#include<complex>
#include<fstream>
#include<iostream>
#include<cstdlib>

using namespace std;



class parameters{}; //declare a general class here, through which all variables except for the one variable are generally passed to the function
// the specific classes for a given later function can be derived from this, e.g. class specific_parameters: public parameters{ double x;... };


#ifndef PI
#define PI 3.1415926535897932384626433
#endif

#ifndef cd
#define cd complex <double>
#endif

#ifndef SIGN
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#endif

/*contains the following small utility functions*/

bool fileExists(const std::string &fileName);
string inttostring(int i);
string doubletostring(double i);
string replace_dot_by_strdot(const std::string &str);
void write_fct_to_file(const double x[], const double f[], int length, const std::string &fileName);
void normalize_dist_trapez(double x[], double P[], int length);
double J_Greiner(double s_lattice_height);
double U_Greiner(double s_lattice_height, double lattice_laser_wavelength);
double block_int(double x[], double f[], int length);
double trapez_int(double x[], double f[], int length);
double lin_interp(double xa[], double ya[], int n, double x);
void splint(double xa[], double ya[], double y2a[], int n, double x, double *y);
void splint_deriv(double xa[], double ya[], double y2a[], int n, double x, double *yd);
void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]);
double calc_J_analytic_approx(double V_lat);
double calc_U_analytic_approx(double V_lat, double qlat_times_as);
double Ridder_find_val(double(*func)(double), double func_val, double x1, double x2, double x_accuracy);
//double Ridder_find_val(double(*func)(double, parameters &), double func_val, double x1, double x2, double x_accuracy, parameters &A);  ///same algorithm as above, but with the possibility to pass an object to contain all additional function parameters
void apply_b_Operator_to_state(cd state[], int N, double &norm);
void apply_b_dag_Operator_to_state(cd state[], int N, double &norm);
double num_recipes_ran1(long &idum);
void read_matrix_from_file(const std::string &fileName, double**&A, int &m, int &n);
void disp(cd z){if(imag(z)>0.0) printf("%f+%fi",real(z),imag(z));else printf("%f%fi",real(z),imag(z));}


template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
}

class tensor4{
public:
	/// class for a tensor of fourth order, i.e. a double-valued quantity with four indices, each index within the specified range 0...(dim-1)
	/// The memory is dynamically allocated and deleted by the constructor / desctructor of the class
	tensor4(int dim_tmp=20);
	tensor4(const std::string &fileName);		//read values from file immediately, size chosen according to data in file
	tensor4(const tensor4& tensor_to_be_copied); //copy constructor, needed because array a is dynamically allocated
	~tensor4();
	void operator= (const tensor4& tensor_to_be_copied);
	double getval(int i1,int i2,int i3,int i4);
	void setval(int i1,int i2,int i3,int i4, double val);
	void disp();
	int dim;
	double ****A;
};


//********************************************************************************
bool fileExists(const std::string &fileName) {
	std::fstream fin;
	fin.open(fileName.c_str(), std::ios::in);
	if (fin.is_open()) {
		fin.close();
		return (true);
	}

	fin.close();
	return (false);
}
//********************************************************************************
string inttostring(int i) // convert int to string
{
	stringstream s;
	s << i;
	return s.str();
}
//********************************************************************************

string doubletostring(double i) // convert double to string
{
	stringstream s;
	s << i;
	return s.str();
}
//********************************************************************************
string replace_dot_by_strdot(const std::string &str)
{
	int pos = str.find(".");
	if(pos==-1)
	{
		return str;
	}
	else
	{
		string out_str=str.substr(0,pos);
		out_str.append("dot");
		out_str.append(str.substr(pos+1));
		return replace_dot_by_strdot(out_str);
	}
}
//********************************************************************************
string replace_minus_by_strminus(const std::string &str)
{
	int pos = str.find("-");
	if(pos==-1)
	{
		return str;
	}
	else
	{
		string out_str=str.substr(0,pos);
		out_str.append("m");
		out_str.append(str.substr(pos+1));
		return replace_dot_by_strdot(out_str);
	}
}
//********************************************************************************
void write_fct_to_file(const double x[], const double f[], int length,
const std::string &fileName) {
	FILE *F = fopen(fileName.c_str(), "w");
	for (int c = 0; c < length; c++)
	fprintf(F, "%g %g\n", x[c], f[c]);
	fclose(F);
}
//********************************************************************************

double trapez_int(double x[], double f[], int length) {
	double sum = 0.0;
	for (int i = 0; i <= (length - 2); i++) {
		sum = sum + (f[i] + f[i + 1]) * (x[i + 1] - x[i]) / 2;
	}

	return (sum);
}
//********************************************************************************

void normalize_dist_trapez(double x[], double P[], int length) {
	double norm = trapez_int(x, P, length);
	for (int c = 0; c < length; c++) {
		P[c] = P[c] / norm;
	}
}

//********************************************************************************
double J_over_U_Greiner(double s_lattice_height,
double lattice_laser_wavelength) {
	/* interpolates the value of J/U for a given lattice strength s=V_lat/E_R in 3D
	from data points extracted from Greiners PhD thesis , pg. 80

	the ratio depends on the lattice laser wavelength, although this is not
	explicitly stated in the thesis
	value J/U for fixed s scales with lambda^3 in three dimensions
	In Greiner's thesis the value plotted is for lambda_lat=852nm

	in units of E_R J/U only scales with lambda though

	The second input is in units of nm
	*/

	// extracted points:
	const int length = 19;
	double Greiner_wavelength = 852; //in nm for the laser this data is plotted in the thesis

	//these two arrays are points extracted from y-logarithmic plot in Geriner's thesis
	double s_given[length] = { 1.01414940130000, 2.40064048610000,
		3.93939309630000, 5.32372789630000, 6.76913424470000,
		8.42905614550000, 10.18093315660000, 11.84022903920000,
		13.56094425790000, 15.55780303800000, 18.26087984760000,
		21.39389201050000, 24.98709707040000, 28.24218277230000,
		31.98848404840000, 37.88376651470000, 42.91883048680000,
		48.35273757740000, 50.01043363580000 };
	double U_over_J_given[length] = { 0.1664756656, 0.4901510852, 1.1859962469,
		2.3872931474, 4.5753484417, 8.4524490669, 15.4249270455,
		25.5171852752, 42.7342995634, 72.456519298, 138.9088442538,
		290.216133622, 613.8890277378, 1148.539296627, 2313.236143328,
		6334.521566543, 13736.97262781, 30160.48873936, 37628.10555965 };
	double ln_U_over_J_given[length];
	for (int i = 0; i < length; i++) {
		ln_U_over_J_given[i] = log(U_over_J_given[i]);
	}

	//these values were calculated with routine below using the fct "spline" to calculate the second derivative
	double second_deriv_for_spline[length] = { -0.5, -0.238606379458,
		0.00191902695201, -0.041487459281, -0.0689085724149,
		0.00368877500038, -0.0394868514368, 0.0138182490967,
		-0.0290768115456, -0.00960069334907, 0.0024072541182,
		-0.011407026038, -0.00393380127385, -0.000447864345618,
		-0.0038038657623, -0.0035736065063, -0.000455014676509,
		-0.00488251609211, 0.00244125804606 };

	/*
	double second_deriv_for_spline[length];
	double yp1=(ln_U_over_J_given[1]-ln_U_over_J_given[0])/(s_given[1]-s_given[0]);
	double ypn=(ln_U_over_J_given[length-1]-ln_U_over_J_given[length-2])/(s_given[length-1]-s_given[length-2]);
	spline(s_given, ln_U_over_J_given, length, yp1, ypn, second_deriv_for_spline);*/

	//the minus in the final exponential is due to the fact that J/U is returned, but the data for U/J is tabulated
	if (s_lattice_height <= s_given[0] || s_lattice_height >= s_given[length
			- 1]) {
		return exp(-lin_interp(s_given, ln_U_over_J_given, length,s_lattice_height)) *lattice_laser_wavelength* lattice_laser_wavelength / (Greiner_wavelength* Greiner_wavelength);
	} else {
		double value;
		double *val_p = &value;
		splint(s_given, ln_U_over_J_given, second_deriv_for_spline, length,
		s_lattice_height, val_p);
		return exp(-value) * lattice_laser_wavelength / Greiner_wavelength;
	}
}

//********************************************************************************
double J_Greiner(double s_lattice_height) {
	/*  in E_R
	interpolates the value of J/E_R for a given lattice strength s=V_lat/E_R in 3D
	from data points extracted from Greiners PhD thesis , pg. 80

	the ratio does not depend on the lattice laser wavelength
	*/

	// extracted points:
	const int length = 499;
	//these two arrays are points extracted from y-logarithmic plot in Geriner's thesis
	double s_given[length] = { 1.2196764e+000, 1.3158992e+000, 1.4121475e+000,
		1.5083704e+000, 1.6045932e+000, 1.7008415e+000, 1.7970898e+000,
		1.8933126e+000, 1.9895354e+000, 2.0857837e+000, 2.1820065e+000,
		2.2782548e+000, 2.3744776e+000, 2.4707005e+000, 2.5669488e+000,
		2.6631971e+000, 2.7594199e+000, 2.8556427e+000, 2.9518910e+000,
		3.0481138e+000, 3.1443621e+000, 3.2405849e+000, 3.3368077e+000,
		3.4330306e+000, 3.5292789e+000, 3.6255017e+000, 3.7217500e+000,
		3.8179728e+000, 3.9141956e+000, 4.0104184e+000, 4.1066667e+000,
		4.2028895e+000, 4.2991379e+000, 4.3953607e+000, 4.4915835e+000,
		4.5878063e+000, 4.6840291e+000, 4.7802774e+000, 4.8765002e+000,
		4.9727485e+000, 5.0689713e+000, 5.1651942e+000, 5.2614170e+000,
		5.3576398e+000, 5.4538626e+000, 5.5500854e+000, 5.6463082e+000,
		5.7425310e+000, 5.8387539e+000, 5.9349767e+000, 6.0311995e+000,
		6.1274223e+000, 6.2236451e+000, 6.3198679e+000, 6.4160907e+000,
		6.5123136e+000, 6.6085364e+000, 6.7047592e+000, 6.8009820e+000,
		6.8972048e+000, 6.9934276e+000, 7.0896504e+000, 7.1858733e+000,
		7.2820961e+000, 7.3783189e+000, 7.4745417e+000, 7.5707645e+000,
		7.6669873e+000, 7.7632101e+000, 7.8594330e+000, 7.9556558e+000,
		8.0518531e+000, 8.1480504e+000, 8.2442732e+000, 8.3404961e+000,
		8.4367189e+000, 8.5329417e+000, 8.6291645e+000, 8.7253873e+000,
		8.8216101e+000, 8.9178329e+000, 9.0140557e+000, 9.1102531e+000,
		9.2064759e+000, 9.3026732e+000, 9.3988960e+000, 9.4950933e+000,
		9.5912907e+000, 9.6875135e+000, 9.7837108e+000, 9.8799336e+000,
		9.9761309e+000, 1.0072354e+001, 1.0168551e+001, 1.0264748e+001,
		1.0360971e+001, 1.0457169e+001, 1.0553366e+001, 1.0649589e+001,
		1.0745811e+001, 1.0842009e+001, 1.0938206e+001, 1.1034429e+001,
		1.1130626e+001, 1.1226824e+001, 1.1323046e+001, 1.1419244e+001,
		1.1515441e+001, 1.1611638e+001, 1.1707836e+001, 1.1804059e+001,
		1.1900256e+001, 1.1996453e+001, 1.2092650e+001, 1.2188848e+001,
		1.2285071e+001, 1.2381268e+001, 1.2477465e+001, 1.2573663e+001,
		1.2669860e+001, 1.2766083e+001, 1.2862280e+001, 1.2958477e+001,
		1.3054675e+001, 1.3150872e+001, 1.3247069e+001, 1.3343267e+001,
		1.3439464e+001, 1.3535661e+001, 1.3631859e+001, 1.3728056e+001,
		1.3824253e+001, 1.3920425e+001, 1.4016648e+001, 1.4112820e+001,
		1.4209017e+001, 1.4305214e+001, 1.4401412e+001, 1.4497609e+001,
		1.4593781e+001, 1.4690004e+001, 1.4786176e+001, 1.4882373e+001,
		1.4978570e+001, 1.5074768e+001, 1.5170965e+001, 1.5267162e+001,
		1.5363359e+001, 1.5459531e+001, 1.5555729e+001, 1.5651926e+001,
		1.5748098e+001, 1.5844295e+001, 1.5940492e+001, 1.6036664e+001,
		1.6132887e+001, 1.6229059e+001, 1.6325256e+001, 1.6421454e+001,
		1.6517625e+001, 1.6613823e+001, 1.6710020e+001, 1.6806192e+001,
		1.6902389e+001, 1.6998587e+001, 1.7094758e+001, 1.7190981e+001,
		1.7287153e+001, 1.7383350e+001, 1.7479548e+001, 1.7575719e+001,
		1.7671891e+001, 1.7768063e+001, 1.7864286e+001, 1.7960458e+001,
		1.8056630e+001, 1.8152801e+001, 1.8249024e+001, 1.8345196e+001,
		1.8441368e+001, 1.8537540e+001, 1.8633763e+001, 1.8729934e+001,
		1.8826106e+001, 1.8922278e+001, 1.9018475e+001, 1.9114673e+001,
		1.9210845e+001, 1.9307016e+001, 1.9403188e+001, 1.9499411e+001,
		1.9595583e+001, 1.9691755e+001, 1.9787927e+001, 1.9884098e+001,
		1.9980270e+001, 2.0076442e+001, 2.0172665e+001, 2.0268837e+001,
		2.0365009e+001, 2.0461180e+001, 2.0557352e+001, 2.0653524e+001,
		2.0749696e+001, 2.0845919e+001, 2.0942090e+001, 2.1038262e+001,
		2.1134434e+001, 2.1230606e+001, 2.1326778e+001, 2.1422950e+001,
		2.1519172e+001, 2.1615344e+001, 2.1711516e+001, 2.1807688e+001,
		2.1903860e+001, 2.2000032e+001, 2.2096203e+001, 2.2192375e+001,
		2.2288547e+001, 2.2384719e+001, 2.2480891e+001, 2.2577063e+001,
		2.2673234e+001, 2.2769432e+001, 2.2865629e+001, 2.2961801e+001,
		2.3057973e+001, 2.3154145e+001, 2.3250316e+001, 2.3346488e+001,
		2.3442660e+001, 2.3538832e+001, 2.3635004e+001, 2.3731176e+001,
		2.3827347e+001, 2.3923519e+001, 2.4019691e+001, 2.4115863e+001,
		2.4212035e+001, 2.4308207e+001, 2.4404378e+001, 2.4500550e+001,
		2.4596722e+001, 2.4692894e+001, 2.4789066e+001, 2.4885238e+001,
		2.4981409e+001, 2.5077581e+001, 2.5173753e+001, 2.5269925e+001,
		2.5366097e+001, 2.5462269e+001, 2.5558440e+001, 2.5654612e+001,
		2.5750784e+001, 2.5846956e+001, 2.5943128e+001, 2.6039300e+001,
		2.6135471e+001, 2.6231643e+001, 2.6327815e+001, 2.6423961e+001,
		2.6520133e+001, 2.6616280e+001, 2.6712451e+001, 2.6808623e+001,
		2.6904795e+001, 2.7000967e+001, 2.7097139e+001, 2.7193311e+001,
		2.7289482e+001, 2.7385654e+001, 2.7481826e+001, 2.7577998e+001,
		2.7674170e+001, 2.7770342e+001, 2.7866488e+001, 2.7962660e+001,
		2.8058806e+001, 2.8154978e+001, 2.8251150e+001, 2.8347322e+001,
		2.8443493e+001, 2.8539665e+001, 2.8635837e+001, 2.8732009e+001,
		2.8828155e+001, 2.8924327e+001, 2.9020499e+001, 2.9116645e+001,
		2.9212817e+001, 2.9308989e+001, 2.9405161e+001, 2.9501333e+001,
		2.9597504e+001, 2.9693651e+001, 2.9789797e+001, 2.9885969e+001,
		2.9982141e+001, 3.0078313e+001, 3.0174459e+001, 3.0270631e+001,
		3.0366777e+001, 3.0462949e+001, 3.0559121e+001, 3.0655293e+001,
		3.0751439e+001, 3.0847585e+001, 3.0943757e+001, 3.1039929e+001,
		3.1136101e+001, 3.1232247e+001, 3.1328393e+001, 3.1424565e+001,
		3.1520737e+001, 3.1616909e+001, 3.1713055e+001, 3.1809227e+001,
		3.1905374e+001, 3.2001545e+001, 3.2097717e+001, 3.2193864e+001,
		3.2290035e+001, 3.2386182e+001, 3.2482354e+001, 3.2578525e+001,
		3.2674672e+001, 3.2770844e+001, 3.2866990e+001, 3.2963162e+001,
		3.3059334e+001, 3.3155505e+001, 3.3251652e+001, 3.3347798e+001,
		3.3443970e+001, 3.3540142e+001, 3.3636314e+001, 3.3732460e+001,
		3.3828606e+001, 3.3924778e+001, 3.4020950e+001, 3.4117096e+001,
		3.4213243e+001, 3.4309414e+001, 3.4405561e+001, 3.4501733e+001,
		3.4597879e+001, 3.4694051e+001, 3.4790197e+001, 3.4886369e+001,
		3.4982515e+001, 3.5078687e+001, 3.5174833e+001, 3.5271005e+001,
		3.5367152e+001, 3.5463323e+001, 3.5559470e+001, 3.5655616e+001,
		3.5751788e+001, 3.5847934e+001, 3.5944081e+001, 3.6040252e+001,
		3.6136399e+001, 3.6232571e+001, 3.6328717e+001, 3.6424889e+001,
		3.6521035e+001, 3.6617182e+001, 3.6713353e+001, 3.6809500e+001,
		3.6905672e+001, 3.7001818e+001, 3.7097964e+001, 3.7194136e+001,
		3.7290282e+001, 3.7386454e+001, 3.7482601e+001, 3.7578747e+001,
		3.7674919e+001, 3.7771065e+001, 3.7867211e+001, 3.7963383e+001,
		3.8059530e+001, 3.8155676e+001, 3.8251848e+001, 3.8347994e+001,
		3.8444140e+001, 3.8540312e+001, 3.8636459e+001, 3.8732605e+001,
		3.8828777e+001, 3.8924923e+001, 3.9021069e+001, 3.9117241e+001,
		3.9213388e+001, 3.9309534e+001, 3.9405706e+001, 3.9501852e+001,
		3.9597999e+001, 3.9694145e+001, 3.9790317e+001, 3.9886463e+001,
		3.9982609e+001, 4.0078756e+001, 4.0174928e+001, 4.0271074e+001,
		4.0367220e+001, 4.0463367e+001, 4.0559513e+001, 4.0655685e+001,
		4.0751831e+001, 4.0847977e+001, 4.0944124e+001, 4.1040270e+001,
		4.1136442e+001, 4.1232588e+001, 4.1328735e+001, 4.1424881e+001,
		4.1521027e+001, 4.1617199e+001, 4.1713345e+001, 4.1809492e+001,
		4.1905638e+001, 4.2001784e+001, 4.2097931e+001, 4.2194103e+001,
		4.2290249e+001, 4.2386395e+001, 4.2482542e+001, 4.2578688e+001,
		4.2674834e+001, 4.2770981e+001, 4.2867153e+001, 4.2963299e+001,
		4.3059445e+001, 4.3155592e+001, 4.3251738e+001, 4.3347910e+001,
		4.3444056e+001, 4.3540202e+001, 4.3636349e+001, 4.3732495e+001,
		4.3828641e+001, 4.3924788e+001, 4.4020960e+001, 4.4117106e+001,
		4.4213252e+001, 4.4309399e+001, 4.4405545e+001, 4.4501691e+001,
		4.4597838e+001, 4.4693984e+001, 4.4790130e+001, 4.4886277e+001,
		4.4982423e+001, 4.5078569e+001, 4.5174716e+001, 4.5270862e+001,
		4.5367008e+001, 4.5463155e+001, 4.5559301e+001, 4.5655447e+001,
		4.5751594e+001, 4.5847766e+001, 4.5943886e+001, 4.6040058e+001,
		4.6136179e+001, 4.6232325e+001, 4.6328472e+001, 4.6424618e+001,
		4.6520765e+001, 4.6616911e+001, 4.6713057e+001, 4.6809178e+001,
		4.6905350e+001, 4.7001471e+001, 4.7097617e+001, 4.7193763e+001,
		4.7289910e+001, 4.7386056e+001, 4.7482177e+001, 4.7578349e+001,
		4.7674470e+001, 4.7770616e+001, 4.7866762e+001, 4.7962909e+001,
		4.8059055e+001, 4.8155176e+001, 4.8251348e+001, 4.8347469e+001,
		4.8443615e+001, 4.8539761e+001, 4.8635882e+001, 4.8732054e+001,
		4.8828175e+001, 4.8924321e+001, 4.9020467e+001, 5.0000000e+001 };
	double J_given[length] = { 1.7988366e-001, 1.7537195e-001, 1.7105184e-001,
		1.6668160e-001, 1.6229162e-001, 1.5814537e-001, 1.5409860e-001,
		1.5003353e-001, 1.4619417e-001, 1.4246205e-001, 1.3881946e-001,
		1.3526427e-001, 1.3180838e-001, 1.2833687e-001, 1.2505810e-001,
		1.2185799e-001, 1.1864343e-001, 1.1560733e-001, 1.1265605e-001,
		1.0977557e-001, 1.0696419e-001, 1.0431150e-001, 1.0164421e-001,
		9.9127479e-002, 9.6596901e-002, 9.4127034e-002, 9.1716420e-002,
		8.9441867e-002, 8.7154805e-002, 8.4996831e-002, 8.2826988e-002,
		8.0709200e-002, 7.8642220e-002, 7.6691906e-002, 7.4788381e-002,
		7.2938026e-002, 7.1132063e-002, 6.9316167e-002, 6.7543835e-002,
		6.5814023e-002, 6.4181846e-002, 6.2588826e-002, 6.1087282e-002,
		5.9621762e-002, 5.8191400e-002, 5.6795353e-002, 5.5432799e-002,
		5.4102933e-002, 5.2804971e-002, 5.1538148e-002, 5.0301717e-002,
		4.9094949e-002, 4.7917131e-002, 4.6767571e-002, 4.5645589e-002,
		4.4550524e-002, 4.3481730e-002, 4.2438577e-002, 4.1420450e-002,
		4.0426749e-002, 3.9456887e-002, 3.8510293e-002, 3.7586408e-002,
		3.6684687e-002, 3.5804600e-002, 3.4945626e-002, 3.4107259e-002,
		3.3289006e-002, 3.2490382e-002, 3.1735405e-002, 3.1022536e-002,
		3.0326768e-002, 2.9647693e-002, 2.8984909e-002, 2.8313539e-002,
		2.7634281e-002, 2.6971318e-002, 2.6324260e-002, 2.5692725e-002,
		2.5076341e-002, 2.4493643e-002, 2.3924925e-002, 2.3387926e-002,
		2.2863810e-002, 2.2369557e-002, 2.1885204e-002, 2.1412472e-002,
		2.0949210e-002, 2.0497062e-002, 2.0053971e-002, 1.9604568e-002,
		1.9180770e-002, 1.8765463e-002, 1.8360119e-002, 1.7962895e-002,
		1.7588816e-002, 1.7208888e-002, 1.6823548e-002, 1.6460163e-002,
		1.6104057e-002, 1.5742879e-002, 1.5402279e-002, 1.5081526e-002,
		1.4755757e-002, 1.4436514e-002, 1.4135873e-002, 1.3841254e-002,
		1.3552538e-002, 1.3270748e-002, 1.2994602e-002, 1.2723989e-002,
		1.2458797e-002, 1.2198917e-002, 1.1945272e-002, 1.1696708e-002,
		1.1453124e-002, 1.1214419e-002, 1.0980496e-002, 1.0752185e-002,
		1.0528448e-002, 1.0309192e-002, 1.0094329e-002, 9.8837700e-003,
		9.6857732e-003, 9.4917429e-003, 9.3015994e-003, 9.1152650e-003,
		8.9326633e-003, 8.7537197e-003, 8.5783607e-003, 8.4130597e-003,
		8.2446561e-003, 8.0859117e-003, 7.9303472e-003, 7.7778992e-003,
		7.6219598e-003, 7.4755584e-003, 7.3315080e-003, 7.1847537e-003,
		7.0464168e-003, 6.9108512e-003, 6.7780013e-003, 6.6421089e-003,
		6.5145283e-003, 6.3840260e-003, 6.2561379e-003, 6.1355851e-003,
		6.0174473e-003, 5.9016760e-003, 5.7927344e-003, 5.6859752e-003,
		5.5765813e-003, 5.4736407e-003, 5.3685857e-003, 5.2652178e-003,
		5.1639205e-003, 5.0646525e-003, 4.9669786e-003, 4.8754379e-003,
		4.7816381e-003, 4.6933718e-003, 4.6068737e-003, 4.5182408e-003,
		4.4348366e-003, 4.3497193e-003, 4.2659688e-003, 4.1838961e-003,
		4.1034675e-003, 4.0243304e-003, 3.9532458e-003, 3.8833026e-003,
		3.8117877e-003, 3.7446188e-003, 3.6785285e-003, 3.6134995e-003,
		3.5439755e-003, 3.4815258e-003, 3.4200789e-003, 3.3596188e-003,
		3.2949795e-003, 3.2369174e-003, 3.1797877e-003, 3.1235755e-003,
		3.0658669e-003, 3.0118842e-003, 2.9587684e-003, 2.9065056e-003,
		2.8550820e-003, 2.8025028e-003, 2.7531188e-003, 2.7045279e-003,
		2.6567173e-003, 2.6096744e-003, 2.5677673e-003, 2.5265331e-003,
		2.4820816e-003, 2.4383439e-003, 2.3953086e-003, 2.3529643e-003,
		2.3113001e-003, 2.2741843e-003, 2.2376646e-003, 2.1982954e-003,
		2.1595585e-003, 2.1214435e-003, 2.0839407e-003, 2.0470400e-003,
		2.0141679e-003, 1.9818236e-003, 1.9469557e-003, 1.9126476e-003,
		1.8788905e-003, 1.8456755e-003, 1.8129939e-003, 1.7838801e-003,
		1.7552339e-003, 1.7270477e-003, 1.6993142e-003, 1.6720259e-003,
		1.6451759e-003, 1.6187571e-003, 1.5927625e-003, 1.5659600e-003,
		1.5383871e-003, 1.5112570e-003, 1.4845626e-003, 1.4582968e-003,
		1.4336782e-003, 1.4106557e-003, 1.3880028e-003, 1.3657137e-003,
		1.3437826e-003, 1.3222036e-003, 1.3009712e-003, 1.2800797e-003,
		1.2595237e-003, 1.2392978e-003, 1.2193967e-003, 1.1998152e-003,
		1.1805481e-003, 1.1615904e-003, 1.1429372e-003, 1.1245834e-003,
		1.1065245e-003, 1.0887555e-003, 1.0712718e-003, 1.0540690e-003,
		1.0371423e-003, 1.0204875e-003, 1.0041001e-003, 9.8797592e-004,
		9.7211063e-004, 9.5650011e-004, 9.4114028e-004, 9.2602709e-004,
		9.1115660e-004, 8.9652491e-004, 8.8212818e-004, 8.6796263e-004,
		8.5402456e-004, 8.4096999e-004, 8.2812505e-004, 8.1612761e-004,
		8.0432284e-004, 7.9270764e-004, 7.8061928e-004, 7.6872505e-004,
		7.5638057e-004, 7.4423433e-004, 7.3228314e-004, 7.2052386e-004,
		7.0895342e-004, 6.9756878e-004, 6.8636696e-004, 6.7534502e-004,
		6.6502173e-004, 6.5486422e-004, 6.4537690e-004, 6.3604193e-004,
		6.2685686e-004, 6.1729764e-004, 6.0789192e-004, 5.9813017e-004,
		5.8852517e-004, 5.7907441e-004, 5.7022271e-004, 5.6151315e-004,
		5.5294345e-004, 5.4493919e-004, 5.3706345e-004, 5.2886690e-004,
		5.2080197e-004, 5.1286656e-004, 5.0463075e-004, 4.9691699e-004,
		4.8971219e-004, 4.8262308e-004, 4.7564781e-004, 4.6878455e-004,
		4.6200412e-004, 4.5494749e-004, 4.4835645e-004, 4.4187126e-004,
		4.3549021e-004, 4.2884922e-004, 4.2264641e-004, 4.1651846e-004,
		4.1048891e-004, 4.0455618e-004, 3.9871872e-004, 3.9295171e-004,
		3.8725430e-004, 3.8164838e-004, 3.7613248e-004, 3.7070516e-004,
		3.6534334e-004, 3.5976310e-004, 3.5455104e-004, 3.4942268e-004,
		3.4437668e-004, 3.3939157e-004, 3.3448652e-004, 3.2964066e-004,
		3.2487261e-004, 3.2018113e-004, 3.1554627e-004, 3.1098584e-004,
		3.0648045e-004, 3.0204740e-004, 2.9768554e-004, 2.9314600e-004,
		2.8890598e-004, 2.8471713e-004, 2.8059554e-004, 2.7654014e-004,
		2.7254987e-004, 2.6860775e-004, 2.6471320e-004, 2.6088119e-004,
		2.5711072e-004, 2.5359989e-004, 2.5013143e-004, 2.4651051e-004,
		2.4313893e-004, 2.3982150e-004, 2.3654410e-004, 2.3312263e-004,
		2.2993678e-004, 2.2661089e-004, 2.2351403e-004, 2.2028104e-004,
		2.1727068e-004, 2.1412799e-004, 2.1120172e-004, 2.0814682e-004,
		2.0530229e-004, 2.0249208e-004, 1.9972700e-004, 1.9699529e-004,
		1.9446049e-004, 1.9180508e-004, 1.8918172e-004, 1.8660049e-004,
		1.8405041e-004, 1.8138824e-004, 1.7890938e-004, 1.7646045e-004,
		1.7405084e-004, 1.7167031e-004, 1.6932800e-004, 1.6701397e-004,
		1.6472786e-004, 1.6247846e-004, 1.6025620e-004, 1.5806964e-004,
		1.5590946e-004, 1.5377535e-004, 1.5167551e-004, 1.4960101e-004,
		1.4767604e-004, 1.4565949e-004, 1.4366728e-004, 1.4181866e-004,
		1.3988209e-004, 1.3796890e-004, 1.3619361e-004, 1.3433385e-004,
		1.3249654e-004, 1.3079166e-004, 1.2900567e-004, 1.2724123e-004,
		1.2560398e-004, 1.2388882e-004, 1.2219437e-004, 1.2062206e-004,
		1.1897493e-004, 1.1734769e-004, 1.1583774e-004, 1.1434601e-004,
		1.1278330e-004, 1.1133324e-004, 1.0990068e-004, 1.0848541e-004,
		1.0700280e-004, 1.0562706e-004, 1.0426792e-004, 1.0292519e-004,
		1.0159866e-004, 1.0029343e-004, 9.9003953e-005, 9.7730037e-005,
		9.6471495e-005, 9.5228142e-005, 9.4004754e-005, 9.2796131e-005,
		9.1602093e-005, 9.0422466e-005, 8.9257074e-005, 8.8110396e-005,
		8.6977557e-005, 8.5858389e-005, 8.4752727e-005, 8.3660408e-005,
		8.2650785e-005, 8.1588977e-005, 8.0539984e-005, 7.9503651e-005,
		7.8479823e-005, 7.7468352e-005, 7.6533455e-005, 7.5609840e-005,
		7.4638487e-005, 7.3678856e-005, 7.2730807e-005, 7.1794199e-005,
		7.0868894e-005, 6.9958447e-005, 6.9058988e-005, 6.8170384e-005,
		6.7292504e-005, 6.6425218e-005, 6.5623590e-005, 6.4831637e-005,
		6.3998751e-005, 6.3175916e-005, 6.2363011e-005, 6.1559917e-005,
		6.0766514e-005, 6.0033177e-005, 5.9308689e-005, 5.8592945e-005,
		5.7885838e-005, 5.7187265e-005, 5.6497122e-005, 5.5815308e-005,
		5.5141722e-005, 5.4476265e-005, 5.3818839e-005, 5.3169347e-005,
		5.2527693e-005, 5.1893783e-005, 5.1267522e-005, 5.0608893e-005,
		4.9997657e-005, 4.9354829e-005, 4.8758257e-005, 4.8168884e-005,
		4.7626551e-005, 4.7051318e-005, 4.6521995e-005, 4.5960562e-005,
		4.5405905e-005, 4.4893446e-005, 4.4352095e-005, 4.3851930e-005,
		4.3357802e-005, 4.2869636e-005, 4.2351857e-005, 4.1875403e-005,
		4.1402790e-005, 4.0903532e-005, 4.0442257e-005, 3.9986549e-005,
		3.9536341e-005, 3.9058821e-005, 3.8619413e-005, 3.8183548e-005,
		3.7723109e-005, 3.7297700e-005, 3.6877426e-005, 3.6462223e-005,
		3.6050367e-005, 3.5645125e-005, 3.5243150e-005, 3.4846027e-005,
		3.4453695e-005, 3.4064527e-005, 3.3708245e-005, 3.0609842e-005 };

	double ln_J_given[length];
	for (int i = 0; i < length; i++) {
		ln_J_given[i] = log(J_given[i]);
	}

	if (s_lattice_height <= 50)
	return exp(lin_interp(s_given, ln_J_given, length, s_lattice_height));
	else
	return 4.0 / sqrt(PI) * pow(s_lattice_height, 0.75) * exp(-2.0 * pow(
	s_lattice_height, 0.5));
}

//********************************************************************************
double U_Greiner(double s_lattice_height, double lattice_laser_wavelength) {
	/*  in E_R interpolates the value of U/E_R for a given lattice strength s=V_lat/E_R in 3D
	from data points extracted from Greiners PhD thesis , pg. 80

	the ratio depends on the lattice laser wavelength, although this is not
	explicitly stated in the thesis
	value J/U for fixed s scales with lambda^3 in three dimensions
	In Greiner's thesis the value plotted is for lambda_lat=852nm
	//in units of E_R, U scales as 1/lambda

	The second input is in units of nm */

	// extracted points:
	double Greiner_wavelength = 852; //in nm for the laser this data is plotted in the thesis
	const int length = 98;
	//these two arrays are points extracted from y-logarithmic plot in Geriner's thesis
	double s_given[length] = { 1.5000000e+000, 2.0000000e+000, 2.5000000e+000,
		3.0000000e+000, 3.5000000e+000, 4.0000000e+000, 4.5000000e+000,
		5.0000000e+000, 5.5000000e+000, 6.0000000e+000, 6.5000000e+000,
		7.0000000e+000, 7.5000000e+000, 8.0000000e+000, 8.5000000e+000,
		9.0000000e+000, 9.5000000e+000, 1.0000000e+001, 1.0500000e+001,
		1.1000000e+001, 1.1500000e+001, 1.2000000e+001, 1.2500000e+001,
		1.3000000e+001, 1.3500000e+001, 1.4000000e+001, 1.4500000e+001,
		1.5000000e+001, 1.5500000e+001, 1.6000000e+001, 1.6500000e+001,
		1.7000000e+001, 1.7500000e+001, 1.8000000e+001, 1.8500000e+001,
		1.9000000e+001, 1.9500000e+001, 2.0000000e+001, 2.0500000e+001,
		2.1000000e+001, 2.1500000e+001, 2.2000000e+001, 2.2500000e+001,
		2.3000000e+001, 2.3500000e+001, 2.4000000e+001, 2.4500000e+001,
		2.5000000e+001, 2.5500000e+001, 2.6000000e+001, 2.6500000e+001,
		2.7000000e+001, 2.7500000e+001, 2.8000000e+001, 2.8500000e+001,
		2.9000000e+001, 2.9500000e+001, 3.0000000e+001, 3.0500000e+001,
		3.1000000e+001, 3.1500000e+001, 3.2000000e+001, 3.2500000e+001,
		3.3000000e+001, 3.3500000e+001, 3.4000000e+001, 3.4500000e+001,
		3.5000000e+001, 3.5500000e+001, 3.6000000e+001, 3.6500000e+001,
		3.7000000e+001, 3.7500000e+001, 3.8000000e+001, 3.8500000e+001,
		3.9000000e+001, 3.9500000e+001, 4.0000000e+001, 4.0500000e+001,
		4.1000000e+001, 4.1500000e+001, 4.2000000e+001, 4.2500000e+001,
		4.3000000e+001, 4.3500000e+001, 4.4000000e+001, 4.4500000e+001,
		4.5000000e+001, 4.5500000e+001, 4.6000000e+001, 4.6500000e+001,
		4.7000000e+001, 4.7500000e+001, 4.8000000e+001, 4.8500000e+001,
		4.9000000e+001, 4.9500000e+001, 5.0000000e+001 };
	double U_given[length] = { 4.2455921e-002, 5.5310302e-002, 6.8164682e-002,
		8.1019063e-002, 9.6062346e-002, 1.1070465e-001, 1.2574102e-001,
		1.4039024e-001, 1.5501209e-001, 1.6963395e-001, 1.8421404e-001,
		1.9861122e-001, 2.1261579e-001, 2.2635391e-001, 2.4009203e-001,
		2.5398908e-001, 2.6756827e-001, 2.8101476e-001, 2.9416077e-001,
		3.0789889e-001, 3.2075327e-001, 3.3360765e-001, 3.4627413e-001,
		3.5867647e-001, 3.7128705e-001, 3.8336935e-001, 3.9522834e-001,
		4.0729987e-001, 4.1948622e-001, 4.3158213e-001, 4.4341674e-001,
		4.5516092e-001, 4.6734727e-001, 4.7864928e-001, 4.9039346e-001,
		5.0169548e-001, 5.1343966e-001, 5.2429951e-001, 5.3560152e-001,
		5.4686986e-001, 5.5776338e-001, 5.6862323e-001, 5.7971965e-001,
		5.9078510e-001, 6.0164494e-001, 6.1250753e-001, 6.2336464e-001,
		6.3422449e-001, 6.4508434e-001, 6.5591693e-001, 6.6591971e-001,
		6.7633739e-001, 6.8658090e-001, 6.9678406e-001, 7.0730954e-001,
		7.1756596e-001, 7.2768959e-001, 7.3832209e-001, 7.4868362e-001,
		7.5879453e-001, 7.6908995e-001, 7.7908220e-001, 7.8882977e-001,
		7.9881570e-001, 8.0901887e-001, 8.1895739e-001, 8.2870496e-001,
		8.3859867e-001, 8.4820009e-001, 8.5826721e-001, 8.6847038e-001,
		8.7778980e-001, 8.8799296e-001, 8.9731239e-001, 9.0751555e-001,
		9.1683498e-001, 9.2618062e-001, 9.3614302e-001, 9.4567576e-001,
		9.5499642e-001, 9.6431584e-001, 9.7363527e-001, 9.8295469e-001,
		9.9245881e-001, 1.0015935e+000, 1.0112540e+000, 1.0202324e+000,
		1.0298642e+000, 1.0396117e+000, 1.0490744e+000, 1.0582219e+000,
		1.0670846e+000, 1.0761489e+000, 1.0854684e+000, 1.0945574e+000,
		1.1034201e+000, 1.1125429e+000, 1.1218623e+000 };

	return (lin_interp(s_given, U_given, length, s_lattice_height)
	* Greiner_wavelength / (lattice_laser_wavelength));
}

//********************************************************************************
double block_int(double x[], double f[], int length) {
	double sum = 0.0;
	for (int i = 0; i <= (length - 2); i++) {
		sum = sum + f[i] * (x[i + 1] - x[i]);
	}

	return (sum);
}

//********************************************************************************
double lin_interp(double xa[], double ya[], int n, double x)

// spline interpolation from numerical recipes in C
{
	void nrerror(char error_text[]);
	int klo, khi, k;
	double h;

	klo = 0;
	khi = n - 1;
	while (khi - klo > 1) {
		k = (khi + klo) >> 1;
		if (xa[k] > x) {
			khi = k;
		} else {
			klo = k;
		}
	}

	h = xa[khi] - xa[klo];
	if (h == 0.0) {
		printf(
		"\n\nBad xa input to routine splint khi=%d  klo=%d  xlo=%g xhi=%g   n=%d \n\n",
		khi, klo, xa[klo], xa[khi], n);
		abort();
	}

	return (ya[klo] + (x - xa[klo]) * (ya[khi] - ya[klo]) / (xa[khi] - xa[klo]));
}

//********************************************************************************
void splint(double xa[], double ya[], double y2a[], int n, double x, double *y)
// spline interpolation from numerical recipes in C
{
	void nrerror(char error_text[]);
	int klo, khi, k;
	double h, b, a;

	klo = 0;
	khi = n - 1;
	while (khi - klo > 1) {
		k = (khi + klo) >> 1;
		if (xa[k] > x) {
			khi = k;
		} else {
			klo = k;
		}
	}

	h = xa[khi] - xa[klo];
	if (h == 0.0) {
		printf(
		"\n\nBad xa input to routine splint khi=%d  klo=%d  xlo=%g  xhi=%g  n=%d \n\n",
		khi, klo, xa[klo], xa[khi], n);
		int te;
		cin >> te;
	}

	a = (xa[khi] - x) / h;
	b = (x - xa[klo]) / h;
	*y = a * ya[klo] + b * ya[khi] + ((a * a * a - a) * y2a[klo] + (b * b * b
	- b) * y2a[khi]) * (h * h) / 6.0;
}

//********************************************************************************

double splint(double xa[], double ya[], double y2a[], int n, double x)
// spline interpolation from numerical recipes in C, same as aboce, but returns the function value
{
	void nrerror(char error_text[]);
	int klo, khi, k;
	double h, b, a;

	klo = 0;
	khi = n - 1;
	while (khi - klo > 1) {
		k = (khi + klo) >> 1;
		if (xa[k] > x) {
			khi = k;
		} else {
			klo = k;
		}
	}

	h = xa[khi] - xa[klo];
	if (h == 0.0) {
		printf(
		"\n\nBad xa input to routine splint khi=%d  klo=%d  xlo=%g  xhi=%g  n=%d \n\n",
		khi, klo, xa[klo], xa[khi], n);
		int te;
		cin >> te;
	}

	a = (xa[khi] - x) / h;
	b = (x - xa[klo]) / h;
	return a * ya[klo] + b * ya[khi] + ((a * a * a - a) * y2a[klo] + (b * b * b
	- b) * y2a[khi]) * (h * h) / 6.0;
}

//********************************************************************************


void splint_deriv(double xa[], double ya[], double y2a[], int n, double x, double *yd)
// spline interpolation from numerical recipes in C, given fct ya(xa) and having calculated y2a with spline-fct, # pts=n, this writes the interpolated
{
	void nrerror(char error_text[]);
	int klo, khi, k;
	double h, b, a;

	klo = 0;
	khi = n - 1;
	while (khi - klo > 1) {
		k = (khi + klo) >> 1;
		if (xa[k] > x) {
			khi = k;
		} else {
			klo = k;
		}
	}

	h = xa[khi] - xa[klo];
	if (h == 0.0) {
		printf("\n\nBad xa input to routine splint\n\n");

		int te;
		cin >> te;
	}

	a = (xa[khi] - x) / h;
	b = (x - xa[klo]) / h;
	*yd = -ya[klo] / h + ya[khi] / h + ((-3.0 * a * a / h + 1 / h) * y2a[klo]
	+ (3.0 * b * b / h - 1 / h) * y2a[khi]) * (h * h) / 6.0;

}

//********************************************************************************
void spline(double x[], double y[], int n, double yp1, double ypn, double y2[])

// -----------------------------------------------------------------------------
//    spline tabulation method from numerical recipes in C, given a fct y(x) with
//    n tabulated pts, x_0<...<x_{n-1},
// -----------------------------------------------------------------------------
// yp1=slope at rist pt, ypn=slope at last pt, calculates 2nd derivative y2
{
	int i, k;
	double p, qn, sig, un;
	double *u = new double[n];

	if (yp1 > 0.99e30) {
		y2[0] = u[0] = 0.0;
	} else {
		y2[0] = -0.5;
		u[0] = (3.0 / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1);
	}

	for (i = 1; i <= n - 2; i++) {
		sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
		p = sig * y2[i - 1] + 2.0;
		y2[i] = (sig - 1.0) / p;
		u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1])
		/ (x[i] - x[i - 1]);
		u[i] = (6.0 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
	}

	if (ypn > 0.99e30) {
		qn = un = 0.0;
	} else {
		qn = 0.5;
		un = (3.0 / (x[n - 1] - x[n - 2])) * (ypn - (y[n - 1] - y[n - 2])
		/ (x[n - 1] - x[n - 2]));
	}

	y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0);
	for (k = n - 2; k >= 1; k--) {
		y2[k] = y2[k] * y2[k + 1] + u[k];
	}

	delete[] u;
}
//********************************************************************************

double calc_J_analytic_approx(double V_lat) {
	return 4 / sqrt(PI) * pow(V_lat, 0.75) * exp(-2 * sqrt(V_lat));
}
/* ----------------------------------------------------------------------------------------------------*/

double calc_U_analytic_approx(double V_lat, double qlat_times_as) {
	return sqrt(8 / PI) * qlat_times_as * pow(V_lat, 0.75);
}
/* ----------------------------------------------------------------------------------------------------*/

//Ridder's Method for finding roots of a function, here: a slight generalization to find the specific function value fct_val
#define MAXIT 60
#define UNUSED (-1.11e30)

double Ridder_find_val(double(*func)(double), double func_val, double x1, double x2, double x_accuracy) {
	//Ridder's algorithm fitting an exponential through three points, iteratively determining and rteturning the value of x1<x<x2 , such that func(x)=func_val
	//use
	int j;
	double ans, fh, fl, fm, fnew, s, xh, xl, xm, xnew;

	fl = (*func)(x1) - func_val;
	fh = (*func)(x2) - func_val;
	if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
		xl = x1;
		xh = x2;
		ans = UNUSED;
		for (j = 1; j <= MAXIT; j++) {
			//printf("run=%d\n", j);
			xm = 0.5 * (xl + xh);
			fm = (*func)(xm) - func_val;
			s = sqrt(fm * fm - fl * fh);
			if (s == 0.0)
			return ans;
			xnew = xm + (xm - xl) * ((fl >= fh ? 1.0 : -1.0) * fm / s);
			if (fabs(xnew - ans) <= x_accuracy)
			return ans;
			ans = xnew;
			fnew = (*func)(ans) - func_val;
			if (fnew == 0.0)
			return ans;
			if (SIGN(fm,fnew) != fm) {
				xl = xm;
				fl = fm;
				xh = ans;
				fh = fnew;
			} else if (SIGN(fl,fnew) != fl) {
				xh = ans;
				fh = fnew;
			} else if (SIGN(fh,fnew) != fh) {
				xl = ans;
				fl = fnew;
			} else
			abort(); // never get here
			if (fabs(xh - xl) <= x_accuracy)
			return ans;
		}
		abort(); // error: "zriddr exceed maximum iterations"
	} else {
		if (fl == 0.0)
		return x1;
		if (fh == 0.0)
		return x2;
		abort(); // never get here: root must be bracketed in zriddr
	}
	return 0.0;
}
#undef MAXIT
#undef UNUSED


/* ----------------------------------------------------------------------------------------------------*/


//Ridder's Method for finding roots of a function, here: a slight generalization to find the specific function value fct_val
#define MAXIT 60
#define UNUSED (-1.11e30)

double Ridder_find_val(double(*func)(double, parameters &), double func_val, double x1, double x2, double x_accuracy, parameters &A) {
	//Ridder's algorithm fitting an exponential through three points, iteratively determining and rteturning the value of x1<x<x2 , such that func(x)=func_val

	int j;
	double ans, fh, fl, fm, fnew, s, xh, xl, xm, xnew;

	fl = (*func)(x1, A) - func_val;
	fh = (*func)(x2, A) - func_val;
	if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
		xl = x1;
		xh = x2;
		ans = UNUSED;
		for (j = 1; j <= MAXIT; j++) {
			//printf("run=%d\n", j);
			xm = 0.5 * (xl + xh);
			fm = (*func)(xm, A) - func_val;
			s = sqrt(fm * fm - fl * fh);
			if (s == 0.0)
			return ans;
			xnew = xm + (xm - xl) * ((fl >= fh ? 1.0 : -1.0) * fm / s);
			if (fabs(xnew - ans) <= x_accuracy)
			return ans;
			ans = xnew;
			fnew = (*func)(ans, A) - func_val;
			if (fnew == 0.0)
			return ans;
			if (SIGN(fm,fnew) != fm) {
				xl = xm;
				fl = fm;
				xh = ans;
				fh = fnew;
			} else if (SIGN(fl,fnew) != fl) {
				xh = ans;
				fh = fnew;
			} else if (SIGN(fh,fnew) != fh) {
				xl = ans;
				fl = fnew;
			} else
			abort(); // never get here
			if (fabs(xh - xl) <= x_accuracy)
			return ans;
		}
		abort(); // error: "zriddr exceed maximum iterations"
	} else {
		if (fl == 0.0)
		return x1;
		if (fh == 0.0)
		return x2;
		abort(); // never get here: root must be bracketed in zriddr
	}
	return 0.0;
}
#undef MAXIT
#undef UNUSED


/* ----------------------------------------------------------------------------------------------------*/
//applies the annihilation operator to a state in the Fock representation, altering the input vector
//the state is of length N, i.e. occupation numbers from 0 to N-1
//the last amplitude is set to zero, as the space is truncated
//normalizes the output state and returns the norm. i.e. b|psi>= a*|phi> with <phi|phi>=1, then a is returned

void apply_b_Operator_to_state(cd state[], int N, double &norm) {

	//for(int m=0; m<N; m++) printf("state[%d]=%f+i%f   norm=%f\n",m,real(state[m]),imag(state[m]),abs(state[m]*state[m]));
	double norm_sq = 0.0;
	for (int m = 1; m < N; m++) {
		state[m - 1] = sqrt(double(m)) * state[m];
		norm_sq += abs(state[m - 1] * state[m - 1]);
	}
	state[N - 1] = 0.0;
	norm = sqrt(norm_sq);
	for (int m = 0; m < N; m++) {
		state[m] = state[m] / norm;
	}
}

/* ----------------------------------------------------------------------------------------------------*/
//applies the creation operator to a state in the Fock representation, altering the input vector
//the state is of length N, i.e. occupation numbers from 0 to N-1
//the first amplitude is set to zero
//normalizes the output state and returns the norm. i.e. b'|psi>= a*|phi> with <phi|phi>=1, then a is returned

void apply_b_dag_Operator_to_state(cd state[], int N, double &norm) {

	printf("\n\nbefore:\n");
	for (int m = 0; m < N; m++)
	printf("%1.14f\n",real(state[m]));

	//for(int m=0; m<N; m++) printf("state[%d]=%f+i%f   norm=%f\n",m,real(state[m]),imag(state[m]),abs(state[m]*state[m]));
	double norm_sq = 0.0;
	for (int m = N-1; m >0; m--) {
		state[m] = sqrt(double(m)) * state[m-1];
		norm_sq += abs(state[m] * state[m]);
	}
	state[0] = 0.0;
	norm = sqrt(norm_sq);

	printf("\n\nafter:\n");
	for (int m = 0; m < N; m++)
	printf("%1.14f\n",real(state[m]));


	for (int m = 0; m < N; m++) {
		state[m] = state[m] / norm;
	}
}

/* ----------------------------------------------------------------------------------------------------*/
//applies the creation operator to a state in the Fock representation, altering the input vector
//the state is of length N, i.e. occupation numbers from 0 to N-1
//the first amplitude is set to zero
//normalizes the output state and returns the norm. i.e. b'|psi>= a*|phi> with <phi|phi>=1, then a is returned

void apply_n_Operator_to_state(cd state[], int N, double &norm) {

	double norm_sq = 0.0;
	for (int m = 0; m < N; m++) {
		state[m] = double(m) * state[m];
		norm_sq += abs(state[m] * state[m]);
	}
	norm = sqrt(norm_sq);

	for (int m = 0; m < N; m++) {
		state[m] = state[m] / norm;
	}
}





/* ----------------------------------------------------------------------------------------------------*/
//random number generator ran1 from numerical recipes in c
//good for up to 100 000 000 random numbers, thereafter use extended version
//returns uniformly distributed numbers between 0 and 1
//initialize with negative idum, do not change within sequence
double num_recipes_ran1(long &idum)
{
	const long IA=16807,IM=2147483647,IQ=127773,IR=2836,NTAB=32;
	const long NDIV=(1+(IM-1)/NTAB);
	const double EPS=3.0e-16,AM=1.0/IM,RNMX=(1.0-EPS);
	static long iy=0;
	static long iv[NTAB];
	int j;
	long k;
	double temp;

	if (idum <= 0 || !iy) {
		if (-idum < 1) idum=1;
		else idum = -idum;
		for (j=NTAB+7;j>=0;j--) {
			k=idum/IQ;
			idum=IA*(idum-k*IQ)-IR*k;
			if (idum < 0) idum += IM;
			if (j < NTAB) iv[j] = idum;
		}
		iy=iv[0];
	}
	k=idum/IQ;
	idum=IA*(idum-k*IQ)-IR*k;
	if (idum < 0) idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}






/* ----------------------------------------------------------------------------------------------------*/
// function to read a matrix from a text file to the double array A
// the m times n matrix A is initialized in this function, i.e. only the pointer has to be passed!
// The matrix dimensions are also written into the integer values m and n
// the matrix in the text file should have the following form
/*
sample_matrix.txt:
****************************
3 by 4 matrix

1.0 2 3 4
4 5 67 5
654 5 5 6
***************************/

void read_matrix_from_file(const std::string &fileName, double**&A, int &m, int &n)
{
	FILE *matrix_file=fopen(fileName.c_str(),"r");

	int i,j,fileExtractions;
	char str [80];
	fileExtractions=fscanf(matrix_file,"%i %s %i",&m,str,&n);	//just to avoid unused warning when compiling with O3 optimization
    if(fileExtractions!=3)
        printf("ERROR in __FILE__ __LINE__");
	fileExtractions=fscanf(matrix_file,"%s",str);
    if(fileExtractions!=1)
        printf("ERROR in __FILE__ __LINE__");
	A = new double* [m];
	for(i=0;i<m;i++)
	{
		A[i] = new double [n];
		for(j=0;j<n;j++)
		{
			fileExtractions=fscanf(matrix_file,"%lg",&A[i][j]);
            if(fileExtractions!=1)
                printf("ERROR in __FILE__ __LINE__");
		}
	}

	fclose(matrix_file);
}




/* ----------------------------------------------------------------------------------------------------*/
/* ----------------------------------------------------------------------------------------------------*/
/* ----------------------------------------------------------------------------------------------------*/
tensor4::tensor4(int dim_tmp){
	//assign 4D array dynamically and initialize to 0
	dim=dim_tmp;
	A=new double*** [dim];
	for(int ct1=0; ct1<dim; ct1++){
		A[ct1]=new double** [dim];
		for(int ct2=0; ct2<dim; ct2++){
			A[ct1][ct2]=new double* [dim];
			for(int ct3=0; ct3<dim; ct3++){
				A[ct1][ct2][ct3]=new double[dim];
				for(int ct4=0; ct4<dim; ct4++){
					A[ct1][ct2][ct3][ct4]=0.0;
				}
			}
		}
	}
}

tensor4::tensor4(const std::string &fileName){
	//constructor reading all data from a file: the file contains a m x 5 matrix, the vlaues are in the 5th column, the prior ones the indices
	double **tmp_A;
	int m_tmp,n_tmp;
	read_matrix_from_file(fileName, tmp_A, m_tmp, n_tmp);
	if(n_tmp!=5){printf("error: matrix read in by tensor4 constructor does not contain 5 columns. Exiting.");abort();}

	//find largest index:
	dim=0;
	for(int ct1=0; ct1<4; ct1++){
	for(int ct2=0; ct2<m_tmp; ct2++){
	if( int(tmp_A[ct2][ct1])>dim ) dim=int(tmp_A[ct2][ct1]);
	}
	}
	dim++; //length is one larger than largest index in C++ notation

	//create structure and initialize all with 0
	A=new double*** [dim];
	for(int ct1=0; ct1<dim; ct1++){
		A[ct1]=new double** [dim];
		for(int ct2=0; ct2<dim; ct2++){
			A[ct1][ct2]=new double* [dim];
			for(int ct3=0; ct3<dim; ct3++){
				A[ct1][ct2][ct3]=new double[dim];
				for(int ct4=0; ct4<dim; ct4++){
					A[ct1][ct2][ct3][ct4]=0.0;
				}
			}
		}
	}

	for(int ct=0; ct<m_tmp; ct++){
	A[int(tmp_A[ct][0])][int(tmp_A[ct][1])][int(tmp_A[ct][2])][int(tmp_A[ct][3])]=tmp_A[ct][4];
	}

	//delete tmp_A
	for(int ct=0; ct<m_tmp; ct++) delete[] tmp_A[ct];
	delete[] tmp_A;
}


tensor4::tensor4(const tensor4& tensor_to_be_copied){ //copy constructor, needed because array a is dynamically allocated
	/*
	dim=tensor_to_be_copied.dim;
	tensor4(dim);	//use normal constructor to initialize, this however fucks up the following part
	for(int ct1=0; ct1<dim; ct1++){			//why can I not use ct<"dim"	 here??????
		for(int ct2=0; ct2<tensor_to_be_copied.dim; ct2++){
			for(int ct3=0; ct3<tensor_to_be_copied.dim; ct3++){
				for(int ct4=0; ct4<tensor_to_be_copied.dim; ct4++){
					A[ct1][ct2][ct3][ct4]=tensor_to_be_copied.A[ct1][ct2][ct3][ct4];	//why does this crash?
				}
			}
		}
	}
	*/
	dim=tensor_to_be_copied.dim;
	A=new double*** [dim];
	for(int ct1=0; ct1<dim; ct1++){
		A[ct1]=new double** [dim];
		for(int ct2=0; ct2<dim; ct2++){
			A[ct1][ct2]=new double* [dim];
			for(int ct3=0; ct3<dim; ct3++){
				A[ct1][ct2][ct3]=new double[dim];
				for(int ct4=0; ct4<dim; ct4++){
					A[ct1][ct2][ct3][ct4]=tensor_to_be_copied.A[ct1][ct2][ct3][ct4];
				}
			}
		}
	}



	}


void tensor4::operator= (const tensor4& tensor_to_be_copied){ //similar to copy constructor, needed because array a is dynamically allocated

for(int ct1=0; ct1<dim; ct1++){
		for(int ct2=0; ct2<dim; ct2++){
			for(int ct3=0; ct3<dim; ct3++){
				delete[] A[ct1][ct2][ct3];
			}
			delete[] A[ct1][ct2];
		}
		delete[] A[ct1];
	}
	delete[] A;


	dim=tensor_to_be_copied.dim;
	A=new double*** [dim];
	for(int ct1=0; ct1<dim; ct1++){
		A[ct1]=new double** [dim];
		for(int ct2=0; ct2<dim; ct2++){
			A[ct1][ct2]=new double* [dim];
			for(int ct3=0; ct3<dim; ct3++){
				A[ct1][ct2][ct3]=new double[dim];
				for(int ct4=0; ct4<dim; ct4++){
					A[ct1][ct2][ct3][ct4]=tensor_to_be_copied.A[ct1][ct2][ct3][ct4];
				}
			}
		}
	}
}

tensor4::~tensor4(){
	for(int ct1=0; ct1<dim; ct1++){
		for(int ct2=0; ct2<dim; ct2++){
			for(int ct3=0; ct3<dim; ct3++){
				delete[] A[ct1][ct2][ct3];
			}
			delete[] A[ct1][ct2];
		}
		delete[] A[ct1];
	}
	delete[] A;
}

void tensor4::disp(){ //similar to copy constructor, needed because array a is dynamically allocated
printf("\n******************************\ndim=%d\nthis tensor=\n",dim);
for(int ct1=0; ct1<dim; ct1++){
		for(int ct2=0; ct2<dim; ct2++){
			for(int ct3=0; ct3<dim; ct3++){
				for(int ct4=0; ct4<dim; ct4++){
				printf("A[%d][%d][%d][%d]=%f\n",ct1,ct2,ct3,ct4,A[ct1][ct2][ct3][ct4]);
				}
			}
		}
	}
printf("******************************\n");
}

double tensor4::getval(int i1,int i2,int i3,int i4){
	if(i1>=0 && i1<dim && i2>=0 && i2<dim && i3>=0 && i3<dim && i4>=0 && i4<dim){
	return A[i1][i2][i3][i4];
	}else{
	printf("error:index out of range in tensor4:getval indices:(%d %d %d %d)\n",i1,i2,i3,i4);return 1E30;
	}
}

void tensor4::setval(int i1,int i2,int i3,int i4,double val){
if(i1>=0 && i1<dim && i2>=0 && i2<dim && i3>=0 && i3<dim && i4>=0 && i4<dim){
	A[i1][i2][i3][i4]=val;
	}else{
	printf("error:index out of range in tensor4:setval\n");
	}
}
/* ----------------------------------------------------------------------------------------------------*/
/* ----------------------------------------------------------------------------------------------------*/

///the factorial function
int factorial(int x) {
int fac = 1;
for (int i=2; i<=x; i++) fac *= i;
return fac;
}
/* ----------------------------------------------------------------------------------------------------*/
/* ----------------------------------------------------------------------------------------------------*/


double associated_legendre_pol(const int l, const int m, const double x)
///Computes the associated Legendre polynomial P_l^m(x) . Here m and l are integers satisfying 0 <= m <= l , while x lies in the range -1 <= x <= 1 .
///taken from the numerical Recipes in C routine "plgndr" in "plgndr.cpp"
{
	int i,ll;
	double fact,pmm,pmmp1,somx2;
    double pll=0;
	if (m < 0 || m > l || fabs(x) > 1.0){
		printf("Bad arguments in routine plgndr");
		return 0;
	}
	pmm=1.0;
	if (m > 0) {
		somx2=sqrt((1.0-x)*(1.0+x));
		fact=1.0;
		for (i=1;i<=m;i++) {
			pmm *= -fact*somx2;
			fact += 2.0;
		}
	}
	if (l == m)
		return pmm;
	else {
		pmmp1=x*(2*m+1)*pmm;
		if (l == (m+1))
			return pmmp1;
		else {
			for (ll=m+2;ll<=l;ll++) {
				pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
				pmm=pmmp1;
				pmmp1=pll;
			}
			return pll;
		}
	}
}


/* ----------------------------------------------------------------------------------------------------*/
/* ----------------------------------------------------------------------------------------------------*/



cd spherical_harmonic(int l, int m, double theta, double phi){
    return sqrt((2.0*l+1.0)*factorial(l-m)  /  (4.0*PI*factorial(l+m) ))* associated_legendre_pol(l, m, cos(theta)) *cd(cos(m*theta),sin(m*theta));
}


/* ----------------------------------------------------------------------------------------------------*/
/* ----------------------------------------------------------------------------------------------------*/


///these functions are independent of phi and purely real
double spherical_harmonic_m0(int l, double theta){
    return sqrt((2.0*double(l)+1.0) /  (4.0*PI))* associated_legendre_pol(l, 0, cos(theta));
}


/* ----------------------------------------------------------------------------------------------------*/
/* ----------------------------------------------------------------------------------------------------*/
#ifndef polint_numerical_recipes_defined
#define polint_numerical_recipes_defined 1
///polynomial interpolation from numerical recipes, taken from file polint.c
template<typename domainType, typename rangeType>
void polint(domainType xa[], rangeType ya[], int n, domainType x, rangeType *y, rangeType *dy)
{
    int i,m,ns=1;
    domainType dif, dift, ho, hp; 
    rangeType w, den;
    rangeType* c=(rangeType*) malloc((n+1)*sizeof(rangeType));
    rangeType* d=(rangeType*) malloc((n+1)*sizeof(rangeType));

    dif=fabs(x-xa[1]);
    for (i=1;i<=n;i++) {
        if ( (dift=fabs(x-xa[i])) < dif) {
            ns=i;
            dif=dift;
        }
        c[i]=ya[i];
        d[i]=ya[i];
    }
    *y=ya[ns--];
    for (m=1;m<n;m++) {
        for (i=1;i<=n-m;i++) {
            ho=xa[i]-x;
            hp=xa[i+m]-x;
            w=c[i+1]-d[i];
            if ( (den=ho-hp) == 0.0) {
                printf("Error in routine polint");
                abort();
            }
            den=w/den;
            d[i]=hp*den;
            c[i]=ho*den;
        }
        *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
    }
    free(c);
    free(d);
}
#endif

/* ----------------------------------------------------------------------------------------------------*/
/* ----------------------------------------------------------------------------------------------------*/
/// trapezoidal integration (nth stage of the refinement) routine taken from numerical recipes in C++, the file trapzd.c
/// n=1 is the crudest approximation. Subsequent calls with n=2,3,... (in that sequential order) will improve the accuracy by adding 2^{n-2} additional interior points

double trapez_int_refined(double (*func)(double), double x_lo, double x_hi, int n, double &s)
{
	double x,tnm,sum,del;
	int it,j;
	if (n == 1) {
		return (s=0.5*(x_hi-x_lo)*( ((*func)(x_lo)) +((*func)(x_hi))  ));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(x_hi-x_lo)/tnm;
		x=x_lo+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) sum += ((*func)(x));
		s=0.5*(s+(x_hi-x_lo)*sum/tnm);
		return s;
	}
}


/* ----------------------------------------------------------------------------------------------------*/
/* ----------------------------------------------------------------------------------------------------*/
/// the same as above, but with the option to pass an additional parameter object
/// trapezoidal integration (nth stage of the refinement) routine taken from numerical recipes in C++, the file trapzd.cpp
/// n=1 is the crudest approximation. Subsequent calls with n=2,3,... (in that sequential order) will improve the accuracy by adding 2^{n-2} additional interior points
double trapez_int_refined(double(*func)(double, parameters &), const double x_lo, const double x_hi, const int n, double &s, parameters &P)
{
	double x,tnm,sum,del;
	int it,j;
	if (n == 1) {
		return (s=0.5*(x_hi-x_lo)*(func(x_lo, P)+func(x_hi, P)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(x_hi-x_lo)/tnm;
		x=x_lo+0.5*del;
		for (sum=0.0,j=0;j<it;j++,x+=del) sum += func(x, P);
		s=0.5*(s+(x_hi-x_lo)*sum/tnm);
		return s;
	}
}




/* ----------------------------------------------------------------------------------------------------*/
/* ----------------------------------------------------------------------------------------------------*/
#ifndef romberg_integration_precision
#define romberg_integration_precision 1E-9
#endif
/// Romberg integration routine taken from numerical recipes in C++, the file qromb.c
double romberg_integral(double (*func)(double), double x_lo, double x_hi)
{
    double s_in_trapez_int_refined; ///initially declared as a static variable in trapez_int_refined - pass this ba reference
	const int JMAX=30, JMAXP=JMAX+1;
	const int K=6;      /// determines quality: if too low the method seems to sometime have thought that it has converged, where it hasn't
	const double EPS=romberg_integration_precision;
	//printf("romberg_integration_precision=%g\n",EPS);
	double ss,dss;
	double s[JMAXP+1],h[JMAXP+1];
	int j;
	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=trapez_int_refined(func,x_lo,x_hi,j, s_in_trapez_int_refined);
		if (j >= K) {
			polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) < EPS*fabs(ss) ||  fabs(ss) < 1E-9 ) return ss;
		}
		s[j+1]=s[j];
		h[j+1]=0.25*h[j];
/* 		printf ("qromb: ss = %le, dss = %le, EPS = %le, j=%d\n",ss,dss,EPS,j); */
	}
	printf("Error: Too many steps in routine qromb\n\n");
	return 0.0;
}


/* ----------------------------------------------------------------------------------------------------*/
/* ----------------------------------------------------------------------------------------------------*/
/// the same as above, but with the option to pass an additional parameter object
/// Romberg integration routine taken from numerical recipes in C++, the file qromb.c
double romberg_integral(double(*func)(double, parameters &), double x_lo, double x_hi, parameters &P)
{
    double s_in_trapez_int_refined; ///initially declared as a static variable in trapez_int_refined - pass this ba reference
	const int JMAX=30, JMAXP=JMAX+1;
	const int K=6;      /// determines quality: if too low the method seems to sometime have thought that it has converged, where it hasn't
	const double EPS=romberg_integration_precision;
	//printf("romberg_integration_precision=%g\n",EPS);
	double ss,dss;
	double s[JMAXP+1],h[JMAXP+1];
	int j;
	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=trapez_int_refined(func,x_lo,x_hi,j, s_in_trapez_int_refined,P);
		if (j >= K) {
			polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) < EPS*fabs(ss) ||  fabs(ss) < 1E-9 ) return ss;
		}
		s[j+1]=s[j];
		h[j+1]=0.25*h[j];
/* 		printf ("qromb: ss = %le, dss = %le, EPS = %le, j=%d\n",ss,dss,EPS,j); */
	}

	printf("Error: Too many steps in routine qromb\n\n");
	return 0.0;
}

#endif // _UTILS_CPP_
