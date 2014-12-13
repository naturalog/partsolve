// g++ ./partsolve.cpp  -O3 -I/opt/AMDAPP/include -std=c++11 -lOpenCL -lboost_system
#include <iostream>
#include <cmath>
#include <string>
#include <cstdlib>
//#include <Eigen/Core>
//#include <unsupported/Eigen/MPRealSupport>
#include <vexcl/vexcl.hpp>
using namespace std;
//using namespace Eigen;

//#define MPREC

#ifdef MPREC
#include <mpreal.h>
const uint prec = 1024; 
using namespace mpfr;
typedef mpreal scalar;
#else
typedef double scalar;
#endif

//typedef Matrix<scalar, Dynamic, Dynamic> mat;
const scalar one = 1;
const scalar two = 2;
const scalar half = one / two;
const scalar four = 4;
const scalar h = 1e-15;

scalar *x;
uint N;
/*
scalar psi(const scalar& t, uint d = 0) {
//	if (d) return (psi(t+h,d-1)-psi(t-h,d-1))/(h*two);
	static const scalar pi = acos(-one);
//	return cos(pi*(one+two)*t);
	scalar pos = 0, neg = 0, z;
	for (uint n = 0; n < N; n++) 
		( (z = cos(t * x[n])) > 0 ) ? 
			( pos += log(z) : neg += log(-z) );
	
	return exp(pos-neg);
}

scalar* samples = 0;
uint size = 1e+6;

scalar numint() {
	static const scalar pi = acos(-one);
	const scalar dx = one/scalar(size);
	scalar res = 0, c = 0, y, t, x;
	cout<<size<<' ';
	if (samples) delete[] samples;
	samples = new scalar[size + 1];
	const scalar prsz = pi / scalar(size);
	for (uint k = 0; k <= size; ++k) {
		x = scalar(k) * prsz;
//		y = psi(x) - c;
//		t = res + y;
//		c = (t - res) - y;
//		res = t;
		res += (samples[k] = psi(x));
	}
	return res * dx * pi;
}
*/
int main(int argc, char** argv) {
#ifdef MPREC
	mpreal::set_default_prec(prec);
#endif
	if (argc != 2) { cout<<"usage: <iters>"<<endl; return 1; }
	vex::Context ctx( vex::Filter::Type(CL_DEVICE_TYPE_GPU) && vex::Filter::DoublePrecision );
	cout << ctx << endl;


	string str;
	getline(cin, str, ',');
	x = new scalar[N = atoi(str.c_str())];
	for (uint k = 0; k < N; k++) {
		getline(cin, str, ',');
		cout << (x[k] = atoi(str.c_str())) << ',';
	}
	cout<<endl;

	scalar result = 0;
	uint max_iters = 10;
	double size = atoi(argv[1]);
	cout<<"iters: "<<size<<endl;
	for (uint iters = 0; iters <= max_iters; iters++) {
		double dx = double(1)/(double(size-1) * max_iters);
		vex::vector<double> X(ctx, size), Y(ctx, size);
		Y = 1;
		vex::Reductor<double, vex::SUM> sum(ctx);

		X = (vex::constants::pi() * (vex::element_index() + double(iters)*size)) * dx;
		for (uint k = 0; k < N; k++)
			if (N - k > 4) {
				Y *= cos(X * x[k]) * cos(X * x[k + 1]) * cos(X * x[k + 2]);
				k += 3;
			} else Y *= cos(X * x[k]);
		Y *= vex::constants::pi() * dx;
		cout<<iters<<' '<<(result += sum(Y))<<endl;
	}
	
//	scalar t = 0;
//	size = 1;
//	while (size *= 10) cout<<numint()<<endl;
//	precalc();
//	for (uint n = 0; n < MAX; n++) cout<<"B["<<n<<"]\t= " << Bn[n] << endl;
//	scalar e = estimator();
//	cout << endl << "estimator: " << e << endl;

	return 0;
}
