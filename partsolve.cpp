#include <iostream>
#include <cmath>
#include <string>
#include <cstdlib>
#include <Eigen/Core>
#include <unsupported/Eigen/MPRealSupport>
using namespace std;
using namespace Eigen;

//#define MPREC

#ifdef MPREC
#include <mpreal.h>
const uint prec = 1024; 
using namespace mpfr;
typedef mpreal scalar;
#else
typedef long double scalar;
#endif

typedef Matrix<scalar, Dynamic, Dynamic> mat;
const scalar one = 1;
const scalar two = 2;
const scalar half = one / two;
const scalar four = 4;
const scalar h = 1e-15;

scalar *x;
uint N;

scalar psi(const scalar& t/*, uint d = 0*/) {
//	if (d) return (psi(t+h,d-1)-psi(t-h,d-1))/(h*two);
	static const scalar pi = acos(-one);
//	return cos(pi*(one+two)*t);
	scalar res = 1;
	for (uint n = 0; n < N; n++) res *= cos(t*x[n]);
	return res;
}

scalar numint(const scalar w) {
	static const scalar pi = acos(-one);
	const scalar dx = one/w;
	scalar res = 0, c = 0, y, t, x;
	cout<<w<<' ';
	for (scalar k = -w; k <= w; ++k) {
		x = k*pi/w;
/*		y = psi(x) - c;
		t = res + y;
		c = (t - res) - y;
		res = t;*/
		res += psi(x);
	}
	return res * dx * pi;
}

int main(int, char**) {
#ifdef MPREC
	mpreal::set_default_prec(prec);
#endif

	string str;
	getline(cin, str);
	x = new scalar[N = atoi(str.c_str())];
	for (uint k = 0; k < N; k++) {
		getline(cin, str);
		cout << (x[k] = atoi(str.c_str())) << ',';
	}
	cout<<endl;
	scalar t = 0;
	uint n = 1;
	while (n *= 2) cout<<numint(n)<<endl;
//	precalc();
//	for (uint n = 0; n < MAX; n++) cout<<"B["<<n<<"]\t= " << Bn[n] << endl;
//	scalar e = estimator();
//	cout << endl << "estimator: " << e << endl;

	return 0;
}
