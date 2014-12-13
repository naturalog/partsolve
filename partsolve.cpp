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
/*		y = psi(x) - c;
		t = res + y;
		c = (t - res) - y;
		res = t;*/
		res += (/*samples[k] =*/ psi(x));
	}
	return res * dx * pi;
}

int main(int, char**) {
#ifdef MPREC
	mpreal::set_default_prec(prec);
#endif

	string str;
	getline(cin, str, ',');
	x = new scalar[N = atoi(str.c_str())];
	for (uint k = 0; k < N; k++) {
		getline(cin, str, ',');
		cout << (x[k] = atoi(str.c_str())) << ',';
	}
	cout<<endl;
	scalar t = 0;
	size = 1;
	while (size *= 10) cout<<numint()<<endl;
//	precalc();
//	for (uint n = 0; n < MAX; n++) cout<<"B["<<n<<"]\t= " << Bn[n] << endl;
//	scalar e = estimator();
//	cout << endl << "estimator: " << e << endl;

	return 0;
}
