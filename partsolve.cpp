#include <iostream>
#include <cmath>
#include <string>
#include <cstdlib>
using namespace std;

#define MPREC

#ifdef MPREC
#include <mpreal.h>
const uint prec = 512;
using namespace mpfr;
typedef mpreal scalar;
#endif

const uint MAX = 1024;
const scalar one = 1;
const scalar two = 2;
const scalar half = one / two;
const scalar four = 4;
#define PM1(x) ( ( x % 2 ) ? -one : one)

scalar Bn[MAX];
scalar facts[MAX];
scalar *x;
uint N;

void calcFactorials() {
	facts[0] = 1;
	for (uint n = 1; n < MAX; ++n) facts[n] = n * facts[n - 1];
}

inline scalar binom(uint n, uint k) { return facts[n] / (facts[k] * facts[n - k]); }

// http://arxiv.org/pdf/1108.0286v3.pdf
void calcBn() {
	Bn[0] = 0;
	Bn[1] = 1;
	for (uint n = 2; n < MAX - 1; ++n) {
		scalar sum = 0;
		for (uint k = 0; k < n; k++) 
			sum += binom(n + 1, k) * Bn[k];
		Bn[n] = -sum / scalar(n + 1);
	}
}

void precalc() { calcFactorials(); calcBn(); }

scalar kappa(const uint m) {
	scalar sum = 0, m4 = pow(four, m);
	for (uint n = 0; n < N; n++) sum += pow(x[n], 2 * n);
	return sum * Bn[m] * m4 * (m4 - one);
}

scalar mu(const uint m) {
	scalar sum = kappa(m);

	for (uint t = 2; t < m - 1; t += 2)
		sum += binom(m - 1, t - 1) * kappa(t) * mu(m - t);

	return sum;
}

scalar estimator(const uint n = N) {
	scalar sum = 1, rem = 1;
	scalar eps = pow(half, n + 1);
	uint m = 2;
	while (rem > eps) {
		sum += (rem = PM1(m) * mu(2 * m) / facts[2 * m + 1]);
		m += 2;
		cout << "m: " << m << "\trem:\t" << rem << "\tsum:\t" << sum << endl;
	}
	return sum;
}

int main(int, char**) {
#ifdef MPREC
	mpreal::set_default_prec(prec);
#endif
	precalc();

	for (uint n = 0; n < MAX; n++) cout<<"B["<<n<<"]\t=" << Bn[n] << endl;

	string str;
	getline(cin, str);
	x = new scalar[N = atoi(str.c_str())];
	for (uint k = 0; k < N; k++) {
		getline(cin, str);
		cout << (x[k] = atoi(str.c_str())) << ',';
	}
	cout<<endl;
	scalar e = estimator();
	cout << endl << "estimator: " << e << endl;


	return 0;
}











