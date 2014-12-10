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

const uint MAX = 128;
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

scalar bn(uint n) {
	scalar sum = 0;
	for (uint k = 1; k <= n + 1; k++) {
		scalar s = 0;
		for (uint j = 1; j <= k; j++) s += pow(j, n);
		sum += s * binom(n + 1, k) * PM1(k) / k;
	}
	return -sum;
}

// http://math.ucsb.edu/~jcs/Bernoulli.pdf
void calcBn() {
	for (uint n = 0; n < MAX - 1; ++n) Bn[n] = bn(n);
	return;
	Bn[0] = 0;
	Bn[1] = -half;
	for (uint n = 2; n < MAX - 1; ++n) {
		if (n % 2) 
			Bn[n] = 0; 
		else {
			scalar sum = 0;
			for (uint k = 0; k < n; k++) 
				sum += binom(n + 1, k) * Bn[k];
			Bn[n] = -sum / scalar(n + 1);
		}
	}
}

scalar mus[MAX];
scalar kappas[MAX];

void precalc() { calcFactorials(); calcBn(); for(uint n = 0; n < MAX; n++) mus[n] = kappas[n] = 0; }

scalar kappa(const uint m) {
	if (kappas[m] != 0) return kappas[m];
	scalar sum = 0, m4 = pow(four, m*2);
	for (uint n = 0; n < N; n++) sum += pow(x[n], 2 * n);
	kappas[m] = sum * Bn[m*2] * m4 * (m4 - one);
	cout<<"kappa["<<m<<"]\t= "<<kappas[m]<<endl;
	return kappas[m];
}


scalar mu(const uint m) {
	if (mus[m] != 0) return mus[m];
	scalar sum = kappa(m);

	for (uint t = 2; t < m - 1; t += 2)
		sum += binom(m - 1, t - 1) * kappa(t) * mu(m - t);

	cout<<"mu["<<m<<"]\t= "<<sum<<endl;

	return mus[m] = sum;
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

//	for (uint n = 0; n < MAX; n++) cout<<"B["<<n<<"]\t= " << Bn[n] << endl;

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











