// g++ ./partsolve.cpp  -O3 -I/opt/AMDAPP/include -std=c++11 -lOpenCL -lboost_system
#include <iostream>
#include <cstdlib>
#include <vexcl/vexcl.hpp>
using namespace std;

int main(int argc, char** argv) {
	if (argc != 3) { cout<<"usage: <iters per batch> <batches>"<<endl; return 1; }
	vex::Context ctx(vex::Filter::Type(CL_DEVICE_TYPE_GPU));
	cout << ctx << endl;

	uint N;
	string str;
	double *x, result = 0, size = atoi(argv[1]);
	uint max_iters = atoi(argv[2]);
        vex::vector<double> X(ctx, size), Y(ctx, size);
        vex::Reductor<double, vex::SUM> sum(ctx);
        double dx = double(1)/(double(size-1) * max_iters);
        cout<<"iters: "<<size<<endl;

	getline(cin, str, ',');
	x = new double[N = atoi(str.c_str())];
	for (uint k = 0; k < N; k++) {
		getline(cin, str, ',');
		cout << (x[k] = atoi(str.c_str())) << ',';
	}
	cout<<endl;

	for (uint iters = 0; iters <= max_iters; iters++) {
		Y = 1;
		X = (vex::constants::pi() * (vex::element_index() + double(iters)*size)) * dx;
		for (uint k = 0; k < N; k++)
			if (N - k > 4) {
				Y *= cos(X * x[k]) * cos(X * x[k + 1]) * cos(X * x[k + 2]);
				k += 3;
			} else Y *= cos(X * x[k]);
		Y *= vex::constants::pi() * dx;
		cout << iters << '/' << max_iters << ":\t" << (result += sum(Y)) <<endl;
	}
	return 0;
}
