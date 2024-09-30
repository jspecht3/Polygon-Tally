/* This is the c++ implementation of the KBasis and ZBasis
 * classes that I previously made in python. This
 * implementation is to check that I know how to code in c++.
 * There should be agreement with the ck coefficients for the
 * transformed zernike basis functions. */

#include <iostream>
#include <string>
#include <Eigen/Dense>
#include "zernike.h"

using Eigen::MatrixXd;
using namespace std;

int main() {
/*
	MatrixXd rho(1,1);
	MatrixXd phi(1,1);

	rho(0,0) = 0.5;
	phi(0,0) = 3.14159;

	MatrixXd coords[2] = {rho, phi};

	MatrixXd m(2,2);
	m(0,0) = 3;
	m(1,0) = 2.5;
	m(0,1) = -1;
	m(1,1) = m(1,0) + m(0,1);

	std::cout << coords[0] << std::endl;
*/
	z1 = Zernike();
	z1.n_nm();
	cout << z1.n_nm_ << endl;

	return 0;
}
