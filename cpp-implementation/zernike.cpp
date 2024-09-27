/* Implementation of the Zernike polynomials */

#include <iostream>
#include <cstdlib>
#include <cmath>
#include "zernike.h"

int fact(int x) {
	if (x <= 1) {
		return 1;
	} else {
		return (x * fact(x - 1));
	}
}

// default constructor
Zernike::Zernike() {
	n_ = 0;
	m_ = 0;

	rho_ = 0;
	phi_ = 0;
}

Zernike::Zernike(double rho, double phi, int n, int m) {
	n_ = n;
	m_ = m;

	rho_ = 0;
	phi_ = 0;
}

// functions
void Zernike::n_nm(){	
	n_nm_ = pow(((2 * (n_ + 1)) / (1 + (m_ == 0))), 0.5);
}

void Zernike::r_nm(){
	float value = 0;

	for(int k = 0; k <= (n_ - abs(m_)) / 2; k++) {
		float top = pow(-1, k) * fact(n_ - k);
		float bot = fact(k) *
			fact(((n_ + abs(m_)) / 2) - k) *
			fact(((n_ - abs(m_)) / 2) - k);
		value += top / bot * pow(rho_, n_ - 2 * k);
	}
	r_nm_ = value;
}

void Zernike::z_nm(){	
	if (m_ >= 0) {
		z_nm_ = n_nm_ * r_nm_ * cos(m_ * phi_);
	} else {
		z_nm_ =  -1 * n_nm_ * r_nm_ * sin(m_ * phi_);
	}
}
