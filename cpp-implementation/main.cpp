/* Radius transform from polygonal coordinates into their 
 * respective unit disk coordinates */

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <array>
#include <math.h>

// converts cartesian to polar coordinates
std::array<double,2> c2p(std::array<double,2> cart) {
	double x = cart[0];
	double y = cart[1];

	double r = std::sqrt(x * x + y * y);
	double theta = std::atan2(y, x);

	std::array<double,2> polar;
	polar[0] = r;
	polar[1] = theta;

	return polar;
}

// intermediate function u_alpha for variable radius calc
double u_alpha(double theta, int num_sides) {
	double alpha = M_PI / num_sides;
	double drop = (theta + alpha) / (2 * alpha);
	int int_drop = drop;

	double u_alpha;
	u_alpha = theta - int_drop * 2 * alpha;

	return u_alpha;
}

// takes the cart coords and returns the variable radius
double r_alpha(
		std::array<double,2> cart,
		double polygon_radius,
		int num_sides) {
	
	// get the polar coords
	std::array<double,2> polar = c2p(cart);
	double r = polar[0];
	double theta = polar[1];

	// calculating
	double alpha = M_PI / num_sides;

	double r_alpha;
	r_alpha = polygon_radius * cos(alpha) / 
		cos(u_alpha(theta, num_sides));

	return r_alpha;
}


int main() {
	// getting the values from the customer
	std::cout << "----- Input -----" << std::endl;

	int num_sides;
	std::cout << "Number of sides: ";
	std::cin >> num_sides;

	double polygon_radius;
	std::cout << "Polygon Radius: ";
	std::cin >> polygon_radius;

	// getting coordinates
	std::array<double,2> cart;
	
	std::cout << "x coord: ";
	std::cin >> cart[0];

	std::cout << "y coord: ";
	std::cin >> cart[1];

	// converting
	std::array<double,2> output;
    output = c2p(cart);

	double r = { output[0] };
	double theta = { output[1] };

	double variable_radius;
	variable_radius = r_alpha(cart, polygon_radius, num_sides);

	// output
	std::cout << "----- Output -----" << std::endl;
	std::cout << "r: " << r << 
		"\ntheta: " << theta << 
		"\nvariable radius: " << variable_radius << 
		"\nrho: " << r / variable_radius << std::endl;

	return 0;
}