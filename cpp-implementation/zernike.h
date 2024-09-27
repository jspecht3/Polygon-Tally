// Class for the Zernike Polynomial
class Zernike {
	public:
		// variables
		int n_;
		int m_;

		double rho_;
		double phi_;

		double n_nm_;
		double r_nm_;
		double z_nm_;

		// constructors
		Zernike();
		Zernike(double rho, double phi, int n, int m);

		// functions
		void n_nm();
		void r_nm();
		void z_nm();

		void print();
};
