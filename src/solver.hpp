/**
 * solver.hpp
 *
 * Author: Viktor Fukala
 * Created on: 2020/11/7
 */
#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <cstdint>

namespace cyrill6
{

class Solver
{
public:
	double *y;
	Solver(const uint32_t n, const double cnv, const double l, const double lam, const double mu, const double g, const double rho);
	void solve();
	~Solver();

private:
	double lam;
	double mu, g, rho;
	double cnv;
	double last_change;
	double dx;
	uint32_t n;

	void iter();
};

}

#endif // SOLVER_HPP
