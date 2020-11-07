/**
 * solver.cpp
 *
 * Author: Viktor Fukala
 * Created on: 2020/11/7
 */
#include "solver.hpp"

#include <algorithm>

#include <eigen3/Eigen/IterativeLinearSolvers>
#include <eigen3/Eigen/SparseCore>

namespace cyrill6
{

Solver::Solver(const uint32_t n, const double cnv, const double l, const double lam, const double mu, const double g, const double rho)
	: n(n), cnv(cnv), dx(l / (n - 1)), lam(lam), mu(mu), g(g), rho(rho), last_change(cnv + 1)
{
	y = new double[n];
	std::fill_n(y, n, 0);
}

Solver::~Solver()
{
	delete [] y;
}

void push_entry(std::vector<Eigen::Triplet<double>> &entries, const uint32_t eq, const uint32_t var, const double val)
{
	entries.push_back(Eigen::Triplet<double>(eq - 1, var - 1, val));
}

void Solver::iter()
{
	Eigen::VectorXd x(n - 2), b(n - 2);
	for (uint32_t i = 1; i < n - 1; i++)
	{
		b[i - 1] = mu * (y[i + 1] + y[i - 1] - 2 * y[i]) / (dx * dx)
			- (rho * g * y[i] + lam) * pow(1 + pow((y[i + 1] - y[i - 1]) / (2 * dx), 2), 1.5);
	}

	std::vector<Eigen::Triplet<double>> entries;
	for (uint32_t i = 1; i < n - 1; i++)
	{
		const double dyi = -2 * mu / (dx * dx)
			-rho * g * pow((1 + pow((y[i + 1] - y[i - 1]) / (2 * dx), 2)), 1.5);
		push_entry(entries, i, i, dyi);
		if (i != n - 2)
		{
			const double dyip = mu / (dx * dx) - (rho * g * y[i] + lam) * 1.5
				* sqrt(1 + pow((y[i + 1] - y[i - 1]) / (2 * dx), 2)) * (y[i + 1] - y[i - 1]) / (2 * dx) / dx;
			push_entry(entries, i, i + 1, dyip);
		}
		if (i != 1)
		{
			const double dyim = mu / (dx * dx) + (rho * g * y[i] + lam) * 1.5
				* sqrt(1 + pow((y[i + 1] - y[i - 1]) / (2 * dx), 2)) * (y[i + 1] - y[i - 1]) / (2 * dx) / dx;
			push_entry(entries, i, i - 1, dyim);
		}
	}

	Eigen::SparseMatrix<double> A(n - 2, n - 2);
	A.setFromTriplets(entries.begin(), entries.end());

	Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
	solver.compute(A);
	x = solver.solve(b);

	double changel1 = 0;
	for (uint32_t i = 0; i < n - 2; i++)
	{
		changel1 += abs(x[i]);
		y[i + 1] -= x[i];
	}
	last_change = changel1;

	// TODO implement solver using algebrasation and Newton solving using Eigen
}

void Solver::solve()
{
	while (last_change > cnv)
	{
		iter();
	}
}

}
