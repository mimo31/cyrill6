#include <cstdlib>
#include <iostream>

#include "solver.hpp"

using std::cout;
using std::endl;

namespace cyrill6
{

void run()
{
	constexpr uint32_t n = 256;
	constexpr double lam = -1;
	constexpr double g = 400;

	Solver s(n, .1, 1, lam, 1, g, 1);
	s.solve();

	for (uint32_t i = 0; i < n; i++)
	{
		cout << i / double(n - 1) << " " << s.y[i] << endl;
	}
}

}

int main()
{
	cyrill6::run();
	return EXIT_SUCCESS;
}
