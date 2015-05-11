#include <src/pidalyse/statistics.hpp>

#include <iostream>
#include <iomanip>

int main()
{
  int min = 2;
  int N = 100;

  int n = 8;
	
	prob_cycle PROB;

  for (int i = 1; i <= n; ++i)
  {
    std::cout << '\t' << i;
  }

  std::cout << "\tsum\n";

  for (int i = min + 1; i <= N / 2; ++i)
  {
    double sum = 0;
    double p;

    std::cout << i;

    for (int j = 1; j <= n; ++j)
    {
      p = PROB(i, N, min, j);
      sum += p;

      std::cout << '\t' << std::fixed << std::setprecision(4) << p;
    }

    std::cout << '\t' << sum << '\n';
  }

  std::cout << "Number of entries in cache: " << PROB.getSize() << '\n';

  double sum = 0;
  for (int i = min + 1; i <= N / 2; ++i)
  {
    sum += PROB.p_X_given_cycle_and_constraint(i, N, min, 4);
  }

  std::cout << "Sum for i = " << 1 << ": " << sum << '\n';

  return 0;
}