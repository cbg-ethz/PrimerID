#ifdef _ALIGNMENT_CPP_

//#include <iostream>
//#include <iomanip>

#include <cmath>
#include <tuple>
#include <map>

//#include <gsl/gsl_randist.h>
//#include <gsl/gsl_cdf.h>

#include <boost/math/distributions/binomial.hpp>

class prob_cycle
{
	public:
		double operator() (int no_mt, int no_total, int no_min, int cycle = 1) const;
		int getSize() const;

		double _prob_cycle(int cycle) const;
		double _prob_X_given_cycle_and_constraint(int no_mt, int no_total, int no_min, int cycle) const;
		double _prob_cycle_given_X_and_constraint(int no_mt, int no_total, int no_min, int cycle) const;
		
	private:
		mutable std::map<std::tuple<int, int, int, int>, double> _cached_results;
};

int prob_cycle::getSize() const
{
	return _cached_results.size();
}

double prob_cycle::_prob_X_given_cycle_and_constraint(int no_mt, int no_total, int no_min, int cycle) const
{
	double p = exp2(-static_cast<double>(cycle));
	
	//return gsl_ran_binomial_pdf(no_mt, p, no_total) / (gsl_cdf_binomial_P(no_total/2, p, no_total) - gsl_cdf_binomial_P(no_min, p, no_total));
	return pdf(boost::math::binomial(no_total, p), no_mt) / (cdf(boost::math::binomial(no_total, p), no_total/2) - cdf(boost::math::binomial(no_total, p), no_min));
}

double prob_cycle::_prob_cycle(int cycle) const
{
	return exp2(static_cast<double>(cycle)-1);
}

double prob_cycle::_prob_cycle_given_X_and_constraint(int no_mt, int no_total, int no_min, int cycle) const
{
	double sum = 0;
	
	for (int j = 1; j < 12; ++j)
	{
		sum+= _prob_X_given_cycle_and_constraint(no_mt, no_total, no_min, j) * _prob_cycle(j);
	}
	
	return _prob_X_given_cycle_and_constraint(no_mt, no_total, no_min, cycle) * _prob_cycle(cycle) / sum;
}

double prob_cycle::operator() (int no_mt, int no_total, int no_min, int cycle) const
{
	std::map<std::tuple<int, int, int, int>, double>::const_iterator it = _cached_results.find(std::forward_as_tuple(no_mt, no_total, no_min, cycle));
	if (it != _cached_results.end())
	{
		return it->second;
	}
	else
	{
		double p = _prob_cycle_given_X_and_constraint(no_mt, no_total, no_min, cycle);
		_cached_results.emplace(std::make_tuple(no_mt, no_total, no_min, cycle), p);
		
		return p;
	}
}

static prob_cycle PCR_prob;

#ifdef UNIT_TEST
int main()
{
	int min = 2;
	int N = 15;
	
	int n = 8;	
	
	for (int i = 1; i <= n; ++i)
	{
		std::cout << '\t' << i;
	}
	std::cout << "\tsum\n";
	
	for(int i = min+1; i <= N/2; ++i)
	{
		double sum = 0;
		double p;
		
		std::cout << i;
		
		for (int j = 1; j <= n; ++j)
		{
			p = PROB(i, N, min, j);
			sum += p;
			
			std::cout << '\t'<< std::fixed << std::setprecision(3) << p;
		}
		
		std::cout << '\t' << sum << '\n';
	}
	
	std::cout << "Number of entries in cache: " << PROB.getSize() << '\n';
	
	
	double sum = 0;
	for(int i = min+1; i <= N/2; ++i)
	{
		sum += PROB._prob_X_given_cycle_and_constraint(i, N, min, 4);
	}
	std::cout << "Sum for i = " << 1 << ": " << sum << '\n';
	
	
	return 0;
}
#endif /*  UNIT_TEST */
#endif /* _ALIGNMENT_CPP_ */