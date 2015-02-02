#ifndef _STATISTICS_HPP_
#define _STATISTICS_HPP_

#include <cmath>
#include <string>
#include <vector>
#include <valarray>
#include <algorithm>
#include <functional>
#include <array>
#include <iomanip>
#include <numeric>
#include <fstream>
#include <cstddef>

#include <gsl/gsl_cdf.h>


#ifdef REVERSE_ENDIAN
#define _ORDER L-1-
#else
#define _ORDER
#endif

inline std::string number_To_DNA(int number, int L)
{
	// map from numeric placeholder to ASCII char
	// 0 -> 'A'
	// 1 -> 'C'
	// 2 -> 'G'
	// 3 -> 'T'
	const static char arr[] = {'A', 'C', 'G', 'T'};
	
	std::string result(L, 'A');
	int j;
	
	for (int i = L-1; i >= 0; --i)
	{
		j = std::pow(4, i);
		result[_ORDER i] = arr[number / j];
		number %= j;
	}
	
	return result;
}

inline int DNA_to_number(const std::string& DNA)
{
	// map from ASCII char to numeric placeholder
	// 'A' -> 0
	// 'C' -> 1
	// 'G' -> 2
	// 'T' -> 3
	const static char arr[] =
		{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
		 -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
		 -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
		 -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0, -1,  1, 
		 -1, -1, -1,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  3};
	const int L = DNA.length();
	
	int result = 0;
	for (int i = 0; i < L; ++i)
	{
		result += std::pow(4, i) * arr[static_cast<unsigned char>(DNA[_ORDER i])];
	}
	
	return result;
}

template <typename T = int>
class DNAvector : public std::valarray<T>
{
	public:
		size_t _L;
		
		DNAvector(int L) :
			std::valarray<T>(static_cast<T>(0), static_cast<int>(std::pow(4, L))),
			_L(L) {}
		
		using std::valarray<T>::operator[];
		using std::valarray<T>::operator+=;
		using std::valarray<T>::operator/=;
		
		// indexing with DNA sequences
		inline T& operator[] (const std::string& str)
		{
			return this->operator[] (DNA_to_number(str));
		}
		
		inline const T& operator[] (const std::string& str) const
		{
			return this->operator[] (DNA_to_number(str));
		}
};

// fix for libc++ bug
#ifdef _LIBCPP_VERSION
template<class _Tp>
struct std::__is_val_expr<DNAvector<_Tp> > : std::true_type {};
#endif /* _LIBCPP_VERSION */

class seq_statistics
{
	private:
		std::vector<DNAvector<int>> _marginal_counts;
		std::vector<DNAvector<int>> _marginal_counts_repl;
		
		std::vector<DNAvector<int>> _pairwise_counts;
		std::vector<DNAvector<int>> _pairwise_counts_repl;
		
		DNAvector<int> _whole_pID_counts;      // occurrences of RT primers
		DNAvector<int> _whole_pID_counts_repl; // occurrences of RT primers, weighted by PCR replicates
		
		std::vector<int> _histogram_of_abundance;
	
	public:
		size_t _L;
		int _unique_counts;
		int _repl_counts;
		
		seq_statistics(int L) :
			_marginal_counts(L, DNAvector<int>(1)),
			_marginal_counts_repl(L, DNAvector<int>(1)),
			
			_pairwise_counts(L-1, DNAvector<int>(2)),
			_pairwise_counts_repl(L-1, DNAvector<int>(2)),
			
			_whole_pID_counts(L),
			_whole_pID_counts_repl(L),
			
			_histogram_of_abundance(10000, 0),
			
			_L(L),
			_unique_counts(0),
			_repl_counts(0) {}
		
		void reset()
		{
			_marginal_counts.assign(_L, DNAvector<int>(1));
			_marginal_counts_repl.assign(_L, DNAvector<int>(1));
			
			_pairwise_counts.assign(_L-1, DNAvector<int>(2));
			_pairwise_counts_repl.assign(_L-1, DNAvector<int>(2));
			
			_whole_pID_counts = DNAvector<int>(1);
			_whole_pID_counts_repl = DNAvector<int>(1);
			
			_histogram_of_abundance.assign(10000, 0);
		}
		
		void show_statistics() const
		{
			const static std::array<char, 4> DNA{{'A','C','G','T'}};
			
			for(size_t i = 0; i < 4; ++i)
			{
				std::cout << std::fixed << std::setprecision(2) << DNA[i] << ": " << _marginal_counts[0][i];
				
			  for(size_t j = 1; j < _L; ++j)
				{
					std::cout << '\t' << _marginal_counts[j][i];
				}
				
				std::cout << '\n';
			}
		}
		
		void addPrimer(const std::string& strPrimer, int replicates)
		{
			// 1.) add to main effects
			for (size_t i = 0; i < _L; ++i)
			{
				++_marginal_counts[i][strPrimer.substr(i, 1)];
				  _marginal_counts_repl[i][strPrimer.substr(i, 1)] += replicates;
			}
			
			// 2.) add to pairwise effects
			for (size_t i = 0; i < _L-1; ++i)
			{
				++_pairwise_counts[i][strPrimer.substr(i, 2)];
				  _pairwise_counts_repl[i][strPrimer.substr(i, 2)] += replicates;
			}
			
			// 3.) add to total RT count
			++_whole_pID_counts[strPrimer];
			  _whole_pID_counts_repl[strPrimer] += replicates;
			
			// 4.) add to abundance
			++_histogram_of_abundance[replicates];
			
			// 5.) keep track of sum
			++_unique_counts;
			_repl_counts += replicates;
		}
		
		void mergestatistics(const seq_statistics& statisticsB)
		{
			// 1.) sum up main effects
			for (size_t i = 0; i < _L; ++i)
			{
				_marginal_counts[i] += statisticsB._marginal_counts[i];
				_marginal_counts_repl[i] += statisticsB._marginal_counts_repl[i];
			}
			
			// 2.) sum up pairwise effects
			for (size_t i = 0; i < _L-1; ++i)
			{
				_pairwise_counts[i] += statisticsB._pairwise_counts[i];
				_pairwise_counts_repl[i] += statisticsB._pairwise_counts_repl[i];
			}
			
			// 3.) sum up total RT count
			_whole_pID_counts += statisticsB._whole_pID_counts;
			_whole_pID_counts_repl += statisticsB._whole_pID_counts_repl;
			
			_unique_counts += statisticsB._unique_counts;
			_repl_counts += statisticsB._repl_counts;
		}
		
		void write_to_csv(const std::string& file_stem) const
		{
			//const static std::array<char, 4> DNA{{'A','C','G','T'}};
			
			std::ofstream output     ((file_stem + "_unique.csv").c_str());
			std::ofstream output_repl((file_stem + "_replicates.csv").c_str());
			
			for(size_t i = 0; i < 4; ++i)
			{
				output      << std::fixed << std::setprecision(5) << _marginal_counts[0][i]      / static_cast<double>(_unique_counts);
				output_repl << std::fixed << std::setprecision(5) << _marginal_counts_repl[0][i] / static_cast<double>(_repl_counts);
				
			  for(size_t j = 1; j < _L; ++j)
				{
					output      << ',' << _marginal_counts[j][i]      / static_cast<double>(_unique_counts);
					output_repl << ',' << _marginal_counts_repl[j][i] / static_cast<double>(_repl_counts);
				}
				
				output      << '\n';
				output_repl << '\n';
			}
			
		  output.close();
			output_repl.close();
		}
		
		void calculate_comprehensive_statistics(const std::string& strLabel) const
		{
			const static std::array<std::string, 4> DNA{{"A","C","G","T"}};
			const static std::array<std::string, 3> trun_DNA{{"C","G","T"}};
			
			std::vector<DNAvector<double>> _marginal_probabilities(_L, DNAvector<double>(1));
			
			std::vector<DNAvector<double>> _pairwise_probabilities(_L-1, DNAvector<double>(2));
			
			// 1.) calculate (MLE) probabilities
			for(size_t i = 0; i < _L; ++i)
			{
				// marginals
				for(size_t j = 0; j < _marginal_counts[i].size(); ++j)
				{
					_marginal_probabilities[i][j] = _marginal_counts[i][j] / static_cast<double>(_unique_counts);
				}
				
				// pairwise
				if (i != _L-1)
				{
					for(size_t j = 0; j < _pairwise_counts[i].size(); ++j)
					{
						_pairwise_probabilities[i][j] = _pairwise_counts[i][j] / static_cast<double>(_unique_counts);
					}
				}
			}
			
			// 2.) print results
			if (!strLabel.empty())
			{
				std::cout << strLabel << ":";
			}
			
			// numbering
			for(size_t i = 0; i < _L; ++i)
				std::cout << "\t\t" << i+1 << "\t";
			std::cout << '\n';
			
			for(size_t i = 0; i < _L; ++i)
				std::cout << "\tp=\tbeta=\t";
			std::cout << '\n';
			
			double beta;
			/* marginals */
			for(const auto& j : DNA)
			{
				std::cout << j << ":";
				for(size_t i = 0; i < _L; ++i)
				{
					beta = std::log(_marginal_probabilities[i][j] / _marginal_probabilities[i]["A"]);
					
					std::cout << "\t" << std::fixed << std::setprecision(4) << _marginal_probabilities[i][j] << "\t" << (beta < 0 ? "" : " ") << beta << "\t";
				}
				std::cout << '\n';
			}

			std::vector<DNAvector<double>> beta_1(_L-1, DNAvector<double>(1));
			std::vector<DNAvector<double>> beta_2(_L-1, DNAvector<double>(1));
			
			// main effects
			// 1st locus
			for(const auto& j : trun_DNA)
			{
				std::cout << "b_1(" << j << "):\t";
				
				for(size_t i = 0; i < _L-1; ++i)
				{
					beta = std::log(_pairwise_probabilities[i][j+"A"] / _pairwise_probabilities[i]["AA"]);
					beta_1[i][j] = beta;
					std::cout << "\t\t" << (beta < 0 ? "" : " ") << beta << "\t";
				}
				
				std::cout << '\n';
			}
			
			// 2nd locus
			for(const auto& j : trun_DNA)
			{
				std::cout << "b_2(" << j << "):\t";
				
				for(size_t i = 0; i < _L-1; ++i)
				{
					beta = std::log(_pairwise_probabilities[i]["A"+j] / _pairwise_probabilities[i]["AA"]);
					beta_2[i][j] = beta;
					std::cout << "\t\t" << (beta < 0 ? "" : " ") << beta << "\t";
				}
				
				std::cout << '\n';
			}
			
			// interaction effects
			for(const auto& j2 : trun_DNA)
			{
				for(const auto& j1 : trun_DNA)
				{
					std::cout << "b_12(" << j1 << j2 << "):";
					
					for(size_t i = 0; i < _L-1; ++i)
					{
						beta = std::log(_pairwise_probabilities[i][j1+j2] / _pairwise_probabilities[i]["AA"]) - beta_1[i][j1] - beta_2[i][j2];
						std::cout << "\t\t" << (beta < 0 ? "" : " ") << beta << "\t";
					}
					
					std::cout << '\n';
				}
			}
			
			// perform independence test
			std::vector<double> LRs(9, 0);
			std::cout << "G:\t";
			for(size_t i = 0; i < _L-1; ++i)
			{
				for(const auto& j2 : DNA)
				{
					for(const auto& j1 : DNA)
					{
						LRs[i] += _pairwise_counts[i][j1+j2]*std::log(_pairwise_probabilities[i][j1+j2] / (_marginal_probabilities[i][j1]*_marginal_probabilities[i+1][j2]));
					}
				}
				
				LRs[i] *= 2;
				std::cout << "\t" << LRs[i] << " (p=" << 1-gsl_cdf_chisq_P(LRs[i], 9) << ")";
			}
			std::cout << "\n\n";
		}
};


// UNIT TEST:
/*
#include <iostream>
int main()
{
	const int L = 2;
	DNAvector<int> LOLvector(L);
	
	std::string DNA1;
	std::string DNA2;
	int number1;
	int number2;
	int TEMP;
	
	for (int i = 0; i < 4*4; ++i)
	{
		TEMP = 2*i;
		LOLvector[i] = TEMP;
		
		DNA1 = number_To_DNA(i, L);
		number1 = DNA_to_number(DNA1);
		DNA2 = number_To_DNA(number1, L);
		number2 = DNA_to_number(DNA2);
		
		std::cout << "i: " << i << "\tDNA1: " << DNA1 << "\tnumber1: " << number1 << "\tDNA2: " << DNA2 << "\tnumber2: " << number2 << "\tTEMP: " << TEMP << "  \tVECTOR: " << LOLvector[DNA2] << '\n';
		
		if ((i != number1) || (number1 != number2))
			std::cout << "ERROR!\n";
		if (DNA1 != DNA2)
			std::cout << "ERROR!\n";
	}
}
*/

#endif /* _STATISTICS_HPP_ */