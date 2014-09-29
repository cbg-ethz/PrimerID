#include <cmath>

#include <proper_read.hpp>

// proper_read
proper_read::proper_read(const std::string& input, const reference& _ref) : _reference(_ref), fullRead(input)
{
	no_of_valid_heterozygous_bases = _ref.no_heterozygous_loci;
	
	// construct heterozygous string
	for (int i = 0; i < _ref.no_heterozygous_loci; ++i)
	{
		const char& allele = fullRead[_ref.heterozygous_loci[i]];

		if (allele != 'N')
			indices_valid_heterozygous.push_back(i);
		else
			--no_of_valid_heterozygous_bases;
		
		heterozygous_loci_string.push_back(allele);
	}
	
	// construct homozygous string
	for (int i = 0; i < _ref.no_homozygous_loci; ++i)
	{
		homozygous_loci_string.push_back(fullRead[_ref.homozygous_loci[i]]);
	}

	// find closest reference strain	
	int best_index = _ref.K;
	int best_ham = 100000;
	int cur_ham;
	
	for (int i = 0; i < _ref.K; ++i)
	{
		cur_ham = hamming_distance(heterozygous_loci_string, _ref.all_reference_strains[i].heterozygous_loci_string);
		
		if (cur_ham < best_ham)
		{
			best_ham = cur_ham;
			best_index = i;
		}
	}
	
	best_reference = best_index;
	hamming_distance_to_best_reference = best_ham;
}

// consensus_read
consensus_read::consensus_read(const std::string& input, const reference& _ref, int _multiplicity) : proper_read(input, _ref), multiplicity(_multiplicity) {}

inline long double emission_probability(char X_i, char Z_i, long double s)
{
	return (X_i == Z_i ? 1-s : s/3);
}

inline long double transition_probability(char Z_i, char Z_i_min_1, int n_steps, long double r, const reference& _ref)
{
	long double p_norecombination = pow(1-r, n_steps);
	long double p_ij = _ref.all_reference_strains[Z_i].frequency + p_norecombination * ((Z_i == Z_i_min_1) - _ref.all_reference_strains[Z_i].frequency);
	
	return p_ij;
}

long double consensus_read::log_prob(long double s, long double r) const
{
	// performs the full HMM calculation for one sequence
	// FORWARD algorithm: P(X) = \sum_{Z} P(X, Z)
	
	if (no_of_valid_heterozygous_bases > 1)
	{
		// probability is a function of r, i.e. the likelihood
		// of this sequence can be used to make inference on r
		std::vector<long double> f_old(_reference.K), f_new(_reference.K);
		
		// 1.) initialize
		for (int j = 0; j < _reference.K; ++j)
			f_new[j] = _reference.all_reference_strains[j].frequency * emission_probability(heterozygous_loci_string[indices_valid_heterozygous[0]], _reference.all_reference_strains[j].heterozygous_loci_string[indices_valid_heterozygous[0]], s);
		
		// 2.) recursion
		long double sum;
		int jump;
		
		for (int i = 1; i < no_of_valid_heterozygous_bases; ++i)
		{
			f_new.swap(f_old);
			jump = indices_valid_heterozygous[i] - indices_valid_heterozygous[i-1];
			
			for (int j = 0; j < _reference.K; ++j)
			{
				sum = 0;
				
				for (int k = 0; k < _reference.K; ++k)
				{
					sum += transition_probability(j, k, jump, r, _reference) * f_old[k];
				}
				f_new[j] = emission_probability(heterozygous_loci_string[indices_valid_heterozygous[i]], _reference.all_reference_strains[j].heterozygous_loci_string[indices_valid_heterozygous[i]], s) * sum;
			}
		}
		
		// 3.) termination
		sum = 0;
		for (int j = 0; j < _reference.K; ++j)
			sum += f_new[j];
		
		return log(sum);
	}
	else
	{
		// either only or just one proper base
		// does not contribute to likelihood
		// except for vertical shift
		return 0;
	}
}