#include <cmath>
#include <tuple>

#include <proper_read.hpp>
#include <iostream>

// proper_read
proper_read::proper_read(const std::string& input, const reference& _reference) : _ref(_reference), no_of_valid_heterozygous_bases(0), fullRead(input)
{
	// 1.) construct heterozygous string
	
	for (int i : _ref.heterozygous_loci)
	{
		no_of_valid_heterozygous_bases += (fullRead[i] != 'N');
		heterozygous_loci_string.push_back(fullRead[i]);
	}
	
	// 2.) construct homozygous string
	for (int i : _ref.homozygous_loci)
		homozygous_loci_string.push_back(fullRead[i]);

	// 3.) find closest reference strain	
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

std::tuple<uint64_t, uint64_t, uint64_t> proper_read::calculate_homozygous_mismatches() const
{
	int mismatches = 0;
	int valid_trials = 0;
	int Ns = 0;
	
	for (int i = 0; i < homozygous_loci_string.length(); ++i)
	{
		if (homozygous_loci_string[i] != 'N')
		{
			++valid_trials;
			mismatches += (homozygous_loci_string[i] != _ref.homozygous_string[i]);
		}
		else
			++Ns;
	}
	
	return std::tuple<uint64_t, uint64_t, uint64_t>(mismatches, valid_trials, Ns);
}

std::tuple<uint64_t, uint64_t, uint64_t> proper_read::hetero_hamming_distance(const proper_read& other_read) const
{
	int mismatches = 0;
	int valid_trials = 0;
	int Ns = 0;
	
	for (int i = 0; i < _ref.no_heterozygous_loci; ++i)
	{
		if ((this->heterozygous_loci_string[i] != 'N') && (other_read.heterozygous_loci_string[i] != 'N'))
		{
			++valid_trials;
			mismatches += (this->heterozygous_loci_string[i] != other_read.heterozygous_loci_string[i]);
		}
		else
			++Ns;
	}
		
	return std::tuple<uint64_t, uint64_t, uint64_t>(mismatches, valid_trials, Ns);
}

// consensus_read
consensus_read::consensus_read(const std::string& input, const reference& _reference, int _multiplicity) : proper_read(input, _reference), multiplicity(_multiplicity) {}

inline long double emission_probability(char X_i, char Z_i, long double s)
{
	if (X_i == 'N')
		return 1;
	else
		return (X_i == Z_i ? 1-s : s/3);
}

inline long double transition_probability(char Z_i, char Z_i_min_1, int n_steps, long double r, const reference& _reference)
{
	long double p_norecombination = pow(1-r, n_steps);
	long double p_ij = _reference.all_reference_strains[Z_i].frequency + p_norecombination * ((Z_i == Z_i_min_1) - _reference.all_reference_strains[Z_i].frequency);
	
	return p_ij;
}

long double consensus_read::log_prob(long double s, long double r) const
{
	// performs the full HMM calculation for one sequence
	// FORWARD algorithm: P(X) = \sum_{Z} P(X, Z)
	
	// probability is a function of r, i.e. the likelihood
	// of this sequence can be used to make inference on r
	std::vector<long double> f_old(_ref.K), f_new(_ref.K);
	
	// 1.) initialize
	for (int j = 0; j < _ref.K; ++j)
		f_new[j] = _ref.all_reference_strains[j].frequency * emission_probability(heterozygous_loci_string[0], _ref.all_reference_strains[j].heterozygous_loci_string[0], s);
	
	// 2.) recursion
	long double sum;
	int jump;
	
	for (int i = 1; i < _ref.no_heterozygous_loci; ++i)
	{
		f_new.swap(f_old);
		jump = _ref.heterozygous_loci[i] - _ref.heterozygous_loci[i-1];
		
		for (int j = 0; j < _ref.K; ++j)
		{
			sum = 0;
			
			for (int k = 0; k < _ref.K; ++k)
			{
				sum += transition_probability(j, k, jump, r, _ref) * f_old[k];
			}
			f_new[j] = emission_probability(heterozygous_loci_string[i], _ref.all_reference_strains[j].heterozygous_loci_string[i], s) * sum;
		}
	}
	
	// 3.) termination
	sum = 0;
	for (int j = 0; j < _ref.K; ++j)
		sum += f_new[j];
	
	return log(sum);
}