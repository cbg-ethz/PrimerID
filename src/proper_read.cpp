#include <cmath>
#include <tuple>
#include <functional>
#include <exception>

#include <src/proper_read.hpp>

hamming_return_type
hamming_distance(const std::string& sA, const std::string& sB)
{
  int mismatches = 0;
  int valid_trials = 0;
  int Ns = 0;

  if (sA.length() != sB.length())
    throw std::range_error("String lengths do not match!\n");

  for (size_t i = 0; i < sA.length(); ++i)
  {
    if ((sA[i] != 'N') && (sB[i] != 'N'))
    {
      ++valid_trials;
      mismatches += (sA[i] != sB[i]);
    }
    else
      ++Ns;
  }

  return hamming_return_type(mismatches, valid_trials, Ns);
}

// proper_read
proper_read::proper_read(const std::string& strDNA, const reference& _reference)
    : _ref(_reference),
			_fullRead(strDNA),
      no_N(0)
{
  // 1.) construct heterozygous string
  for (int i : _ref.heterozygous_loci)
  {
    no_N += (_fullRead[i] == 'N');
    heterozygous_loci_string.push_back(_fullRead[i]);
  }

  no_of_valid_heterozygous_bases = _ref.no_heterozygous_loci - no_N;

	/*
  // 2.) construct homozygous string
  for (int i : _ref.homozygous_loci)
    homozygous_loci_string.push_back(_fullRead[i]);
	*/

  // 3.) find closest reference strain
  int best_index = _ref.K;
  int best_ham = 100000;

  int cur_ham;
  int valid_trials;
  int Ns;

  for (int i = 0; i < _ref.K; ++i)
  {
    std::tie(cur_ham, valid_trials, Ns) = hamming_distance(
        heterozygous_loci_string,
        _ref.all_reference_strains[i].heterozygous_loci_string);

    if (cur_ham < best_ham)
    {
      best_ham = cur_ham;
      best_index = i;
    }
  }

  // 4.) assign RECOMBINANT if hamming distance to known is too large
  hamming_distance_to_best_reference = best_ham;
  if (hamming_distance_to_best_reference >= 2)
  {
    // suspect, possibly recombinant
    best_reference = _reference.K;
  }
  else
  {
    // probably no recombinant
    best_reference = best_index;
  }
}

/*
hamming_return_type
proper_read::calculate_homozygous_mismatches() const
{
  int mismatches = 0;
  int valid_trials = 0;
  int Ns = 0;

  for (int i : _ref.homozygous_loci)
  {
    if (_fullRead[i] != 'N')
    {
      ++valid_trials;
      mismatches += (_fullRead[i] != _ref.all_reference_strains[0].DNA[i]);
    }
    else
      ++Ns;
  }

  return hamming_return_type(mismatches, valid_trials, Ns);
}
*/

hamming_return_type
proper_read::hetero_hamming_distance(const std::string& other_string) const
{
  return hamming_distance(this->heterozygous_loci_string, other_string);
}

hamming_return_type
proper_read::hetero_hamming_distance(const proper_read& other_read) const
{
  return hetero_hamming_distance(other_read.heterozygous_loci_string);
}

// consensus_read
consensus_read::consensus_read(std::string&& strDNA, const reference& _reference, int _multiplicity)
	: proper_read(strDNA, _reference),
		_fullConsensus(std::move(strDNA)),
		multiplicity(_multiplicity) {}

hamming_return_type consensus_read::calculate_homozygous_mismatches() const
{
  int mismatches = 0;
  int valid_trials = 0;
  int Ns = 0;

  for (int i : _ref.homozygous_loci)
  {
    if (_fullConsensus[i] != 'N')
    {
      ++valid_trials;
      mismatches += (_fullConsensus[i] != _ref.all_reference_strains[0].DNA[i]);
    }
    else
      ++Ns;
  }

  return hamming_return_type(mismatches, valid_trials, Ns);
}

inline long double emission_probability(char X_i, char Z_i, long double s)
{
  if (X_i == 'N')
    return 1;
  else
    return (X_i == Z_i ? 1 - s : s / 3);
}

inline long double transition_probability(char Z_i, char Z_i_min_1, int n_steps,
                                          long double r,
                                          const reference& _reference)
{
  long double p_norecombination = pow(1 - r, n_steps);
  long double p_ij = _reference.all_reference_strains[Z_i].frequency + p_norecombination * ((Z_i == Z_i_min_1) - _reference.all_reference_strains[Z_i].frequency);

  return p_ij;
}

double consensus_read::log_prob(double s, double r) const
{
  // performs the full HMM calculation for one sequence
  // FORWARD algorithm: P(X) = \sum_{Z} P(X, Z)
	// complexity: O(L*K^2), i.e., approx 25*L

  // probability is a function of r, i.e. the likelihood
  // of this sequence can be used to make inference on r
  std::vector<long double> f_old(_ref.K), f_new(_ref.K);

  // 1.) initialize
  for (int j = 0; j < _ref.K; ++j)
    f_new[j] = _ref.all_reference_strains[j].frequency * emission_probability(
                                                             heterozygous_loci_string[0],
                                                             _ref.all_reference_strains[j].heterozygous_loci_string[0], s);

  // 2.) recursion
  long double sum;
  int n_jump;

  for (int i = 1; i < _ref.no_heterozygous_loci; ++i)
  {
    f_new.swap(f_old);
    n_jump = _ref.heterozygous_loci[i] - _ref.heterozygous_loci[i - 1];

    for (int j = 0; j < _ref.K; ++j)
    {
      sum = 0;

      for (int k = 0; k < _ref.K; ++k)
      {
        sum += transition_probability(j, k, n_jump, r, _ref) * f_old[k];
      }

      f_new[j] = emission_probability(
                     heterozygous_loci_string[i],
                     _ref.all_reference_strains[j].heterozygous_loci_string[i], s) * sum;
    }
  }

  // 3.) termination
  sum = 0;
  for (int j = 0; j < _ref.K; ++j)
    sum += f_new[j];

  return log(sum);
}