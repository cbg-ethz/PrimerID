#ifndef _PROPER_READ_HPP_
#define _PROPER_READ_HPP_

#include <vector>
#include <string>

#include <src/reference.hpp>

struct proper_read
{
  // reference stuff
  const reference&   _ref;
	const std::string& _fullRead;

  int best_reference;
  int hamming_distance_to_best_reference;

  int no_N;
  int no_of_valid_heterozygous_bases;

  /* actual sequence data */
  // heterozygous
  std::string heterozygous_loci_string;
  // std::vector<int> indices_valid_heterozygous;

  // homozygous
  //std::string homozygous_loci_string;

  /* member functions*/
  proper_read(const std::string& strDNA, const reference& _reference);

  //hamming_return_type calculate_homozygous_mismatches() const;
  hamming_return_type hetero_hamming_distance(const std::string& other_string) const;
  hamming_return_type hetero_hamming_distance(const proper_read& other_read) const;
};

struct consensus_read : public proper_read
{
	const std::string  _fullConsensus;
	int multiplicity;
	
  using proper_read::proper_read;
  consensus_read(std::string&& input, const reference& _reference, int _multiplicity);
	
	hamming_return_type calculate_homozygous_mismatches() const;

  // s:  substitution rate
  // r: recombination rate
  double log_prob(double s, double r) const;
};

#endif /* _PROPER_READ_HPP_ */