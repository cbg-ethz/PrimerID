#ifndef _PROPER_READ_HPP_
#define _PROPER_READ_HPP_

#include <vector>
#include <string>

#include <reference.hpp>

struct proper_read
{
	// reference stuff
	const reference& _reference;
	
	int best_reference;
	int hamming_distance_to_best_reference;
	int no_of_valid_heterozygous_bases;
	
	/* actual sequence data */
	std::string fullRead;
	
	// heterozygous
	std::string heterozygous_loci_string;
	std::vector<int> indices_valid_heterozygous;
	
	// homozygous
	std::string homozygous_loci_string;
	
	/* member functions*/
	proper_read(const std::string& input, const reference& _ref);
};

struct consensus_read : public proper_read
{
	using proper_read::proper_read;
	consensus_read(const std::string& input, const reference& _ref, int _multiplicity);
	
	// s:  substitution rate
	// r: recombination rate
	long double log_prob(long double s, long double r) const;
	int multiplicity;
};

#endif /* _PROPER_READ_HPP_ */