#ifndef _ALIGNMENT_HPP_
#define _ALIGNMENT_HPP_

#include <string>
#include <vector>
#include <map>
#include <tuple>

#include <cstdint>

#include <reference.hpp>
#include <proper_read.hpp>

class alignment
{
	private:
		reference _reference;
		
		// extracted primerID with a map to all PCR replicates
		std::map<std::string, std::vector<proper_read> > raw_primerID_map;
		
		// primerID map with collisions/PCR-mutations-in-PrimerID removed
		std::map<std::string, std::vector<proper_read>* > collision_free_primerID_map;
		
		// primerID map to consensus sequence
		std::map<std::string, consensus_read> consensus_primerID_map;
		
		std::string call_consensus_and_remove_collisions(const std::vector<proper_read>& reads, int minDisplay, const std::string& PrimerID = "");
		
	public:
		int number_singletons;
		int number_collisions;
		int number_indecisive;
		
		std::string input_fileName;
		std::string just_fileName;
		
		int _min_current_coverage;
		double _min_majority_fraction;
		
	// member functions:
		alignment(const std::string& fileName);
		
		void remove_primerID_collisions(int minC, double minMajorFraction, bool report = true, int minDisplay = 0);
		
		// RT stuff
		std::tuple<uint64_t, uint64_t, uint64_t> calculate_RT_mismatches() const;
		
		long double LogLik(long double s, long double r) const;
		
		// PCR stuff
		std::tuple<double, uint64_t, uint64_t> calculate_PCR_mismatches() const;
		
		// variance/overdispersion stuff
		void display_raw_and_primerID_counts() const;
		
		// chao estimator
		void calculate_effective_RNA_number() const;
		
		// diagnostics
		void show_primerIDs_with_min_coverage(int minC) const;
		void write_consensus_to_fasta() const;
		void write_to_fasta() const;
};


class alignments
{
	private:
		std::vector<alignment> collections_alignments;
		
		int _min_current_coverage;
		double _min_majority_fraction;
		
		bool is_collisions_removed;
		
		long double calculate_RT_LogLik_recombination(double s, double r) const;
		long double neg_LogLik_recombination(double s, double r) const;
		
	public:
		double _RT_sub_rate;
		double _RT_recomb_rate;
		
		double _PCR_sub_rate;
		
		alignments(const std::vector<std::string>& inputFiles);
		
		void remove_primerID_collisions(int minC = 10, double minMajorFraction = 0.9, bool printOut = false);
		void calculate_effective_population_size();
		
		// RT stuff
		double estimate_RT_substitution_rate(bool report = false);
		
		double estimate_RT_recombination_rate(bool report = false);
		double estimate_RT_recombination_rate(double s, bool report = false);
		
		// PCR stuff
		double estimate_PCR_substitution_rate(bool report = false) const;
		double estimate_PCR_recombinant_fraction(bool report = false) const;
		
		// variance/overdispersion stuff
		void display_raw_and_primerID_counts() const;
		
		// chao estimator
		void estimate_effective_RNA_number(bool report = false) const;
		
		// diagnostics
		void show_primerIDs_with_min_coverage(int minC = 100) const;
		void write_all_consensus_to_fasta() const;
		void write_all_to_fasta() const;
};

#endif /* _ALIGNMENT_HPP_ */