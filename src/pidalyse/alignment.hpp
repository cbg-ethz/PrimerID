#ifndef _ALIGNMENT_HPP_
#define _ALIGNMENT_HPP_

#include <string>
#include <vector>
#include <map>
#include <tuple>

#include <cstddef>
#include <cstdint>

#include "reference.hpp"
#include "proper_read.hpp"
#include "statistics.hpp"

class alignment
{
private:
  reference _reference;
	
	// store raw reads from SAM data
	struct raw_paired_read {
		const std::string DNA;
		const std::string PrimerID;
		int no_inserts;
		int no_dels;
		int len_overhang;
		int len_pID;
		
		raw_paired_read(
			std::string&& strDNA,
			std::string&& strPrimerID,
			int i_no_inserts,
			int i_no_dels,
			int i_len_overhang,
			int i_len_pID) :
		DNA(std::move(strDNA)),
		PrimerID(std::move(strPrimerID)),
		no_inserts(i_no_inserts),
		no_dels(i_no_dels),
		len_overhang(i_len_overhang),
		len_pID(i_len_pID) {}
	};
	
	// collection of raw PAIRED reads
	std::vector<raw_paired_read> _raw_reads;
	
	// map of primerID -> PCR ensembl
  std::map<std::string, std::vector<proper_read>  > raw_primerID_map;

  // map of primerID -> collision-free PCR ensembl
  std::map<std::string, std::vector<proper_read*> > collision_free_primerID_map;

  // map of primerID -> consensus sequence of collision-free PCR ensembl
  std::map<std::string, consensus_read> consensus_primerID_map;

	std::string call_consensus_and_remove_collisions(std::vector<proper_read>& reads, int minDisplay, const std::string& PrimerID, std::vector<proper_read*>& filtered_reads);
	
public:
	// primerID statistics
	seq_statistics _pid_stats;
	
  int number_singletons;
  int number_collisions;
  int number_indecisive;

  int mismatches;
  int valid_trials;
  int Ns;

  std::string input_fileName;
  std::string just_fileName;
	std::string fileStem;

  size_t _min_current_coverage;
  double _min_majority_fraction;

  // member functions:
  alignment(const std::string& fileName);

	void filtering_QA();
  void remove_primerID_collisions(int minC, double minMajorFraction, bool report = true, int minDisplay = 0);

  // RT stuff
  hamming_return_type calculate_RT_mismatches() const;
  double LogLik(double s, double r) const;

  // PCR stuff
  d_hamming_return_type calculate_PCR_mismatches() const;

  // variance/overdispersion stuff
  void display_raw_and_primerID_counts() const;

  // chao estimator
  void calculate_effective_RNA_number() const;

  // diagnostics
  void show_recombination_patterns() const;
  void show_primerIDs_with_min_coverage(size_t minC) const;
  void show_clone_frequencies() const;

  void write_consensus_to_fasta() const;
  void write_to_fasta() const;
  void write_raw_to_fasta() const;
};

class alignments
{
public:
  size_t _n;

  std::vector<alignment> collections_alignments;

  int _min_current_coverage;
  double _min_majority_fraction;

  bool is_collisions_removed;

  double calculate_RT_LogLik_recombination(double s, double r) const;
  double neg_LogLik_recombination(double s, double r) const;

public:
	// primerID statistics
	seq_statistics _pid_stats;
	
  double _RT_sub_rate;
  double _RT_sub_rate_CI_lower;
  double _RT_sub_rate_CI_upper;

  double _RT_recomb_rate;
  double _RT_recomb_rate_CI_lower;
  double _RT_recomb_rate_CI_upper;
  double _RT_LogLik;

  double _PCR_sub_rate;
  double _PCR_sub_rate_CI_lower;
  double _PCR_sub_rate_CI_upper;

  alignments(const std::vector<std::string>& inputFiles);

	void filtering_QA();
  void remove_primerID_collisions(int minC = 10, double minMajorFraction = 0.9,
                                  bool printOut = false);
  void calculate_effective_population_size();

  // RT stuff
  double estimate_RT_substitution_rate(bool report = false);

  double estimate_RT_recombination_rate(bool report = false);
  double estimate_RT_recombination_rate(double s, bool report = false);

  // PCR stuff
  double estimate_PCR_substitution_rate(bool report = false);
  double estimate_PCR_recombinant_fraction(bool report = false);

  // variance/overdispersion stuff
  void display_raw_and_primerID_counts() const;
	
	// bias
	void write_prob_to_csv();
	double calculate_RT_bias_pvalue() const;

  // chao estimator
  void estimate_effective_RNA_number(bool report = false) const;

  // diagnostics
  void show_recombination_patterns() const;
  void show_primerIDs_with_min_coverage(int minC = 100) const;
  void show_clone_frequencies() const;

  void plot_RT_recombination_LogLik(double lower = 1E-6, int n = 1000) const;

  void write_all_consensus_to_fasta() const;
  void write_all_to_fasta() const;
  void write_raw_to_fasta() const;
};

#endif /* _ALIGNMENT_HPP_ */