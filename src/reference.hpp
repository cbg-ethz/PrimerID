#ifndef _REFERENCE_HPP_
#define _REFERENCE_HPP_

#include <vector>
#include <string>

enum referenceType
{
  START_3223,
  START_3236
};

typedef std::tuple<uint64_t, uint64_t, uint64_t> hamming_return_type;
hamming_return_type hamming_distance(const std::string& sA, const std::string& sB);

struct record
{
  std::string DNA;
  std::string name;

  std::string heterozygous_loci_string;

  // statistics:
  int raw_counts;
  int PID_counts;

  double frequency;

  record(const std::string& _DNA, const std::string& _name);
};

class reference
{
public:
  // global parameters
  static const int min_valid = 22;
  static const int max_mismatches = 0;

  // member variables
  std::string referenceFile;
  std::string just_fileName;

  referenceType reference_variant;

  int K;
  int genome_length;

  int raw_total;
  int PID_total;

  std::vector<record> all_reference_strains;
  bool freq_initialised;

  // loci information
  std::vector<int> included_loci;
  std::vector<int> not_included_loci;

  // heterozygous information
  int no_heterozygous_loci;
  std::vector<int> heterozygous_loci;

  // homozygous information
  int no_homozygous_loci;
  std::vector<int> homozygous_loci;
  std::string homozygous_string;

  // information for truncating
  int replace_start;
  int replace_length;

  int PID_start;
  int PID_length;

  int overhang_start;
  int overhang_min_length;

  // member functions:
  reference(const std::string& fileName, referenceType variant);
  reference() = default;

  void display_hamming_distance() const;
  void display_strains_hetero() const;
  void display_strains_verbose() const;

  void assign_counts(const std::string& read, bool consensus = false);
  void normalise_counts();
  void reset_consensus_counts();
};

extern reference ref3223;
extern reference ref3236;

#endif /* _REFERENCE_HPP_ */