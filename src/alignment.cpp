#include <iostream>
#include <iomanip>
#include <fstream>

#include <alignment.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/math/tools/minima.hpp>
#include <boost/math/tools/roots.hpp>

// custom functors
#include <prob_cycle.hpp>
#include <kmeans.hpp>

// alignment
struct SAMentry
{
  std::string QNAME;
  int FLAG;
  std::string RNAME;
  int POS;
  int MAPQ;
  std::string CIGAR;
  std::string RNEXT;
  int PNEXT;
  int TLEN;
  std::string SEQ;
  std::string QUAL;

  SAMentry() = default;

  SAMentry(
      const std::string& _QNAME,
      const std::string& _FLAG,
      const std::string& _RNAME,
      const std::string& _POS,
      const std::string& _MAPQ,
      const std::string& _CIGAR,
      const std::string& _RNEXT,
      const std::string& _PNEXT,
      const std::string& _TLEN,
      const std::string& _SEQ,
      const std::string& _QUAL)
      : QNAME(_QNAME), FLAG(stoi(_FLAG)), RNAME(_RNAME), POS(stoi(_POS)), MAPQ(stoi(_MAPQ)), CIGAR(_CIGAR), RNEXT(_RNEXT), PNEXT(stoi(_PNEXT)), TLEN(stoi(_TLEN)), SEQ(_SEQ), QUAL(_SEQ)
  {
  }
};

typedef std::pair<std::string, std::pair<SAMentry, SAMentry>> raw_sequence_pair;
typedef std::map<std::string, std::pair<SAMentry, SAMentry>> raw_sequence_pairs;

bool construct_sequence(const SAMentry& READ, std::string& DNA, std::string& PrimerID, const reference& ref)
{
  const std::string& read_sequence = READ.SEQ;
  const std::string& CIGAR = READ.CIGAR;

  int curPos = 0;
  int curOffset = 0;
  int offsetRef = READ.POS;
  int curLen;
  int lengthPID = 0;

  bool valid = true;

  DNA.clear();
  PrimerID.clear();

  // 1.) tokenize CIGAR
  for (int i = 0; i < CIGAR.length(); ++i)
  {
    if (isalpha(CIGAR[i]))
    {
      curLen = stoi(CIGAR.substr(curPos, i - curPos));
      curPos = i + 1;

      switch (CIGAR[i])
      {
      case 'M':
        DNA.append(read_sequence, curOffset, curLen);
        curOffset += curLen;
        offsetRef += curLen;
        break;

      case 'D':
        DNA.append(curLen, 'N');
        offsetRef += curLen;
        break;

      case 'I':
        if (offsetRef == ref.PID_start)
        {
          valid = false;
          lengthPID += curLen;
        }

        curOffset += curLen;
        break;

      case 'S':
        curOffset += curLen;
        break;
      }
    }
  }

  // 2.) check all strings
  if (READ.POS > 260)
  {
    // right sequence
    if (offsetRef - ref.overhang_start >= ref.overhang_min_length)
    {
      // overhang is good
      PrimerID = DNA.substr(ref.PID_start - READ.POS, ref.PID_length);
      lengthPID += std::count(PrimerID.begin(), PrimerID.end(), 'N');

      if (lengthPID > 0)
      {
        // contains some N's, PrimerID shorter than 10 nucleotides
        valid = false;
      }
    }
    else
    {
      // overhang too short
      valid = false;
    }

    if (READ.POS + DNA.length() > ref.replace_start)
      DNA.erase(ref.replace_start - READ.POS);

    DNA.append(ref.PID_start - READ.POS - DNA.length(), 'N');
  }

  return valid;
}

void merge_reads(const std::string& left_Read, const std::string& right_Read, int left_start, int right_start, const reference& ref, std::string& new_final_read)
{
  int global_start, global_end;

  global_start = left_start;
  global_end = right_start + right_Read.length();

  new_final_read.clear();
  new_final_read.append(std::max(left_start - 1, 0), 'N');
  new_final_read.append(left_Read);
  new_final_read.append(std::max(static_cast<int>(right_start - 1 - new_final_read.length()), 0), 'N');
  new_final_read.append(right_Read);
  new_final_read.append(std::max(static_cast<int>(ref.genome_length - new_final_read.length()), 0), 'N');
}

alignment::alignment(const std::string& fileName)
    : input_fileName(fileName)
{
  int slash_pos = input_fileName.find_last_of('/');
  slash_pos = (slash_pos == std::string::npos ? 0 : slash_pos + 1);
  just_fileName = input_fileName.substr(slash_pos);

  if (input_fileName.find("3223") != std::string::npos)
    this->_reference = ref3223;
  else
    this->_reference = ref3236;

  // 1.) read raw SAM input
  raw_sequence_pairs sam_data;

  std::ifstream input;
  input.open(fileName.c_str());
  raw_sequence_pairs::iterator it;

  if (input.is_open())
  {
    std::string temp, QNAME;
    std::vector<std::string> SplitVec;
    int line_no = 0;

    while (input.good())
    {
      ++line_no;
      if (line_no % 100000 == 0)
        std::cout << "Loaded " << line_no << " lines\n";

      getline(input, temp);
      temp.erase(std::remove(temp.begin(), temp.end(), '\n'), temp.end());

      if ((temp.length() == 0) || (temp[0] == '@'))
        continue;

      boost::split(SplitVec, temp, boost::is_any_of("\t"), boost::token_compress_on);

      QNAME = SplitVec[0].substr(0, SplitVec[0].find(' '));
      it = sam_data.find(QNAME);
      if (it == sam_data.end())
      {
        // not found
        sam_data.insert(raw_sequence_pair(QNAME, std::pair<SAMentry, SAMentry>(
                                                     SAMentry(
                                                         QNAME,
                                                         SplitVec[1],
                                                         SplitVec[2],
                                                         SplitVec[3],
                                                         SplitVec[4],
                                                         SplitVec[5],
                                                         SplitVec[6],
                                                         SplitVec[7],
                                                         SplitVec[8],
                                                         SplitVec[9],
                                                         SplitVec[10]),
                                                     SAMentry())));
      }
      else
      {
        // found
        it->second.second = SAMentry(
            QNAME,
            SplitVec[1],
            SplitVec[2],
            SplitVec[3],
            SplitVec[4],
            SplitVec[5],
            SplitVec[6],
            SplitVec[7],
            SplitVec[8],
            SplitVec[9],
            SplitVec[10]);

        // swap elements such that first element is left sequence and second element is right sequence
        if (it->second.first.POS > 260)
        {
          std::swap(it->second.first, it->second.second);
        }
      }
    }
  }
  else
  {
    std::cout << fileName << " is not readable!\n";
    exit(EXIT_FAILURE);
  }

  input.close();

  // 2.) convert raw SAM data to proper sequencing reads
  int accepted = 0;

  std::string primerID;
  std::string leftRead, rightRead;
  std::string temp_final_read;
  std::pair<std::map<std::string, std::vector<proper_read>>::iterator, bool> insert_it;
  bool valid;

  for (const auto& i : sam_data)
  {
    // 1.) check whether reads have a matching mate
    if ((i.second.second.SEQ.length() == 0) || (i.second.first.SEQ.length() == 0))
      continue;

    // 2.) parse left and right read and check whether primerID is valid
    construct_sequence(i.second.first, leftRead, primerID, _reference);
    valid = construct_sequence(i.second.second, rightRead, primerID, _reference);

    // 3.) merge left and right read
    merge_reads(leftRead, rightRead, i.second.first.POS, i.second.second.POS, _reference, temp_final_read);

    _reference.assign_counts(temp_final_read);

    if (!(valid))
      continue;

    // 4.) add to map
    insert_it = raw_primerID_map.emplace(primerID, std::vector<proper_read>());
    insert_it.first->second.emplace_back(temp_final_read, _reference);
    ++accepted;
  }

  std::cout << input_fileName << ": " << sam_data.size() << " pairs, " << accepted << " passed initial QA (" << std::fixed << std::setprecision(1) << static_cast<double>(accepted) / sam_data.size() * 100 << "%)\n";
}

/*
void __verbose(const std::string& PrimerID, const std::vector<proper_read>& reads, int valid_reads)
{
	std::cout << PrimerID << " (" << reads.size() << ", valid: " << valid_reads << ") - \n";
	for (const proper_read& i : reads)
	{
		std::cout << i.heterozygous_loci_string;
		std::cout << '\t';
		if (i.hamming_distance_to_best_reference >= 2)
			std::cout << "RECOMBINANT";
		else
			std::cout << i._ref.all_reference_strains[i.best_reference].name;

		std::cout << " (" << i.hamming_distance_to_best_reference << " / " << i.no_of_valid_heterozygous_bases << ")\n";
	}

	std::cout << "\n";
}
*/

std::string alignment::call_consensus_and_remove_collisions(const std::vector<proper_read>& reads, int minDisplay, const std::string& PrimerID)
{
  static kmeans clustering;

  if (reads.size() < _min_current_coverage)
  {
    // indecisive - coverage too low
    return std::string("1");
  }

  clustering(reads, _min_majority_fraction);
  if (clustering._frac1 < _min_majority_fraction)
  {
    // collision
    return std::string("0");
  }
  else
  {
    // no collision
    return clustering._full_consensus1;
  }
}

void alignment::remove_primerID_collisions(int minC, double minMajorFraction, bool report, int minDisplay)
{
  _min_current_coverage = minC;
  _min_majority_fraction = minMajorFraction;

  std::string consensus;
  number_collisions = 0;
  number_indecisive = 0;
  number_singletons = 0;

  collision_free_primerID_map.clear();
  consensus_primerID_map.clear();

  std::pair<std::map<std::string, consensus_read>::iterator, bool> it;

  for (auto& i : raw_primerID_map)
  {
    consensus = call_consensus_and_remove_collisions(i.second, minDisplay, i.first);

    switch (consensus[0])
    {
    case '0':
      // collision
      ++number_collisions;
      break;

    case '1':
      // indecisive
      ++number_indecisive;
      break;

    default:
      // OK, proper consensus
      ++number_singletons;

      _reference.assign_counts(consensus, true);

      collision_free_primerID_map.emplace(i.first, &(i.second));
      it = consensus_primerID_map.emplace(std::piecewise_construct,
                                          std::forward_as_tuple(i.first),
                                          std::forward_as_tuple(consensus, _reference, i.second.size()));
      break;
    }
  }

  _reference.normalise_counts();

  if (report)
  {
    std::cout << "Indecisive PrimerIDs: " << number_indecisive << '\n';
    std::cout << " Collision PrimerIDs: " << number_collisions << " (" << std::fixed << std::setprecision(1) << static_cast<double>(number_collisions) / (number_singletons + number_collisions) * 100 << "%)\n";
    std::cout << "     Total PrimerIDs: " << number_singletons + number_collisions << '\n';
  }
}

void alignment::show_primerIDs_with_min_coverage(int minC) const
{
  for (const auto& i : raw_primerID_map)
  {
    if (i.second.size() >= minC)
    {
      std::cout << i.first << '\n' << std::string(10, '-') << '\n';
      for (const auto& j : i.second)
        std::cout << j.heterozygous_loci_string << '\n';

      std::cout << '\n';
    }
  }
}

std::tuple<uint64_t, uint64_t, uint64_t> alignment::calculate_RT_mismatches() const
{
  int mismatches = 0;
  int valid_trials = 0;
  int Ns = 0;
  std::tuple<uint64_t, uint64_t, uint64_t> temp;

  for (const std::pair<const std::string, consensus_read>& i : consensus_primerID_map)
  {
    temp = i.second.calculate_homozygous_mismatches();

    mismatches += std::get<0>(temp);
    valid_trials += std::get<1>(temp);
    Ns += std::get<2>(temp);
  }

  return std::tuple<uint64_t, uint64_t, uint64_t>(mismatches, valid_trials, Ns);
}

double alignment::LogLik(double s, double r) const
{
  double total_likelihood = 0;
  double temp;

  if (_reference.freq_initialised == false)
  {
    std::cout << "Cannot calculate Log Likelihood without estimated reference frequencies!\n";
    exit(EXIT_FAILURE);
  }

  std::map<std::string, double> string_to_prob_cache;
  std::map<std::string, double>::const_iterator it;

  for (const auto& i : consensus_primerID_map)
  {
    it = string_to_prob_cache.find(i.second.heterozygous_loci_string);
    if (it == string_to_prob_cache.end())
    {
      // new sequence
      temp = i.second.log_prob(s, r);
      string_to_prob_cache.emplace(i.second.heterozygous_loci_string, temp);
      total_likelihood += temp;

      // std::cout << i.second.heterozygous_loci_string << '\t' << temp << '\n';
    }
    else
    {
      total_likelihood += it->second;
    }
  }

  return total_likelihood;
}

std::tuple<double, uint64_t, uint64_t> alignment::calculate_PCR_mismatches() const
{
  double mismatches = 0;
  uint64_t trials = 0;
  uint64_t numberN = 0;

  std::map<char, int> current_loci;

  // bool PCRerror;
  // int where;

  int cov;
  // double fracMajor;
  ranked_DNA_list ranker;
  double p;

  for (std::map<std::string, std::vector<proper_read>*>::const_iterator i = collision_free_primerID_map.begin(); i != collision_free_primerID_map.end(); ++i)
  {
    cov = i->second->size();
    // PCRerror = false;

    for (int j = 0; j < _reference.no_homozygous_loci; ++j)
    {
      // 1.) perform pileup
      current_loci = { { 'A', 0 }, { 'C', 0 }, { 'G', 0 }, { 'T', 0 }, { 'N', 0 } };

      for (int k = 0; k < cov; ++k)
      {
        ++current_loci[(*(i->second))[k].homozygous_loci_string[j]];
      }

      // 2.) perform ranking
      ranker(current_loci);

      if (ranker.no_N >= 2)
        continue;

      ++trials;

      if (ranker.no_mt <= 2)
        continue;

      p = PCR_prob(ranker.no_mt, ranker.no_mt + ranker.no_wt, 2);
      mismatches += p;

      /*
			   if (p > 0.995)
			   {
			        PCRerror = true;
			        where = j;
			   }
			 */
    }

    /*
		   if ((PCRerror) && (cov > 10))
		   {
		        std::cout << i->first << " (" << where << ")\n";

		        for (std::vector<proper_read>::const_iterator x = i->second->begin(); x != i->second->end(); ++x)
		        {
		                std::cout << x->homozygous_loci_string << '\n';
		        }
		        std::cout << '\n';
		   }
		 */
  }

  return std::tuple<double, uint64_t, uint64_t>(mismatches, trials, numberN);
}

void alignment::calculate_effective_RNA_number() const
{
  int eff_N = raw_primerID_map.size() + round(0.5 * number_singletons * number_singletons / number_collisions);

  std::cout << "        Number RNAs: " << eff_N;
}

void alignment::show_recombination_patterns() const
{
  for (const auto& i : consensus_primerID_map)
  {
    if (i.second.best_reference == i.second._ref.K)
    {
      // recombinant
      std::cout << i.first << " (" << i.second.multiplicity << ")\n";

      for (const auto& j : i.second._ref.all_reference_strains)
      {
        std::cout << j.name << '\t';

        for (int k = 0; k < i.second._ref.no_heterozygous_loci; ++k)
        {
          if ((i.second.heterozygous_loci_string[k] == 'N') || (i.second.heterozygous_loci_string[k] == j.heterozygous_loci_string[k]))
            std::cout << i.second.heterozygous_loci_string[k];
          else
            std::cout << ' ';
        }

        std::cout << '\n';
      }

      std::cout << '\n';
    }
  }
}

void alignment::display_raw_and_primerID_counts() const
{
  for (int i = 0; i < _reference.K; ++i)
  {
    std::cout << _reference.all_reference_strains[i].name << '\t';
    std::cout << _reference.all_reference_strains[i].PID_counts << '\t';
    std::cout << _reference.all_reference_strains[i].raw_counts << '\n';
  }

  std::cout << "-----------------------\n";
  std::cout << '\t' << _reference.PID_total << '\t' << _reference.raw_total;
}

void alignment::show_clone_frequencies() const
{
  _reference.display_strains_hetero();
}

void alignment::write_consensus_to_fasta() const
{
  int startPos = input_fileName.find_last_of('/') + 1;
  int endPos = input_fileName.find_last_of(".sam") - 3;

  std::string output_fileName = input_fileName.substr(startPos, endPos - startPos) + "_cons.fasta";

  std::ofstream output(output_fileName.c_str());
  for (const auto& i : consensus_primerID_map)
  {
    output << '>' << i.first << '_' << i.second.multiplicity << '_' << (i.second.best_reference == i.second._ref.K ? "RECOMB" : i.second._ref.all_reference_strains[i.second.best_reference].name) << '_' << i.second.hamming_distance_to_best_reference << '\n' << i.second.fullRead << '\n';
  }

  output.close();
}

void alignment::write_to_fasta() const
{
  int startPos = input_fileName.find_last_of('/') + 1;
  int endPos = input_fileName.find_last_of(".sam") - 3;

  std::string output_fileName = input_fileName.substr(startPos, endPos - startPos) + ".fasta";

  std::ofstream output(output_fileName.c_str());
  for (const auto& i : collision_free_primerID_map)
  {
    int k = 0;
    for (const proper_read& j : *(i.second))
    {
      ++k;
      output << '>' << i.first << '_' << k << '\n' << j.fullRead << '\n';
    }
  }

  output.close();
}

void alignment::write_raw_to_fasta() const
{
  int startPos = input_fileName.find_last_of('/') + 1;
  int endPos = input_fileName.find_last_of(".sam") - 3;

  std::string output_fileName = input_fileName.substr(startPos, endPos - startPos) + ".fasta";

  std::ofstream output(output_fileName.c_str());
  for (const auto& i : raw_primerID_map)
  {
    int k = 0;
    for (const proper_read& j : i.second)
    {
      ++k;
      output << '>' << i.first << '_' << k << '\n' << j.fullRead << '\n';
    }
  }

  output.close();
}

// alignments
alignments::alignments(const std::vector<std::string>& inputFiles)
    : _n(inputFiles.size()), _min_current_coverage(-1), _min_majority_fraction(-1), is_collisions_removed(false)
{
  collections_alignments.reserve(_n);

  for (const std::string& i : inputFiles)
  {
    std::cout << "Loading " << i << '\n';
    collections_alignments.emplace_back(i);
  }
}

void alignments::remove_primerID_collisions(int minC, double minMajorFraction, bool printOut)
{
  _min_current_coverage = minC;
  _min_majority_fraction = minMajorFraction;

  for (alignment& i : collections_alignments)
  {
    i.remove_primerID_collisions(_min_current_coverage, _min_majority_fraction, printOut, 40);
  }
}

/* *** RT stuff *** */
double alignments::estimate_RT_substitution_rate(bool report)
{
  uint64_t mismatches = 0;
  uint64_t valid_trials = 0;
  uint64_t Ns = 0;
  std::tuple<uint64_t, uint64_t, uint64_t> temp;

  for (const alignment& i : collections_alignments)
  {
    temp = i.calculate_RT_mismatches();

    mismatches += std::get<0>(temp);
    valid_trials += std::get<1>(temp);
    Ns += std::get<2>(temp);
  }

  _RT_sub_rate = static_cast<double>(mismatches) / valid_trials;
  _RT_sub_rate_CI_lower = boost::math::binomial_distribution<>::find_lower_bound_on_p(valid_trials, mismatches, 0.025);
  _RT_sub_rate_CI_upper = boost::math::binomial_distribution<>::find_upper_bound_on_p(valid_trials, mismatches, 0.025);

  if (report)
  {
    std::cout << std::string(50, '=') << '\n';
    std::cout << "RT substitution rate" << '\n';
    std::cout << "-----------------------\n";
    std::cout << '\n';

    std::cout << "           mt bases: " << mismatches << '\n';
    std::cout << "        Total bases: " << valid_trials << '\n';
    std::cout << "          'N' bases: " << Ns << '\n';
    std::cout << '\n';

    std::cout << "   Est. RT sub rate: " << std::scientific << std::setprecision(2) << _RT_sub_rate << '\n';
    std::cout << "             95% CI: [" << _RT_sub_rate_CI_lower << ", " << _RT_sub_rate_CI_upper << "]\n";

    std::cout << std::string(50, '=') << '\n';
  }

  return _RT_sub_rate;
}

double alignments::calculate_RT_LogLik_recombination(double s, double r) const
{
  double total_likelihood = 0;

  for (const alignment& i : collections_alignments)
    total_likelihood += i.LogLik(s, r);

  return total_likelihood;
}

double alignments::neg_LogLik_recombination(double s, double r) const
{
  return -calculate_RT_LogLik_recombination(s, r);
}

double alignments::estimate_RT_recombination_rate(bool report)
{
  return estimate_RT_recombination_rate(_RT_sub_rate, report);
}

double alignments::estimate_RT_recombination_rate(double s, bool report)
{
  std::pair<double, double> result;

  // 1.) calculate MLE
  result = boost::math::tools::brent_find_minima([&](double r)->double
  {return this->neg_LogLik_recombination(s, r);
                                                 },
                                                 1E-10, 0.9, 1000);

  _RT_recomb_rate = result.first;
  _RT_LogLik = -result.second;

  // 2.) calculate 95% CI
  uintmax_t max_iter = 1000;
  boost::math::tools::eps_tolerance<double> tol(10000);

  // lower CI boundary
  max_iter = 1000;
  result = boost::math::tools::toms748_solve([&](double r)->double
  {return this->calculate_RT_LogLik_recombination(s, r) - _RT_LogLik + 1.920729;
                                             },
                                             _RT_recomb_rate / 50,
                                             _RT_recomb_rate,
                                             tol,
                                             max_iter);

  _RT_recomb_rate_CI_lower = (result.first + result.second) / 2;

  // upper CI boundary
  max_iter = 1000;
  result = boost::math::tools::toms748_solve([&](double r)->double
  {return this->calculate_RT_LogLik_recombination(s, r) - _RT_LogLik + 1.920729;
                                             },
                                             _RT_recomb_rate,
                                             _RT_recomb_rate * 50,
                                             tol,
                                             max_iter);

  _RT_recomb_rate_CI_upper = (result.first + result.second) / 2;

  if (report)
  {
    // 3.) display
    std::cout << std::string(50, '=') << '\n';
    std::cout << "RT recombination rate" << '\n';
    std::cout << "-----------------------\n";
    std::cout << '\n';

    std::cout << "Est. RT recomb rate: " << std::scientific << std::setprecision(2) << _RT_recomb_rate << '\n';
    std::cout << "     Log-Likelihood: " << std::fixed << std::setprecision(2) << _RT_LogLik << '\n';
    std::cout << "             95% CI: [" << std::scientific << std::setprecision(2) << _RT_recomb_rate_CI_lower << ", " << _RT_recomb_rate_CI_upper << "]\n";

    std::cout << std::string(50, '=') << '\n';
  }

  return _RT_recomb_rate;
}

/* *** PCR stuff *** */
double alignments::estimate_PCR_substitution_rate(bool report)
{
  double mismatches = 0;
  uint64_t valid_trials = 0;
  uint64_t Ns = 0;
  std::tuple<double, uint64_t, uint64_t> temp;

  for (const alignment& i : collections_alignments)
  {
    temp = i.calculate_PCR_mismatches();

    mismatches += std::get<0>(temp);
    valid_trials += std::get<1>(temp);
    Ns += std::get<2>(temp);
  }

  _PCR_sub_rate = mismatches / valid_trials;
  _PCR_sub_rate_CI_lower = boost::math::binomial_distribution<>::find_lower_bound_on_p(valid_trials, mismatches, 0.025);
  _PCR_sub_rate_CI_upper = boost::math::binomial_distribution<>::find_upper_bound_on_p(valid_trials, mismatches, 0.025);

  if (report)
  {
    std::cout << std::string(50, '=') << '\n';
    std::cout << "PCR substitution rate" << '\n';
    std::cout << "-----------------------\n";
    std::cout << '\n';

    std::cout << "           mt bases: " << mismatches << '\n';
    std::cout << "        Total bases: " << valid_trials << '\n';
    std::cout << "          'N' bases: " << Ns << '\n';
    std::cout << '\n';

    std::cout << "  Est. PCR sub rate: " << std::scientific << std::setprecision(2) << _PCR_sub_rate << '\n';
    std::cout << "             95% CI: [" << _PCR_sub_rate_CI_lower << ", " << _PCR_sub_rate_CI_upper << "]\n";

    std::cout << std::string(50, '=') << '\n';
  }

  return _PCR_sub_rate;
}

/* *** variance/overdispersion stuff *** */
void alignments::display_raw_and_primerID_counts() const
{
  std::cout << std::string(50, '=') << '\n';
  std::cout << "Overdispersion/Variance" << '\n';
  std::cout << "-----------------------";

  for (const alignment& i : collections_alignments)
  {
    std::cout << '\n' << i.just_fileName << '\n';
    i.display_raw_and_primerID_counts();
    std::cout << '\n';
  }

  std::cout << std::string(50, '=') << '\n';
}

/* *** chao estimator *** */
void alignments::estimate_effective_RNA_number(bool report) const
{
  if (report)
  {
    std::cout << std::string(50, '=') << '\n';
    std::cout << "RNA numbers" << '\n';
    std::cout << "-----------------------";

    for (const alignment& i : collections_alignments)
    {
      std::cout << '\n' << i.just_fileName << '\n';
      i.calculate_effective_RNA_number();
      std::cout << '\n';
    }

    std::cout << std::string(50, '=') << '\n';
  }
}

/* *** diagnostics *** */
void alignments::show_recombination_patterns() const
{
  for (const alignment& i : collections_alignments)
  {
    i.show_recombination_patterns();
  }
}

void alignments::show_primerIDs_with_min_coverage(int minC) const
{
  for (const alignment& i : collections_alignments)
  {
    std::cout << i.input_fileName << '\n';

    i.show_primerIDs_with_min_coverage(minC);
  }
}

void alignments::show_clone_frequencies() const
{
  for (const alignment& i : collections_alignments)
  {
    i.show_clone_frequencies();
  }
}

void alignments::plot_RT_recombination_LogLik(double lower, int n) const
{
  std::pair<double, double> result;
  double CI_LogLik = calculate_RT_LogLik_recombination(_RT_sub_rate, lower);

  uintmax_t max_iter = 1000;
  boost::math::tools::eps_tolerance<double> tol(10000);

  // upper plot boundary
  result = boost::math::tools::toms748_solve([&](double r)->double
  {return this->calculate_RT_LogLik_recombination(_RT_sub_rate, r) - CI_LogLik;
                                             },
                                             _RT_recomb_rate,
                                             _RT_recomb_rate * 100,
                                             tol,
                                             max_iter);

  double upper = (result.first + result.second) / 2;

  std::cout << "Drawing Plot from " << lower << " to " << upper << "\n";
  double factor = 1.0 / (n - 1) * (log(upper) - log(lower));
  factor = exp(factor);

  std::vector<double> x;
  x.reserve(n + 10);

  std::vector<double> y;
  y.reserve(n + 10);

  double temp;

  int I = 1;
  for (double i = lower; i <= upper; i *= factor, ++I)
  {
    x.emplace_back(i);

    temp = calculate_RT_LogLik_recombination(_RT_sub_rate, i);
    y.emplace_back(temp);

    if (I % 50 == 0)
      std::cout << "Iteration " << std::fixed << std::setprecision(0) << I << "\tr: " << std::scientific << std::setprecision(4) << i << "\tLogLik: " << temp << '\n';
  }

  // write R file
  std::ofstream Routput("LogLik.R");

  Routput << "x = c(" << std::scientific << std::setprecision(4) << x[0];
  for (int i = 1; i < x.size(); ++i)
    Routput << ", " << x[i];

  Routput << ")\n";

  Routput << "y = c(" << std::fixed << std::setprecision(4) << y[0];
  for (int i = 1; i < y.size(); ++i)
    Routput << ", " << y[i];

  Routput << ")\n\n";

  Routput << "require(\"sfsmisc\")\n";
  Routput << "require(\"extrafont\")\n";
  Routput << "loadfonts()\n";
  Routput << "\n";
  Routput << "pdf(\"tmp.pdf\", family=\"CM Roman\")\n";
  Routput << "par(mar = c(3, 5.3, 3, 1))\n";
  Routput << "par(xpd = TRUE)\n";
  Routput << "\n";
  Routput << "plot(x, y, type=\"n\", log=\"x\", main = \"\", axes=FALSE, xlab=\"\", ylab=\"\")\n";
  Routput << "lines(x, y, type=\"l\")\n";
  Routput << "\n";
  Routput << "mtext(\"RT HMM LogLikelihood\", side=3, line = 1, cex = 1.8, font = 2)\n";
  Routput << "\n";
  Routput << "# x-axis\n";
  Routput << "eaxis(1)\n";
  Routput << "mtext(expression(italic(r)), side=1, line = 2, cex = 1.5)\n";
  Routput << "\n";
  Routput << "# y-axis\n";
  Routput << "eaxis(2)\n";
  Routput << "mtext(\"LogLik\", side=2, line = 4, cex = 1.5)\n";
  Routput << "\n";
  Routput << "YMIN = par(\"usr\")[3]\n";
  Routput << "YMAX = par(\"usr\")[4]\n";
  Routput << "\n";
  Routput << "LOGLIKMAX = " << std::fixed << std::setprecision(4) << _RT_LogLik << " -1.920729\n";
  Routput << "X_CI_LOW = " << std::scientific << std::setprecision(4) << _RT_recomb_rate_CI_lower << "\n";
  Routput << "X_CI_HIGH = " << std::scientific << std::setprecision(4) << _RT_recomb_rate_CI_upper << "\n";
  Routput << "r = " << std::scientific << std::setprecision(2) << _RT_recomb_rate << "\n";
  Routput << "\n";
  Routput << "MIDP = sqrt(X_CI_LOW*X_CI_HIGH)\n";
  Routput << "\n";
  Routput << "segments(X_CI_LOW, YMIN, X_CI_LOW, LOGLIKMAX)\n";
  Routput << "segments(X_CI_HIGH, YMIN, X_CI_HIGH, LOGLIKMAX)\n";
  Routput << "\n";
  Routput << "segments(r, YMIN, r, YMIN + 0.93*(LOGLIKMAX-YMIN), lty=3)\n";
  Routput << "segments(r, YMIN + 0.97*(LOGLIKMAX-YMIN), r, YMIN + 1.05*(LOGLIKMAX-YMIN), lty=3)\n";
  Routput << "\n";
  Routput << "EXP = floor(log10(r))\n";
  Routput << "BASE = signif(r, 3) / 10^EXP\n";
  Routput << "\n";
  Routput << "text(r, YMIN + 1.05*(LOGLIKMAX-YMIN), bquote(hat(italic(r)) == .(BASE) %*% 10^.(EXP)), pos = 3, cex = 1.0)\n";
  Routput << "\n";
  Routput << "arrows(X_CI_LOW, YMIN + 0.95*(LOGLIKMAX-YMIN), MIDP/1.28, YMIN + 0.95*(LOGLIKMAX-YMIN), code = 1, length = 0.07)\n";
  Routput << "arrows(MIDP*1.28, YMIN + 0.95*(LOGLIKMAX-YMIN), X_CI_HIGH, YMIN + 0.95*(LOGLIKMAX-YMIN), code = 2, length = 0.07)\n";
  Routput << "\n";
  Routput << "text(x= MIDP, y=YMIN + 0.95*(LOGLIKMAX-YMIN), labels=\"95% CI\")\n";
  Routput << "dev.off()\n";
  Routput << "\n";
  Routput << "embed_fonts(\"tmp.pdf\", outfile=\"LogLik.pdf\")\n";
  Routput << "file.remove(\"tmp.pdf\")\n";

  Routput.close();
}

void alignments::write_all_consensus_to_fasta() const
{
  for (const alignment& i : collections_alignments)
  {
    i.write_consensus_to_fasta();
  }
}

void alignments::write_all_to_fasta() const
{
  for (const alignment& i : collections_alignments)
  {
    i.write_to_fasta();
  }
}

void alignments::write_raw_to_fasta() const
{
  for (const alignment& i : collections_alignments)
  {
    i.write_raw_to_fasta();
  }
}