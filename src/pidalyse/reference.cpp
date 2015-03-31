#include <algorithm>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>

#include "reference.hpp"

// record
record::record(const std::string& _DNA, const std::string& _name)
    : DNA(_DNA),
      name(_name),
      raw_counts(0),
      PID_counts(0),
      frequency(-1) {};

// reference
reference::reference(const std::string& fileName, referenceType variant)
    : referenceFile(fileName),
      reference_variant(variant),
      raw_total(0),
      PID_total(0),
      freq_initialised(false),
      no_heterozygous_loci(0),
      no_homozygous_loci(0)
{
  size_t slash_pos = referenceFile.find_last_of('/');
  slash_pos = (slash_pos == std::string::npos ? 0 : slash_pos + 1);
  just_fileName = referenceFile.substr(slash_pos);

  // 1.) load FASTA file
  std::ifstream input;
  input.open(referenceFile.c_str());

  if (input.is_open())
  {
    std::string temp, id;

    while (input.good())
    {
      getline(input, temp);
      temp.erase(std::remove(temp.begin(), temp.end(), '\n'), temp.end());

      if (temp.length() == 0)
        continue;

      if (temp[0] == '>')
      {
        // header
        id = temp.substr(1, std::string::npos);
      }
      else
      {
        // sequence
        all_reference_strains.emplace_back(record(temp, id));
      }
    }
  }
  else
  {
    std::cout << fileName << " is not readable!\n";
    exit(EXIT_FAILURE);
  }

  input.close();

  // 2.) initialise property variable
  // this depends on the locus under consideration
  replace_start = 496;
  if (reference_variant == START_3223)
    replace_length = 14;
  else
    replace_length = 27;

  PID_length = 10;
  PID_start = replace_start + replace_length;

  overhang_start = PID_start + PID_length;
  overhang_min_length = 23;

  // 3.) set loci to include
  // hard-coded at the moment, will be replaced at a later stage
	// 0-based indexing used here, makes usage with C++ easier
  if (reference_variant == START_3223)
  {
    included_loci = { 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246,
                      294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454, 455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494};
  }
  else
  {
    included_loci = { 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246,
                      312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454, 455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494};
  }

  // determine complement of included_loci
  for (int i = 0; i < replace_start - 1; ++i)
  {
    if (!std::binary_search(included_loci.begin(), included_loci.end(), i))
    {
      not_included_loci.push_back(i);
    }
  }

  // 4.) determine homozygous/heterozygous loci
  K = all_reference_strains.size();
  genome_length = all_reference_strains[0].DNA.length();

  char allele;
  bool homozygous;

  // decide whether locus i is homo or heterozygous
  for (int i : included_loci)
  {
    allele = all_reference_strains[0].DNA[i];
    homozygous = true;

    for (int j = 1; j < K; ++j)
    {
      homozygous = (allele == all_reference_strains[j].DNA[i]);
      if (!(homozygous))
        break;
    }

    if (homozygous)
    {
      homozygous_loci.push_back(i);
      homozygous_string.push_back(allele);
    }
    else
    {
      heterozygous_loci.push_back(i);
    }
  }

  no_heterozygous_loci = heterozygous_loci.size();
  no_homozygous_loci = homozygous_loci.size();

  // construct heterozygous strings
  for (record& i : all_reference_strains)
  {
    i.frequency = 0;

    for (int j : heterozygous_loci)
      i.heterozygous_loci_string.push_back(i.DNA[j]);
  }
}

void reference::display_strains_hetero() const
{
  std::cout << std::string(50, '=') << '\n';
  std::cout << "Reference strains" << '\n';
  std::cout << "-----------------------\n";

  std::cout << "Reference: " << just_fileName << '\n';
  std::cout << '\n';

  std::cout << "Heterozygous loci: ";

  // show positions of heterozygous loci
  for (int j : heterozygous_loci)
    std::cout << j << ' ';

  std::cout << '\n';

  // show heterozygous strings
  for (const record& i : all_reference_strains)
  {
    std::cout << i.name << '\t' << i.heterozygous_loci_string;
    if (i.frequency != 0)
      std::cout << '\t' << std::fixed << std::setprecision(1) << i.frequency * 100 << '\t' << i.PID_counts;

    std::cout << '\n';
  }

  if (PID_total)
  {
    std::cout << "Total PIDs used for frequency estimation: " << PID_total << '\n';
  }

  std::cout << "Number homozygous loci:   " << no_homozygous_loci << '\n';
  std::cout << "Number heterozygous loci: " << no_heterozygous_loci << '\n';
  std::cout << std::string(50, '=') << '\n';
}

void reference::display_strains_verbose() const
{
  const static int width = 60;

  std::cout << std::string(50, '=') << '\n';
  std::cout << "Reference strains (verbose)" << '\n';
  std::cout << "-----------------------\n";

  std::cout << "Reference: " << just_fileName << '\n';
  std::cout << '\n';

  for (int i = 0; i < genome_length; i += width)
  {
    // print 5VM DNA
    for (const record& j : all_reference_strains)
      std::cout << j.name << '\t' << j.DNA.substr(i, width) << '\n';

    // print annotation
    std::cout << std::string(7 - std::to_string(i + 1).length(), ' ') << i + 1;
    std::cout << '\t';
    for (int j = i; j < std::min(i + width, genome_length); ++j)
    {
      if (j < replace_start - 1)
      {
        // in main read
        if (std::binary_search(included_loci.begin(), included_loci.end(), j))
        {
          // found
          if (std::binary_search(heterozygous_loci.begin(), heterozygous_loci.end(), j))
          {
            // heterozygous locus
            std::cout << '*';
          }
          else
          {
            // homozygous locus
            std::cout << '.';
          }
        }
        else
        {
          // not found
          std::cout << ' ';
        }
      }
      else
      {
        // remainder, i.e. right primer
        std::cout << 'x';
      }
    }

    std::cout << ' ' << std::min(i + width, genome_length) << "\n\n";
  }

  std::cout << std::string(50, '=') << '\n';
}

void reference::assign_counts(const std::string& read, bool consensus)
{
  int min = 1000;
  int index = K;

  int hamming;
  int valid_trials;

  for (int i = 0; i < K; ++i)
  {
    hamming = 0;
    valid_trials = 0;

    for (int j : heterozygous_loci)
    {
      if (read[j] != 'N')
      {
        ++valid_trials;
        hamming += (read[j] != all_reference_strains[i].DNA[j]);
      }
    }

    // minimum number of valid bases to call a reference
    if (valid_trials < reference::min_valid)
      break;

    if (hamming < min)
    {
      min = hamming;
      index = i;
    }
  }

  if (index != K)
  {
    if (min <= reference::max_mismatches)
    {
      if (consensus)
        ++all_reference_strains[index].PID_counts;
      else
        ++all_reference_strains[index].raw_counts;
    }
  }
}

void reference::normalise_counts()
{
  PID_total = 0;
  raw_total = 0;

  for (const record& i : all_reference_strains)
  {
    PID_total += i.PID_counts;
    raw_total += i.raw_counts;
  }

  for (record& i : all_reference_strains)
    i.frequency = static_cast<double>(i.PID_counts) / PID_total;

  freq_initialised = true;
}

void reference::reset_consensus_counts()
{
  PID_total = 0;
  raw_total = 0;

  for (record& i : all_reference_strains)
    i.PID_counts = 0;

  freq_initialised = false;
}