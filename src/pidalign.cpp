/**
 * Copyright (c) 2014-2015 Jochen Singer, David Seifert
 *
 * This file is part of PIDalign
 *
 * PIDalign is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or any later version.
 *
 * PIDalign is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * PIDalign. If not, see <http://www.gnu.org/licenses/>.
 */

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/bam_io.h>
#include <seqan/modifier.h>

#include <src/ednafull.hpp>

#include <iostream>
#include <fstream>
#include <algorithm>

using namespace seqan;

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class AppOptions
// --------------------------------------------------------------------------

// This struct stores the options from the command line.
//
// You might want to rename this to reflect the name of your app.

struct AppOptions
{
    // The first (and only) argument of the program is stored here.
    StringSet<CharString> patternFileNames;
    CharString refFileName;
    CharString outputFileName;
    int numThreads;
    unsigned bufferSize;
    int minScore;

    AppOptions() :
        numThreads(1),
        bufferSize(100000),
        minScore(minValue<int>())
    {}
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

ArgumentParser::ParseResult
parseCommandLine(AppOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    ArgumentParser parser("pidalign");
    // Set short description, version, and date.
    setShortDescription(parser, "Needleman-Wunsch (Gotoh variant) aligner that produces a SAM output");
    setVersion(parser, "0.1");
    setDate(parser, "Februrary 2015");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fB<R1.fastq>\\fP [\\fB<R2.fastq>\\fP]");
    //addDescription(parser, "pid.");

    ArgParseArgument patternArg(ArgParseArgument::INPUT_FILE, "IN", true);
    setValidValues(patternArg, "fastq fq");
    addArgument(parser, patternArg);
		
    addOption(parser, ArgParseOption("r", "referenceFileName", "Name of the reference file.", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "r");
    setValidValues(parser, "r", "fasta fa fna");

    addOption(parser, ArgParseOption("o", "outputFileName", "Name of the output file.", ArgParseArgument::OUTPUT_FILE, "OUT"));
    setRequired(parser, "o");
    setValidValues(parser, "o", "sam");

    addOption(parser, ArgParseOption("t", "threads", "The number of threads to be used.", ArgParseArgument::INTEGER));
    addOption(parser, ArgParseOption("bf", "bufferSize", "The number of hits stored in a buffer before writing them to disk.", ArgParseArgument::INTEGER));
    addOption(parser, ArgParseOption("s", "minScore", "MinScore of alignment in order to be considered.", ArgParseArgument::INTEGER));

    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBpidalign\\fP \\fB-r\\fP \\fIreference.fasta\\fP\\fP \\fB-o\\fP \\fIoutput.sam\\fP \\fP\\fIR*.fastq\\fP",
                "Use \\fIreference.fasta\\fP to align \\fP\\fIR*.fastq\\fP with the SAM alignment in \\fIoutput.sam\\fP.");

    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    getArgumentValue(options.refFileName, parser, 0);
    resize(options.patternFileNames, getArgumentValueCount(parser, 1));
    for (unsigned i =0 ; i < getArgumentValueCount(parser, 1); ++i)
    {
        getArgumentValue(options.patternFileNames[i], parser, 1, i);
    }
    getOptionValue(options.outputFileName, parser, "o");
    getOptionValue(options.numThreads, parser, "th");
    getOptionValue(options.bufferSize, parser, "bf");
    getOptionValue(options.minScore, parser, "s");

    return ArgumentParser::PARSE_OK;
}

template <typename TStream, typename TAlign>
void writeLocal(TStream & stream,
           TAlign & align,
           String<char> const & patternId,
           String<char> const & refId,
           String<char> & cigar,
           int score)
{
    int clipBegin = row(align, 1)._array[0];
    int clipEnd = length(row(align, 1)) - row(align, 1)._array[length(row(align, 1)._array) - 1];
    setClippedBeginPosition(row(align, 0), clipBegin);
    setClippedBeginPosition(row(align, 1), clipBegin);
    setClippedEndPosition(row(align, 0), clipEnd);
    setClippedEndPosition(row(align, 1), clipEnd);
    getCigarString(cigar,row(align, 0), row(align, 1), 1000);
    stream << patternId << "\t";                                //QNAME
    stream << (unsigned short)0x0002 << "\t";                        //FLAG
    stream << refId << "\t";   //RNAME
    stream << row(align, 1)._array[0] + 1 << "\t";    //POS
    stream << score << "\t";
    stream << cigar << "\t";                                         //CIGAR
    stream << "*" << "\t";                                           //RNEXT
    stream << 0 << "\t";                                             //PNEXT
    stream << 0 << "\t";                                             //TLen
    stream << source(row(align, 1)) << "\t";                                //SEQ
    stream << "*" << "\n";                                           //QUAL

}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{
    typedef String<Iupac> TRefSeq;

    // Parse the command line.
    ArgumentParser parser;
    AppOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    std::cout << "primerIdIdentifyer\n"
              << "==================\n\n";

    // Preparation of out stream including SAM header
    std::ofstream stream(toCString(options.outputFileName), std::ios::out | std::ios::app);
    stream << "@HD\tVN:1.5\tSO:coordinate\n";
    // Reference
    SeqFileIn seqInRef(toCString(options.refFileName));
    StringSet<String<char> > refIds;
    StringSet<TRefSeq> refSeqs;
    readRecords(refIds, refSeqs, seqInRef);

    for (unsigned i = 0; i < length(refSeqs); ++i)
        stream << "@SQ\tSN:" << refIds[i] << "\tLN:" << length(refSeqs[i]) << "\n";
    stream.close();


    // Alignment helpers
    
    Ednafull scoringScheme(-1, -20);
    AlignConfig<true, false, false, true> alignConfig;

    unsigned counter = 0;
    for (unsigned fileCounter = 0; fileCounter < length(options.patternFileNames); ++fileCounter)
    {
        SeqFileIn seqInPattern(toCString(options.patternFileNames[fileCounter]));

        // This is the parallelization starting point
        SEQAN_OMP_PRAGMA(parallel num_threads(options.numThreads))
        for (; ;)
        {
            // Preparing the align Objects
            String<Align<TRefSeq> > alignObjs, alignObjsRevComp;
            resize(alignObjs, length(refSeqs));
            resize(alignObjsRevComp, length(refSeqs));

            int maxScore, maxScoreRevComp, maxId, maxIdRevComp;

            for (unsigned i = 0; i < length(refSeqs); ++i)
            {
                resize(rows(alignObjs[i]), 2);
                resize(rows(alignObjsRevComp[i]), 2);
                assignSource(row(alignObjs[i], 0), refSeqs[i]);
                assignSource(row(alignObjsRevComp[i], 0), refSeqs[i]);
            }

            // Pattern
            StringSet<CharString> patternIds;
            StringSet<String<Dna5> > patternSeqs;
            String<Dna5> revComPatternSeq;

            // Cigar String
            CharString cigar;

            std::ostringstream localStream;
            // Reading the pattern
            SEQAN_OMP_PRAGMA(critical (read_chunk))
            {
                if (!atEnd(seqInPattern))
                {
                    /*if(*/readRecords(patternIds, patternSeqs, seqInPattern, options.bufferSize);/* != 0)*/
                    //{
                    //    std::cout << "ERROR: Could not read samples!\n";
                    //}
                }
                counter++;
            }

            if (length(patternIds) == 0)
            {
                //std::cerr << "TEST" << std::endl;
                break;
            }

           // int lDiag = -30;
           // int uDiag = 30;
            for (unsigned i = 0; i < length(patternSeqs); ++i)
            {
                // Determine if right or left aligned
                String<Dna5> left= prefix(patternSeqs[i], 20);
                String<Dna5> right = suffix(patternSeqs[i], length(patternSeqs[i]) - 20);

                maxScore = minValue<int>();
                maxId = -1;
                for (unsigned j = 0; j < length(alignObjs); ++j)
                {
                    assignSource(row(alignObjs[j], 1), patternSeqs[i]);
                    int result = globalAlignment(alignObjs[j], scoringScheme, alignConfig, maxScore);//, lDiag, uDiag);
                    if (maxScore < result)
                    {
                        maxScore = result;
                        maxId = j;
                    }
                }
                revComPatternSeq = patternSeqs[i];
                reverseComplement(revComPatternSeq);
                maxScoreRevComp = minValue<int>();
                maxIdRevComp = -1;
                for (unsigned j = 0; j < length(alignObjsRevComp); ++j)
                {
                    assignSource(row(alignObjsRevComp[j], 1), revComPatternSeq);
                    int result = globalAlignment(alignObjsRevComp[j], scoringScheme, alignConfig, maxScoreRevComp);//, lDiag, uDiag);
                    if (maxScoreRevComp < result)
                    {
                        maxScoreRevComp = result;
                        maxIdRevComp = j;
                    }
                }
                if (maxScore > maxScoreRevComp)
                {
                    int normScore = maxScore / (length(revComPatternSeq) * 5.0) * 255.0;
                    if (normScore > options.minScore)
                        writeLocal(localStream, alignObjs[maxId], patternIds[i], refIds[maxId], cigar, normScore);
                }
                else
                {
                    int normScore = maxScoreRevComp / (length(revComPatternSeq) * 5.0) * 255.0;
                    if(normScore > options.minScore)
                        writeLocal(localStream, alignObjsRevComp[maxIdRevComp], patternIds[i], refIds[maxIdRevComp], cigar, normScore);
                }
            }
            SEQAN_OMP_PRAGMA(critical (write_chunk))
            {   //std::ostringstream testStream;
                //testStream << counter;
                //String<char> test = options.outputFileName;
                //append(test, testStream.str());
                //std::cerr << test << std::endl;
                std::ofstream stream(toCString(options.outputFileName), std::ios::out | std::ios::app);
                stream << localStream.str();
                stream.close();
            }
        }
    }

    return 0;
}


