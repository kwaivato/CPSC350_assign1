// logger.h
// A declared program used to read and write files for a DNA reader program.
// @author Vidal M. Arroyo
// @author arroy118@mail.chapman.edu
// @version 1.0

#include <fstream>

#include "geneCalculator.h"
#include "genePool.h"

using namespace std;

// Used to read in DNA-read files, compute statistical features, generate 1000 DNA strings from these features, and output
// this in a results file. Must use with the geneCalculator, genePool, and boxMuller classes.
class Logger
{
  public:
    Logger();
    ~Logger();
    // A gene calculator was used to determine the statistical parameters of a given DNA-read file.
    GeneCalculator getGeneCalculator();
    // A gene pool was used to generate DNA reads from a gaussian distribution based on the statistical parameters of
    // a given DNA-read file.
    GenePool getGenePool();
    void setGeneCalculator();
    void setGenePool();
    // The Out method is used to create a results file and write the statistical features + generated DNA reads from the
    // first file read-in.
    void Out(string f);
    // The AppendOut method is used to append results to the existing results file and write the statistical features +
    // generated DNA reads from the n + 1 file read-ins, given that n is any integer greater than 0.
    void AppendOut(string f);
  private:
    ifstream m_input_stream;
    ofstream m_output_stream;
    // The m_character variable is used to hold characters when performing string operations in the logger class.
    char m_character;
    // The m_single_line variable is used to hold a string for each line in a given DNA-read file.
    string m_single_line;
    // The m_single_line_cleaned variable is used to hold a 'cleaned' string for each line in a given DNA-read file, where
    // cleaned is defined as uppercasing each variable and replacing non-ACGT variables (and new line characters) with x's.
    string m_single_line_cleaned;
    // The m_entire_file_read variable is a long string that holds all of a DNA-read file's DNA reads and documents
    // spaces between reads with x's.
    string m_entire_file_read;
    GeneCalculator m_gene_calculator;
    GenePool m_gene_pool;
};
