// logger.cpp
// An implemented program used to read and write files for a DNA reader program.
// @author Vidal M. Arroyo
// @author arroy118@mail.chapman.edu
// @version 1.0

#include "logger.h"

using namespace std;

Logger::Logger()
{
  m_single_line = "";
  m_single_line_cleaned = "";
  m_entire_file_read = "";
  GeneCalculator m_gene_calculator;
  GenePool m_gene_pool;
}


Logger::~Logger()
{

}


GeneCalculator Logger::getGeneCalculator()
{
  return m_gene_calculator;
}

GenePool Logger::getGenePool()
{
  return m_gene_pool;
}

void Logger::setGeneCalculator()
{
  GeneCalculator m_gene_calculator;
}

void Logger::setGenePool()
{
  GenePool m_gene_pool;
}

// The Out method creates a results file titled 'VidalArroyo.txt' and appends both statistical results along with
// new DNA reads to this file.
void Logger::Out(string f)
{
  m_input_stream.open(f);
  // This if/else statement is used to catch any files that do not exist. If this isn't here, then the program will enter
  // an infinte loop and fail.
  if (m_input_stream.fail())
  {
    cout << "File does not exist!" << endl;
  }
  else
  {
    // The while loop while run as long as a DNA-read file has not been completely read.
    while(getline(m_input_stream, m_single_line))
    {
      // Since it is not guaranteed that the file is in all upper-caps or lower-caps, I use this for loop to make all
      // the base pairs upper case. This is in line with standard molecular biology notation.
      for (int i = 0; i < m_single_line.size(); ++i)
      {
        if (toupper(m_single_line.at(i)) == 'A')
        {
          // The m_single_line_cleaned string represents each line of DNA strings in the file, but now all the characters
          // are upper-case.
          m_single_line_cleaned += 'A';
        }
        else if (toupper(m_single_line.at(i)) == 'C')
        {
          m_single_line_cleaned += 'C';
        }
        else if (toupper(m_single_line.at(i)) == 'G')
        {
          m_single_line_cleaned += 'G';
        }
        else if (toupper(m_single_line.at(i)) == 'T')
        {
          m_single_line_cleaned += 'T';
        }
        // If a base pair is not present, then the character is skipped.
        else
        {
          continue;
        }
      }
      // The m_entire_file_read will hold all of the strings in the file for analysis. The x's are used to denote
      // miscellaneous characters and/or \n characters, which denote a new DNA read. These x's will allow us to
      // compute read length statistical parameters later on, since without them it is impossible for us to tell
      // when a given read begins or ends.
      m_entire_file_read += (m_single_line_cleaned + 'x');
      m_single_line_cleaned = "";
    }

    m_input_stream.close();

    // This line computes all of the statistical parameters that underlie the read lengths of a given DNA-read file.
    // See the geneCalculator class for implementation details.
    m_gene_calculator.ComputeStats(m_entire_file_read);

    m_output_stream.open("VidalArroyo.txt");

    // The following lines input my name and student ID to the top of the file.
    m_output_stream << "Name: Vidal Arroyo" << endl;
    m_output_stream << "Student ID: 2253413\n" << endl;

    // Now, the logger will print all statistical results to the results file.
    m_output_stream << "Here are the summary statistics of " << f << "\n" << endl;
    m_output_stream << "Sum (Σ) = " << m_gene_calculator.getSum() << endl;
    m_output_stream << "Mean (μ) = " << m_gene_calculator.getMean() << endl;
    m_output_stream << "Variance (σ^2) = " << m_gene_calculator.getVariance() << endl;
    m_output_stream << "Standard Deviation (σ) = " << m_gene_calculator.getStandardDeviation() << endl;
    m_output_stream << "P(A) = " << m_gene_calculator.getProbabilityA() << endl;
    m_output_stream << "P(C) = " << m_gene_calculator.getProbabilityC() << endl;
    m_output_stream << "P(G) = " << m_gene_calculator.getProbabilityG() << endl;
    m_output_stream << "P(T) = " << m_gene_calculator.getProbabilityT() << endl;
    m_output_stream << "P(AA) = " << m_gene_calculator.getProbabilityAA() << endl;
    m_output_stream << "P(AC) = " << m_gene_calculator.getProbabilityAC() << endl;
    m_output_stream << "P(AG) = " << m_gene_calculator.getProbabilityAG() << endl;
    m_output_stream << "P(AT) = " << m_gene_calculator.getProbabilityAT() << endl;
    m_output_stream << "P(CA) = " << m_gene_calculator.getProbabilityCA() << endl;
    m_output_stream << "P(CC) = " << m_gene_calculator.getProbabilityCC() << endl;
    m_output_stream << "P(CG) = " << m_gene_calculator.getProbabilityCG() << endl;
    m_output_stream << "P(CT) = " << m_gene_calculator.getProbabilityCT() << endl;
    m_output_stream << "P(GA) = " << m_gene_calculator.getProbabilityGA() << endl;
    m_output_stream << "P(GC) = " << m_gene_calculator.getProbabilityGC() << endl;
    m_output_stream << "P(GG) = " << m_gene_calculator.getProbabilityGG() << endl;
    m_output_stream << "P(GT) = " << m_gene_calculator.getProbabilityGT() << endl;
    m_output_stream << "P(TA) = " << m_gene_calculator.getProbabilityTA() << endl;
    m_output_stream << "P(TC) = " << m_gene_calculator.getProbabilityTC() << endl;
    m_output_stream << "P(TG) = " << m_gene_calculator.getProbabilityTG() << endl;
    m_output_stream << "P(TT) = " << m_gene_calculator.getProbabilityTT() << endl;

    m_output_stream << "\nNow, we will generate 1000 DNA strings that follow a Gaussian distribution with the same mean and variance as calculated above." << endl;

    // The genePool class is used to generate DNA strings that follow the Gaussian distribution of the read DNA-read file.
    // The for loop is simply used to generate 1000 of these strings and output them to a results file.
    for (int i = 0; i < 1000; ++i)
    {
      m_gene_pool.setDNA(m_gene_calculator.getMean(), m_gene_calculator.getStandardDeviation(), m_gene_calculator.getProbabilityA(), m_gene_calculator.getProbabilityC(), m_gene_calculator.getProbabilityG());
      // Since it is possible for the box-muller transform to give us a number less than 1 (leading to a 0 length of a DNA-read),
      // we must use this if statement to catch null strings. By decrementing i, we will ensure that we generate exactly
      // 1000 strings.
      if (m_gene_pool.getDNA() == "")
      {
        --i;
        continue;
      }
      // Each string is output here, with its index number to show that it is the nth string.
      m_output_stream << (i + 1) << ": " << m_gene_pool.getDNA() << endl;
    }

    m_output_stream << "\n";

    m_output_stream.close();
    cout << "Done! Results were saved to VidalArroyo.txt.\n";
  }
}

// The AppendOut method is nearly indistinguisgable from the Out method, except for the fact that this method appends
// to the 'VidalArroyo.txt' results file and it does not put the name and student ID in the file again (to avoid
// redundancy). As a result of this high similarity between the two methods, no extra comments will be added here.
void Logger::AppendOut(string f)
{
  m_input_stream.clear();
  m_input_stream.open(f);
  if (m_input_stream.fail())
  {
    cout << "File does not exist!" << endl;
  }
  else
  {
    while(getline(m_input_stream, m_single_line))
    {
      for (int i = 0; i < m_single_line.size(); ++i)
      {
        if (toupper(m_single_line.at(i)) == 'A')
        {
          m_single_line_cleaned += 'A';
        }
        else if (toupper(m_single_line.at(i)) == 'C')
        {
          m_single_line_cleaned += 'C';
        }
        else if (toupper(m_single_line.at(i)) == 'G')
        {
          m_single_line_cleaned += 'G';
        }
        else if (toupper(m_single_line.at(i)) == 'T')
        {
          m_single_line_cleaned += 'T';
        }
        else
        {
          continue;
        }
      }
      m_entire_file_read += (m_single_line_cleaned + 'x');
      m_single_line_cleaned = "";
    }

    m_input_stream.close();

    m_gene_calculator.ComputeStats(m_entire_file_read);

    m_output_stream.open("VidalArroyo.txt", std::ios_base::app);

    m_output_stream << "Here are the summary statistics of " << f << "\n" << endl;
    m_output_stream << "Sum (Σ) = " << m_gene_calculator.getSum() << endl;
    m_output_stream << "Mean (μ) = " << m_gene_calculator.getMean() << endl;
    m_output_stream << "Variance (σ^2) = " << m_gene_calculator.getVariance() << endl;
    m_output_stream << "Standard Deviation (σ) = " << m_gene_calculator.getStandardDeviation() << endl;
    m_output_stream << "P(A) = " << m_gene_calculator.getProbabilityA() << endl;
    m_output_stream << "P(C) = " << m_gene_calculator.getProbabilityC() << endl;
    m_output_stream << "P(G) = " << m_gene_calculator.getProbabilityG() << endl;
    m_output_stream << "P(T) = " << m_gene_calculator.getProbabilityT() << endl;
    m_output_stream << "P(AA) = " << m_gene_calculator.getProbabilityAA() << endl;
    m_output_stream << "P(AC) = " << m_gene_calculator.getProbabilityAC() << endl;
    m_output_stream << "P(AG) = " << m_gene_calculator.getProbabilityAG() << endl;
    m_output_stream << "P(AT) = " << m_gene_calculator.getProbabilityAT() << endl;
    m_output_stream << "P(CA) = " << m_gene_calculator.getProbabilityCA() << endl;
    m_output_stream << "P(CC) = " << m_gene_calculator.getProbabilityCC() << endl;
    m_output_stream << "P(CG) = " << m_gene_calculator.getProbabilityCG() << endl;
    m_output_stream << "P(CT) = " << m_gene_calculator.getProbabilityCT() << endl;
    m_output_stream << "P(GA) = " << m_gene_calculator.getProbabilityGA() << endl;
    m_output_stream << "P(GC) = " << m_gene_calculator.getProbabilityGC() << endl;
    m_output_stream << "P(GG) = " << m_gene_calculator.getProbabilityGG() << endl;
    m_output_stream << "P(GT) = " << m_gene_calculator.getProbabilityGT() << endl;
    m_output_stream << "P(TA) = " << m_gene_calculator.getProbabilityTA() << endl;
    m_output_stream << "P(TC) = " << m_gene_calculator.getProbabilityTC() << endl;
    m_output_stream << "P(TG) = " << m_gene_calculator.getProbabilityTG() << endl;
    m_output_stream << "P(TT) = " << m_gene_calculator.getProbabilityTT() << endl;

    m_output_stream << "\nNow, we will generate 1000 DNA strings that follow a Gaussian distribution with the same mean and variance as calculated above." << endl;

    for (int i = 0; i < 1000; ++i)
    {
      m_gene_pool.setDNA(m_gene_calculator.getMean(), m_gene_calculator.getStandardDeviation(), m_gene_calculator.getProbabilityA(), m_gene_calculator.getProbabilityC(), m_gene_calculator.getProbabilityG());
      if (m_gene_pool.getDNA() == "")
      {
        --i;
        continue;
      }
      m_output_stream << (i + 1) << ": " << m_gene_pool.getDNA() << endl;
    }

    m_output_stream << "\n";

    m_output_stream.close();
    cout << "Done! Results were appended to VidalArroyo.txt.\n";
  }
}
