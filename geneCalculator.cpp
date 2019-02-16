// geneCalculator.cpp
// An implemented calcultor program used to compute statistical features for a DNA reader program.
// @author Vidal M. Arroyo
// @author arroy118@mail.chapman.edu
// @version 1.0

#include "geneCalculator.h"

using namespace std;

// In order for this class to properly analyze a DNA-read distribution, separate DNA reads must be separated by x's.

GeneCalculator::GeneCalculator()
{
  // To make things easy, all statistical features (and variables used to compute them) are initiated to 0.
  m_a_count = 0;
  m_c_count = 0;
  m_g_count = 0;
  m_t_count = 0;
  m_number_of_reads = 0;
  m_read_length = 0;
  m_sum = 0;
  m_mean = 0.0;
  m_temporary_variance = 0.0;
  m_sum_variance = 0.0;
  m_variance = 0.0;
  m_standard_deviation = 0.0;
  m_probability_a = 0.0;
  m_probability_c = 0.0;
  m_probability_g = 0.0;
  m_probability_t = 0.0;
  m_probability_aa = 0.0;
  m_probability_ac = 0.0;
  m_probability_ag = 0.0;
  m_probability_at = 0.0;
  m_probability_ca = 0.0;
  m_probability_cc = 0.0;
  m_probability_cg = 0.0;
  m_probability_ct = 0.0;
  m_probability_ga = 0.0;
  m_probability_gc = 0.0;
  m_probability_gg = 0.0;
  m_probability_gt = 0.0;
  m_probability_ta = 0.0;
  m_probability_tc = 0.0;
  m_probability_tg = 0.0;
  m_probability_tt = 0.0;
}

GeneCalculator::~GeneCalculator()
{

}

// As previously detailed, the ComputeStats method acts as a universal 'mutator' for all statistical features of a DNA-read
// file. This was done for both computational and mathematical simplicity, since some statistical measurements are used
// to compute others (I.E. the mean is used to compute the variance, which is used to compute the standard deviation, etc.)
void GeneCalculator::ComputeStats(string entire_file_read)
{
  // This calculator starts by finding out how many total and individual base pairs (A, C, G, and T) are in the DNA-read
  // file. The Logger class creates a curated 'entire_file_read' string that works very well with this class.
  for (int i = 0; i < entire_file_read.size(); ++i)
  {
    if (entire_file_read.at(i) == 'A')
    {
      // Throughout the calculator, pre-increments are used instead of post-increments for the slight benefits in performance.
      ++m_a_count;
      ++m_sum;
    }
    else if (entire_file_read.at(i) == 'C')
    {
      ++m_c_count;
      ++m_sum;
    }
    else if (entire_file_read.at(i) == 'G')
    {
      ++m_g_count;
      ++m_sum;
    }
    else if (entire_file_read.at(i) == 'T')
    {
      ++m_t_count;
      ++m_sum;
    }
    else
    {
      continue;
    }
  }

  // Here, a for loop is used to find out how many individual DNA 'reads' are in a DNA-read file. This method works under
  // the assumption that DNA reads are separatd by x's, along with characters that are not base pairs. The logger class
  // does this, which is why this calculator works best (and possibly only) with the logger class.
  for(int i = 0; i < entire_file_read.size(); ++i)
  {
    // If a character is an x, we don't care about it. So we continue.
    if(entire_file_read.at(i) == 'x')
    {
      continue;
    }

    // The following if/else-if statements are used to check if if we have reached the end of a read or the end of the file.
    // This first if statement checks if we are at the end of the file. If we are, then this was our last read. If not, then
    // we proceed. This statement MUST be checked first because, if not, then we will encounter an error at the end of every
    // file.
    if(i == (entire_file_read.size() - 1))
    {
      m_number_of_reads++;
    }
    // This is used to check if we reach the end of a read. As previously stated, this assumes that every read is
    // separed by x's. If this is the case, then we have another read, and the number of reads is increased.
    else if (entire_file_read.at(i+1) == 'x')
    {
      m_number_of_reads++;
    }
    // If we reach this statement, it means that we don't have an x, don't have an x ahead of us, and are not at the end
    // of the file. This means we are in the middle of a read. So we continue!
    else
    {
      continue;
    }

  }

  // The mean length is computed here, dividing the sum of lengths by the number of reads.
  m_mean = static_cast<double>(m_sum)/static_cast<double>(m_number_of_reads);


  // Now that we have the mean we can compute the variance. Since this algorithm is very similar to that of the
  // number of reads algorithm above, I will only comment where deviations exist.
  for(int i = 0; i < entire_file_read.size(); ++i)
  {
    if(entire_file_read.at(i) == 'x')
    {
      continue;
    }

    // If the entire_file_read gets past the previous if statement, then it means that it has an A, C, G, or T. So, we want
    // to increment the read length. The read length will be used to compute the variance for a given read.
    ++m_read_length;

    if(i == (entire_file_read.size() - 1))
    {
      // The variance is computed as the sum of square differences between each length and the mean length. So, when designing
      // this algorithm, I decided it would be easiest to compute the variance for each read (m_temporary_variance) and then
      // add them together (m_sum_variance). Then, we could divide m_sum_variance by n to give us the real variance.
      m_temporary_variance = pow((m_read_length - m_mean),2);
      m_sum_variance += m_temporary_variance;
      m_read_length = 0;
    }
    else if (entire_file_read.at(i+1) == 'x')
    {
      m_temporary_variance = pow((m_read_length - m_mean),2);
      m_sum_variance += m_temporary_variance;
      m_read_length = 0;
    }
    else
    {
      continue;
    }

  }

  // As detailed above, the variance can be computed by dividing the sum of variance by the number of reads.
  m_variance = m_sum_variance/static_cast<double>(m_number_of_reads);


  // The standard deviation, which is simply the square root of the variance, is computed here.
  m_standard_deviation = sqrt(m_variance);

  // The probability of each nucleotide is computed here by diving the count of each nucleotide by the sum of all
  // nucleotides.
  m_probability_a = static_cast<double>(m_a_count)/static_cast<double>(m_sum);
  m_probability_c = static_cast<double>(m_c_count)/static_cast<double>(m_sum);
  m_probability_g = static_cast<double>(m_g_count)/static_cast<double>(m_sum);
  m_probability_t = static_cast<double>(m_t_count)/static_cast<double>(m_sum);

  // In this entire for-loop, all of the bi-gram probabilities are computed. I will only detail the first couple lines
  // to show how this works.
  for(int i = 0; i < entire_file_read.size(); ++i)
  {
    // This first checks whether the algorithm is at the end of the file. If so, then it breaks.
    if(i == (entire_file_read.size() - 1))
    {
      break;
    }
    // This checks if the first nucleotide is an 'A'. This is similarly done for C, G, and T further below.
    else if (entire_file_read.at(i) == 'A')
    {
      // This checks whether we are at the last nucleotide, which (if so) means that we cannot compute a bigram probability
      // since there is only one nucleotide left. So, it breaks.
      if(i == (entire_file_read.size() - 2))
      {
        break;
      }
      // This checks if the second nucleotide is an 'A'. This is similarly done for C, G, and T further below.
      else if (entire_file_read.at(i+1) == 'A')
      {
        // After doing a proof by induction, I was able to show that the total number of bigrams in a DNA read file
        // is equal to the sum of all nucleotides divided by the number of reads. If you would like me to prove this to you,
        // I would be happy to in office hours.
        m_probability_aa += 1/(static_cast<double>(m_sum - m_number_of_reads));
      }
      else if (entire_file_read.at(i+1) == 'C')
      {
        m_probability_ac += 1/(static_cast<double>(m_sum - m_number_of_reads));
      }
      else if (entire_file_read.at(i+1) == 'G')
      {
        m_probability_ag += 1/(static_cast<double>(m_sum - m_number_of_reads));
      }
      else if (entire_file_read.at(i+1) == 'T')
      {
        m_probability_at += 1/(static_cast<double>(m_sum - m_number_of_reads));
      }
      else
      {
        continue;
      }
    }
    else if (entire_file_read.at(i) == 'C')
    {
      if(i == (entire_file_read.size() - 2))
      {
        break;
      }
      else if (entire_file_read.at(i+1) == 'A')
      {
        m_probability_ca += 1/(static_cast<double>(m_sum - m_number_of_reads));
      }
      else if (entire_file_read.at(i+1) == 'C')
      {
        m_probability_cc += 1/(static_cast<double>(m_sum - m_number_of_reads));
      }
      else if (entire_file_read.at(i+1) == 'G')
      {
        m_probability_cg += 1/(static_cast<double>(m_sum - m_number_of_reads));
      }
      else if (entire_file_read.at(i+1) == 'T')
      {
        m_probability_ct += 1/(static_cast<double>(m_sum - m_number_of_reads));
      }
      else
      {
        continue;
      }
    }
    else if (entire_file_read.at(i) == 'G')
    {
      if(i == (entire_file_read.size() - 2))
      {
        break;
      }
      else if (entire_file_read.at(i+1) == 'A')
      {
        m_probability_ga += 1/(static_cast<double>(m_sum - m_number_of_reads));
      }
      else if (entire_file_read.at(i+1) == 'C')
      {
        m_probability_gc += 1/(static_cast<double>(m_sum - m_number_of_reads));
      }
      else if (entire_file_read.at(i+1) == 'G')
      {
        m_probability_gg += 1/(static_cast<double>(m_sum - m_number_of_reads));
      }
      else if (entire_file_read.at(i+1) == 'T')
      {
        m_probability_gt += 1/(static_cast<double>(m_sum - m_number_of_reads));
      }
      else
      {
        continue;
      }
    }
    else if (entire_file_read.at(i) == 'T')
    {
      if(i == (entire_file_read.size() - 2))
      {
        break;
      }
      else if (entire_file_read.at(i+1) == 'A')
      {
        m_probability_ta += 1/(static_cast<double>(m_sum - m_number_of_reads));
      }
      else if (entire_file_read.at(i+1) == 'C')
      {
        m_probability_tc += 1/(static_cast<double>(m_sum - m_number_of_reads));
      }
      else if (entire_file_read.at(i+1) == 'G')
      {
        m_probability_tg += 1/(static_cast<double>(m_sum - m_number_of_reads));
      }
      else if (entire_file_read.at(i+1) == 'T')
      {
        m_probability_tt += 1/(static_cast<double>(m_sum - m_number_of_reads));
      }
      else
      {
        continue;
      }
    }
    else
    {
      continue;
    }
  }
}

int GeneCalculator::getSum()
{
  return m_sum;
}

double GeneCalculator::getMean()
{
  return m_mean;
}

double GeneCalculator::getVariance()
{
  return m_variance;
}

double GeneCalculator::getStandardDeviation()
{
  return m_standard_deviation;
}

double GeneCalculator::getProbabilityA()
{
  return m_probability_a;
}

double GeneCalculator::getProbabilityC()
{
  return m_probability_c;
}

double GeneCalculator::getProbabilityG()
{
  return m_probability_g;
}

double GeneCalculator::getProbabilityT()
{
  return m_probability_t;
}

double GeneCalculator::getProbabilityAA()
{
  return m_probability_aa;
}

double GeneCalculator::getProbabilityAC()
{
  return m_probability_ac;
}

double GeneCalculator::getProbabilityAG()
{
  return m_probability_ag;
}

double GeneCalculator::getProbabilityAT()
{
  return m_probability_at;
}

double GeneCalculator::getProbabilityCA()
{
  return m_probability_ca;
}

double GeneCalculator::getProbabilityCC()
{
  return m_probability_cc;
}

double GeneCalculator::getProbabilityCG()
{
  return m_probability_cg;
}

double GeneCalculator::getProbabilityCT()
{
  return m_probability_ct;
}

double GeneCalculator::getProbabilityGA()
{
  return m_probability_ga;
}

double GeneCalculator::getProbabilityGC()
{
  return m_probability_gc;
}

double GeneCalculator::getProbabilityGG()
{
  return m_probability_gg;
}

double GeneCalculator::getProbabilityGT()
{
  return m_probability_gt;
}

double GeneCalculator::getProbabilityTA()
{
  return m_probability_ta;
}

double GeneCalculator::getProbabilityTC()
{
  return m_probability_tc;
}

double GeneCalculator::getProbabilityTG()
{
  return m_probability_tg;
}

double GeneCalculator::getProbabilityTT()
{
  return m_probability_tt;
}
