// geneCalculator.h
// A declared calcultor program used to compute statistical features for a DNA reader program.
// @author Vidal M. Arroyo
// @author arroy118@mail.chapman.edu
// @version 1.0

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;

// Given a string of DNA reads separated by characters, this file can compute statistical features for a given DNA-read file.
class GeneCalculator
{
  public:
    GeneCalculator();
    ~GeneCalculator();
    int getSum();
    double getMean();
    double getVariance();
    double getStandardDeviation();
    double getProbabilityA();
    double getProbabilityC();
    double getProbabilityG();
    double getProbabilityT();
    double getProbabilityAA();
    double getProbabilityAC();
    double getProbabilityAG();
    double getProbabilityAT();
    double getProbabilityCA();
    double getProbabilityCC();
    double getProbabilityCG();
    double getProbabilityCT();
    double getProbabilityGA();
    double getProbabilityGC();
    double getProbabilityGG();
    double getProbabilityGT();
    double getProbabilityTA();
    double getProbabilityTC();
    double getProbabilityTG();
    double getProbabilityTT();
    // Given the interdependence of statistical parameters on one-another (I.E. the variance is needed to compute the
    // standard deviation, etc.), I created this method to act as a universal mutator of all statistical parameters at
    // the same time.
    void ComputeStats(string entire_file_read);
  private:
    // The count variables are used to keep track of how many A's, C's, G's, and T's were present in the DNA read file.
    // This will come into use later when we need to compute the probability of each base pair.
    int m_a_count;
    int m_c_count;
    int m_g_count;
    int m_t_count;
    // The number of reads variable stores the number of DNA 'lines' that were in the file, given that each line represents
    // a single gene entity. The term read is used commonly in molecular biology applications, which is why I used the
    // same term here.
    int m_number_of_reads;
    int m_read_length;
    int m_sum;
    double m_mean;
    // To simplify computations, I created a) temporary and b) sum variance variables to be able to a) determine the variance
    // of a single DNA read's length from the mean, and then b) add this to the previous variances of other reads. Using
    // this technique, I was able to greatly simplify the computation of the variance of the sample.
    double m_temporary_variance;
    double m_sum_variance;
    double m_variance;
    double m_standard_deviation;
    double m_probability_a;
    double m_probability_c;
    double m_probability_g;
    double m_probability_t;
    // The bigram probabilities were computed with a different sum from the A, G, C, and T probabilities. I detail what
    // this denominator is (and why it was used) in the geneCalculator.cpp file.
    double m_probability_aa;
    double m_probability_ac;
    double m_probability_ag;
    double m_probability_at;
    double m_probability_ca;
    double m_probability_cc;
    double m_probability_cg;
    double m_probability_ct;
    double m_probability_ga;
    double m_probability_gc;
    double m_probability_gg;
    double m_probability_gt;
    double m_probability_ta;
    double m_probability_tc;
    double m_probability_tg;
    double m_probability_tt;
};
