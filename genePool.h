// genePool.h
// A declared program used to generate DNA strings from a gaussian DNA distribution for a DNA reader program.
// @author Vidal M. Arroyo
// @author arroy118@mail.chapman.edu
// @version 1.0

#include <string>

#include "boxMuller.h"

using namespace std;

// This class can be used to generate DNA strings from the mean and standard deviation of another file's DNA string
// distribution. Must use with the boxMuller class.
class GenePool{
  public:
    GenePool();
    ~GenePool();
    string getDNA();
    // This method will create a random DNA string that comes from the distrubutin of a DNA-read file with a given
    // mean and variance of the DNA read-length. It will generate the strings with the same base-pair frequencies as seen
    // in a DNA-read file. Since it is assumed that P(A) + P(C) + P(G) + P(T) = 1, this function does not take in
    // probability_t as a method for simplicity and computational efficiency.
    void setDNA(double mean, double standard_deviation, double probability_a, double probability_c, double probability_g);
  private:
      // The Box-Muller transform is used to generate the length of our DNA string so that it stays stochastic.
      BoxMuller m_box_muller_transform;
      string m_dna;
      // A random number will be used to determine the base pair composition of the DNA string that this class generates.
      double m_random_num;
      // This variable will be used to store base-pairs for our DNA string.
      char m_base_pair;
};
