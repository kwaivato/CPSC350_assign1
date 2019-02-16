// genePool.cpp
// An implemented program used to generate DNA strings from a gaussian DNA distribution for a DNA reader program.
// @author Vidal M. Arroyo
// @author arroy118@mail.chapman.edu
// @version 1.0

#include "genePool.h"

using namespace std;

// The default values set in this constructor are set to permit optimal computation of DNA strings.
GenePool::GenePool()
{
  BoxMuller m_box_muller_transform;
  m_dna = "";
  m_random_num = 0;
  m_base_pair = ' ';
}

GenePool::~GenePool()
{

}

void GenePool::setDNA(double mean, double standard_deviation, double probability_a, double probability_c, double probability_g)
{
  m_dna = "";
  // The box-muller transform is used here to give us a randomized length, D,  of a DNA string that follows the mean and
  // variance of a previous DNA-read file. See the box-muller transform class for technical details.
  m_box_muller_transform.setC(mean);
  m_box_muller_transform.setD(mean,standard_deviation);
  // Here, we cast D to be an int because the length of a DNA string can only be an integer. Due to the box-muller transform,
  // D is a normally distributed double.
  for (int i = 0; i < static_cast<int>(m_box_muller_transform.getD()); ++i)
  {
    // This random integer, m_rand_num, will be used to give us a probability that we can use to determine the base pair
    // composition of a given strand of DNA.
    m_random_num = (RAND_MAX - rand())/static_cast<double>(RAND_MAX);
    // Due to the nature of rand(), m_random_num will be a random probability between 0 and 1. If this probability is less
    // than the probability of A, then A is chosen as the base pair. The higher the probability of A, the more likely this
    // is to happen. We apply similar methods to both C and G, with the caveat that we add the previous probabilities
    // since this is the only way another base pair can be chosen. If we don't do this, then we will get more T's than
    // we should.
    if (m_random_num <= probability_a){
      m_base_pair = 'A';
      m_dna += m_base_pair;
    }else if (m_random_num <= (probability_a + probability_c)){
      m_base_pair = 'C';
      m_dna += m_base_pair;
    }else if (m_random_num <= (probability_a + probability_c + probability_g)){
      m_base_pair = 'G';
      m_dna += m_base_pair;
    // The else (and not else if) statement is used because the probability of A, C, G, and T must add up to 1. Additionally,
    // if the probability of A, C, and G are very low, then their sum will be low and this path will be chosen more often -
    // which is exactly what we want.
    }else{
      m_base_pair = 'T';
      m_dna += m_base_pair;
    }
  }
}

string GenePool::getDNA()
{
  return m_dna;
}
