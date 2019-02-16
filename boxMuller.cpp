// boxMuller.cpp
// An implemented program used to generate random lengths of DNA strings from a gaussian DNA distribution for a DNA reader program.
// @author Vidal M. Arroyo
// @author arroy118@mail.chapman.edu
// @version 1.0

#include "boxMuller.h"

using namespace std;

BoxMuller::BoxMuller()
{

}

BoxMuller::~BoxMuller()
{

}

double BoxMuller::getD()
{
  return m_d;
}

void BoxMuller::setC(double mean)
{
  // The goal of this line was to create a truly random set of DNA read lengths. Since seeding rand() with srand(int x)
  // for an x that I choose beforehand will always seed the lengths to be the same, I wanted to modify this and seed the
  // lengths based on the unique characters of a given DNA-read file. This was my rational for seeding it with an
  // int version of the mean + rand(), the latter which will always generate a random number. Linearly transforming
  // rand() by adding the mean to it will allow each file to generate a unique set of strings for its statistical
  // character.
  srand(rand() + static_cast<int>(mean));
  // This is used to set A as a randomly distributed number between 0 and 1.
  m_a = (RAND_MAX - rand( ))/ static_cast<double>(RAND_MAX);
  // This is used to set B as a randomly distributed number between 0 and 1.
  m_b = (RAND_MAX - rand( ))/ static_cast<double>(RAND_MAX);
  // This is used to compute a random variable C from A and B using the Box-Muller transform.
  m_c = (sqrt(-2*log(m_a)))*cos(2*M_PI*m_b);
}

void BoxMuller::setD(double mean, double standard_deviation)
{
  // Given that C was computed as a normal Gaussian with mean 0 and variance 1, we can conver C to a normal random variable
  // D with the mean and variance of the DNA-read lengths. In this transformation, we use the standard deviation - which is
  // the square root of the variance.
  m_d = (standard_deviation*m_c) + mean;
}
