// boxMuller.h
// A declared program used to generate random lengths of DNA strings from a gaussian DNA distribution for a DNA reader program.
// @author Vidal M. Arroyo
// @author arroy118@mail.chapman.edu
// @version 1.0

#include <cmath>
#include <cstdlib>
#include <string>

using namespace std;

// This class mathematically implements the box-muller transform given some mean and standard deviation.
class BoxMuller{
  public:
    BoxMuller();
    ~BoxMuller();
    double getD();
    // The setC mutator will set the value for the random variable C, which is a constant used in the Box-Muller transform.
    void setC(double mean);
    // The setD mutator will set the value for D, which is a normal random variable generated through the Box-Muller
    // transform.
    void setD(double mean, double standard_deviation);
  private:
    // A is random variable that is uniformly distributed between [0,1) that is used in the Box-Muller transform.
    double m_a;
    // B is random variable that is uniformly distributed between [0,1) that is used in the Box-Muller transform.
    double m_b;
    // C is a random variable computed from A and B that is used to compute a normal random variable D.
    double m_c;
    // In this case, D is a double because I want this class to be applicable to other program that use the Box-Muller
    // transform.
    double m_d;
};
