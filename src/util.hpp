#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include <fstream>
#include <boost/regex.hpp>
#include <vector>
#include <string>
#include <cmath>
#include "Sequence.hpp"
#include "SequenceEntry.hpp"
#include <boost/numeric/ublas/matrix.hpp>
namespace tops {

  typedef std::vector <double> DoubleVector;
  typedef std::vector <int> IntVector;
  typedef boost::numeric::ublas::matrix<double> Matrix;
  typedef std::vector <std::string> StringVector;
  typedef std::map <std::string, std::string> StringMap;

  //! Remove spaces from the ends of the string 
  void trim_spaces (std::string & s);

  //! Split the string by using a separator 
  void split_regex (const std::string & s, std::vector <std::string> & result, const boost::regex & re);

  //! Calculates the value of log(exp(log_a) + exp(log_b)) 
  double log_sum( double log_a, double log_b);

  //! Divides the a by b
  double safe_division(double a, double b);

  //! Returns true if a is close to b with a given tolerance
  bool close(double a, double b, double tolerance) ;

  ///! Calculates the value of D mod d
  int mod(int D, int d);

  //! Reads a list of sequences from a file
  void readSequencesFromFile(SequenceEntryList & s, AlphabetPtr alphabet, std::string file_name) ;

  //! Reads a list of sequences from a file
  void readSequencesFromFile(SequenceList & s, AlphabetPtr alphabet, std::string file_name) ;


  /*   Sheather and Jones bandwidth */
  double sj_bandwidth(const DoubleVector &data);
  double kernel_density_estimation(double x, double bw, const DoubleVector &data);
  double kernel_density_estimation_gaussian(double x, double bw, const DoubleVector &data);

  /* Epanechnikov kernel */
  double epanechnikov(double x, double h);
  
// code from R-1.7.0/src/appl/bandwidths.c  
#define abs9(a) (a > 0 ? a:-a)
  void band_den_bin(int n, int nb, double *d, const DoubleVector &x,  DoubleVector &cnt);
  void band_phi6_bin(int n, int nb, double d, DoubleVector &x, double h, double *u);
  void band_phi4_bin(int n, int nb, double d, DoubleVector x, double h, double *u);  

  double mean(const DoubleVector &data);
  
  double var(const DoubleVector &data);  
  /* quantile */
  double quantile (DoubleVector data, double q);

    /* interquantile */
  double iqr (const DoubleVector &data);
  
  


}
#endif
