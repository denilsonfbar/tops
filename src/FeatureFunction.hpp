/*
 *       FeatureFunction.hpp
 *
 *       Copyright 2011 Andre Yoshiaki Kashiwabara <akashiwabara@usp.br>
 *                      Ígor Bonadio <ibonadio@ime.usp.br>
 *                      Vitor Onuchic <vitoronuchic@gmail.com>
 *                      Alan Mitchell Durham <aland@usp.br>
 *                 2022 Denilson Fagundes Barbosa <denilsonfbar@gmail.com>
 *
 *       This program is free software; you can redistribute it and/or modify
 *       it under the terms of the GNU  General Public License as published by
 *       the Free Software Foundation; either version 3 of the License, or
 *       (at your option) any later version.
 *
 *       This program is distributed in the hope that it will be useful,
 *       but WITHOUT ANY WARRANTY; without even the implied warranty of
 *       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *       GNU General Public License for more details.
 *
 *       You should have received a copy of the GNU General Public License
 *       along with this program; if not, write to the Free Software
 *       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *       MA 02110-1301, USA.
 */

#ifndef FEATURE_FUNCTION_HPP
#define FEATURE_FUNCTION_HPP

#define VERBOSE 0
#define VERBOSE_DETAILS 0

namespace tops {

  // Interface for Semi-Markov CRF features functions  
  class DLLEXPORT FeatureFunction {
    public:
      FeatureFunction() {}
      virtual ~FeatureFunction() {}

      virtual double ff(int y_tm1, int t, int b_t, int e_t, int y_t, const Sequence& x) const = 0;

      void print_parameters(int y_tm1, int t, int b_t, int e_t, int y_t, const Sequence& x) const {
        std::cout << "t: " << t << "\tb: " << b_t << "\te: " << e_t << "\ty: " << y_t << "\ty-1: " << y_tm1 << "\tx[t]: " << x[t];
      }
  };
  typedef boost::shared_ptr<FeatureFunction> FeatureFunctionPtr;
  
}

#endif
