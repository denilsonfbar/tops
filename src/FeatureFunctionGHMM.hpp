/*
 *       FeatureFunctionGHMM.hpp
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

#ifndef FEATURE_FUNCTION_GHMM_HPP
#define FEATURE_FUNCTION_GHMM_HPP

#include "GeneralizedHiddenMarkovModel.hpp"

namespace tops {

  class DLLEXPORT FeatureFunctionGHMMTransitions : public FeatureFunction {
    private:
      GeneralizedHiddenMarkovModelPtr _ghmm_model;
    public:
      FeatureFunctionGHMMTransitions(GeneralizedHiddenMarkovModelPtr ghmm_model) : _ghmm_model(ghmm_model) {}
      virtual ~FeatureFunctionGHMMTransitions() {}

      virtual double ff(int y_tm1, int t, int b_t, int e_t, int y_t, const Sequence& x) const {
        double value;
          if (t == 0){
            #if VERBOSE_DETAILS
              cout << "init :\t"; print_parameters( y_tm1,t,b_t,e_t,y_t,x);
            #endif
            value = _ghmm_model->getInitialProbabilities()->log_probability_of(y_t);
          }
          else{
            #if VERBOSE_DETAILS
              cout << "trans:\t"; print_parameters( y_tm1,t,b_t,e_t,y_t,x);
            #endif
            value = _ghmm_model->getAllStates()[y_tm1]->transition()->log_probability_of(y_t);
          }

        #if VERBOSE_DETAILS
          std::cout << "\tlog: " << value << "\treal: " << exp(value) << std::endl;
        #endif

        return value;
      }
  };
  typedef boost::shared_ptr<FeatureFunctionGHMMTransitions> FeatureFunctionGHMMTransitionsPtr;


  class DLLEXPORT FeatureFunctionGHMMObservations : public FeatureFunction {
    private:
      GeneralizedHiddenMarkovModelPtr _ghmm_model;
    public:
      FeatureFunctionGHMMObservations(GeneralizedHiddenMarkovModelPtr ghmm_model) : _ghmm_model(ghmm_model) {}
      virtual ~FeatureFunctionGHMMObservations() {}

      virtual double ff(int y_tm1, int t, int b_t, int e_t, int y_t, const Sequence& x) const {

        // TODO: para executar initialize_prefix_sum_arrays(x);
        if (t == 0 && b_t == 0 && y_t == 0){  //  && e_t == 0
          Matrix v;
          Sequence states;
          double score_ghmm_viterbi = _ghmm_model->_viterbi(x, states, v);  // Inefficient GHMM Viterbi algorithm

          #if VERBOSE_DETAILS
            std::cout << "--> Executed initialize_prefix_sum_arrays(x)" << std::endl;
          #endif
        }

        double value = _ghmm_model->getAllStates()[y_t]->observation()->prefix_sum_array_compute(b_t,e_t);

        #if VERBOSE_DETAILS
          std::cout << "obser:\t"; print_parameters( y_tm1,t,b_t,e_t,y_t,x);
          std::cout << "\tlog: " << value << "\treal: " << exp(value) << std::endl;
        #endif

        return value;
      }
  };
  typedef boost::shared_ptr<FeatureFunctionGHMMObservations> FeatureFunctionGHMMObservationsPtr;


  class DLLEXPORT FeatureFunctionGHMMDurations : public FeatureFunction {
    private:
      GeneralizedHiddenMarkovModelPtr _ghmm_model;
    public:
      FeatureFunctionGHMMDurations(GeneralizedHiddenMarkovModelPtr ghmm_model) : _ghmm_model(ghmm_model) {}
      virtual ~FeatureFunctionGHMMDurations() {}

      virtual double ff(int y_tm1, int t, int b_t, int e_t, int y_t, const Sequence& x) const {

        int d = e_t - b_t + 1;
        double value = _ghmm_model->getAllStates()[y_t]->duration_probability(d);

        #if VERBOSE_DETAILS
          std::cout << "durat:\t"; print_parameters( y_tm1,t,b_t,e_t,y_t,x);
          std::cout << "\tlog: " << value << "\treal: " << exp(value) << std::endl;
        #endif

        return value;
      }
  };
  typedef boost::shared_ptr<FeatureFunctionGHMMDurations> FeatureFunctionGHMMDurationsPtr;

}

#endif
