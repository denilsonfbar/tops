/*
 *       SemiMarkovCRFModel.cpp
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

#include <boost/numeric/ublas/matrix.hpp>
#include "Alphabet.hpp"
#include "SemiMarkovCRFModel.hpp"
#include "ProbabilisticModelParameter.hpp"
#include "Symbol.hpp"
#include <iostream>
#include <cmath>
#include <sstream>
#include <vector>
#include <iterator>
#include <stdio.h>

namespace tops {

  ProbabilisticModelParameters SemiMarkovCRFModel::parameters() const {
    ProbabilisticModelParameters answer;
    /*
    int nstates = _states.size();
    answer.add("model_name", StringParameterValuePtr(new StringParameterValue(model_name().c_str())));
    answer.add("state_names", _state_names->getParameterValue());
    answer.add("observation_symbols", alphabet()->getParameterValue());
    std::map <std::string, double> trans;
    std::stringstream out;
    out << getStateName(0) << "|" << getStateName(0) ;
    trans[out.str()] =  exp(getState(0)->transitions()->log_probability_of(0)) ;
    for(int i = 0; i < nstates; i++)
      for(int j = 0; j < nstates; j++)
        if((i != 0) || (j != 0)) {
          std::stringstream out2;
          out2 << getStateName(j) << "|" << getStateName(i) ;
          trans[out2.str()] = exp(getState(i)->transitions()->log_probability_of(j));
        }
    answer.add("transitions", DoubleMapParameterValuePtr(new DoubleMapParameterValue(trans)));

    std::map <std::string, double> emission;
    std::stringstream out3;
    out3 <<   alphabet()->getSymbol(0)->name() << "|" << getStateName(0) ;
    emission[out3.str()] =  exp(getState(0)->emission()->log_probability_of(0));
    for(int i = 0; i < nstates; i++)
      for(int j = 0; j < (int)alphabet()->size(); j++)
        if((i != 0) || (j != 0)){
          std::stringstream out4;
          out4 <<alphabet()->getSymbol(j)->name() << "|" <<  getStateName(i) ;
          emission[out4.str()] =  exp(getState(i)->emission()->log_probability_of(j));
        }
    answer.add("emission_probabilities", DoubleMapParameterValuePtr(new DoubleMapParameterValue(emission)));
    double sum = 0;
    std::vector <double> probs;
    probs.resize(nstates);
    for(int i = 0; i < nstates; i++){
      probs[i] = exp(_initial_probability->log_probability_of(i));
      sum += probs[i];
    }
    std::map <std::string, double> initial;
    std::stringstream out5;
    out5 <<  getStateName(0);
    initial[out5.str()] =  probs[0]/sum;
    for(int i = 0; i < nstates; i++)
      for(int j = 0; j < (int)alphabet()->size(); j++)
        if((i != 0) || (j != 0)){
          std::stringstream out6;
          out6 << getStateName(i);
          initial[out6.str()] = probs[i]/sum;
        }
    answer.add("initial_probabilities", DoubleMapParameterValuePtr(new DoubleMapParameterValue(initial)));
    */
    return answer;
  }

  std::string SemiMarkovCRFModel::str() const {
    std::stringstream out;
    /*
    int nstates = _states.size();
    out << "model_name = \"" << model_name() << "\"" << std::endl;
    out << "state_names = (" ;
    out << "\"" << getStateName(0) << "\"";
    for(int i = 1; i < (int)getStateNames()->size(); i++)
      out << ",\"" << getStateName(i) << "\"";
    out << ")" << std::endl;

    out << "observation_symbols = (" ;
    out << "\"" << alphabet()->getSymbol(0)->name() << "\"";
    for(int i = 1; i < (int)alphabet()->size(); i++)
      out << ",\"" << alphabet()->getSymbol(i)->name() << "\"";
    out << ")" << std::endl;

    out << "transitions = (" ;
    out << "\"" << getStateName(0) << "\" | \"" << getStateName(0) << "\": " << exp(getState(0)->transitions()->log_probability_of(0)) ;
    for(int i = 0; i < nstates; i++)
      for(int j = 0; j < nstates; j++)
        if((i != 0) || (j != 0))
          out << ";\n \"" << getStateName(j) << "\" | \"" << getStateName(i) << "\": " << exp(getState(i)->transitions()->log_probability_of(j));
    out << ")" << std::endl;

    out << "emission_probabilities = (" ;
    out << "\"" <<  alphabet()->getSymbol(0)->name() << "\" | \"" << getStateName(0) << "\": " << exp(getState(0)->emission()->log_probability_of(0)) ;
    for(int i = 0; i < nstates; i++)
      for(int j = 0; j < (int)alphabet()->size(); j++)
        if((i != 0) || (j != 0))
          out << ";\n \"" << alphabet()->getSymbol(j)->name() << "\" | \"" <<  getStateName(i) << "\": " << exp(getState(i)->emission()->log_probability_of(j));
    out << ")" << std::endl;

    double sum = 0;
    std::vector <double> probs;
    probs.resize(nstates);
    for(int i = 0; i < nstates; i++){
      probs[i] = exp(_initial_probability->log_probability_of(i));
      sum += probs[i];
    }

    out << "initial_probabilities = (";
    out << "\"" << getStateName(0) << "\":  " << probs[0]/sum ;
    for(int i = 1; i < nstates; i++)
      out << ";\n \"" << getStateName(i) << "\": " << probs[i]/sum;
    out << ")" << std::endl;
    */
    return out.str();
  }

  std::string SemiMarkovCRFModel::getStateName(int state) const {
    // return getState(state)->getName()->name();
    return "";
  }

  AlphabetPtr SemiMarkovCRFModel::getStateNames() const {
    return _state_names;
  }


  void SemiMarkovCRFModel::initialize(const ProbabilisticModelParameters& parameters) {

    ProbabilisticModelParameterValuePtr state_names = parameters.getMandatoryParameterValue("state_names");
    ProbabilisticModelParameterValuePtr observation_symbols = parameters.getMandatoryParameterValue("observation_symbols");
    AlphabetPtr states = AlphabetPtr(new Alphabet());
    AlphabetPtr observations = AlphabetPtr(new Alphabet());
    states->initializeFromVector(state_names->getStringVector());
    observations->initializeFromVector(observation_symbols->getStringVector());
    setStates(states);
    setObservationSymbols(observations);

    // Creating a GHMM model for extract features:
    string model_name = "model/ghmm.model";
    ProbabilisticModelCreatorClient creator;
    ProbabilisticModelPtr model = creator.create(model_name);
    DecodableModelPtr dec_model = boost::dynamic_pointer_cast<DecodableModel> (model);
    GeneralizedHiddenMarkovModelPtr ghmm_model = boost::dynamic_pointer_cast<GeneralizedHiddenMarkovModel> (dec_model);
    // std::cout << ghmm_model->str();

    // Setup features functions:
    FeatureFunctionPtr feature_function_ghmm_transitions = FeatureFunctionGHMMTransitionsPtr(new FeatureFunctionGHMMTransitions(ghmm_model));
    _features_functions.push_back(feature_function_ghmm_transitions);
    _weights.push_back(1.0);

    FeatureFunctionPtr feature_function_ghmm_observations = FeatureFunctionGHMMObservationsPtr(new FeatureFunctionGHMMObservations(ghmm_model));
    _features_functions.push_back(feature_function_ghmm_observations);
    _weights.push_back(1.0);

    FeatureFunctionPtr feature_function_ghmm_durations = FeatureFunctionGHMMDurationsPtr(new FeatureFunctionGHMMDurations(ghmm_model));
    _features_functions.push_back(feature_function_ghmm_durations);
    _weights.push_back(1.0);
  }

   void SemiMarkovCRFModel::setStates(AlphabetPtr state_names) {
    _n_states = state_names->size();
    _state_names = state_names;
  }

  void SemiMarkovCRFModel::setObservationSymbols(AlphabetPtr obs) {
    tops::ProbabilisticModel::setAlphabet(obs);
  }

  double SemiMarkovCRFModel::sumFeatures(int y_tm1, int t, int b_t, int e_t, int y_t, const Sequence& x) const {

    double factor, weighted_factor, sum_weighted_factors = 0.0;
    int K = _features_functions.size();

    for (int k = 0; k < K; k++){
      factor = _features_functions[k]->ff(y_tm1, t, b_t, e_t, y_t, x);
      weighted_factor = _weights[k] * factor;
      sum_weighted_factors += weighted_factor;
    }

    #if VERBOSE_DETAILS
        std::cout << "Sum:\tlog: " << sum_weighted_factors << "\treal: " << exp(sum_weighted_factors) << std::endl;
    #endif

    return sum_weighted_factors;
  }

  double SemiMarkovCRFModel::viterbi(const Sequence& sequence, Sequence& path, Matrix& viterbi) const {

    typedef boost::numeric::ublas::matrix<int> MatrixInt;
    int nstates = _n_states;
    int size =  sequence.size();

    Matrix gamma(nstates, size, -HUGE);
    Matrix psi(nstates, size, 0 );
    IntMatrix psilen(nstates, size, 0);

    #if VERBOSE_DETAILS
      std::cout << "\nInicialization" << std::endl;
    #endif

    for (int k = 0; k < nstates; k++)
      gamma(k,0) = sumFeatures(-1, 0, 0, 0, k, sequence);

    for (int i = 1; i < size; i++){
      #if VERBOSE_DETAILS
        std::cout << "\n----Position: " << i << std::endl;
      #endif

      for (int k = 0; k < nstates; k++){
        #if VERBOSE_DETAILS
          std::cout << "\n---State: " << k << std::endl;
        #endif

        gamma(k,i) = -HUGE;

        for(int d = 1; d <= i; d++){
          #if VERBOSE_DETAILS
            std::cout << "\n--Distance: " << d << std::endl;
          #endif

          int b_t = i - d + 1;

          #if VERBOSE_DETAILS
            std::cout << "\n-Predess: " << 0 << std::endl;
          #endif

          double gmax = gamma(0, i-d) + sumFeatures(0, i, b_t, i, k, sequence);
          int pmax = 0;

          for(int p = 1; p < nstates; p++){
            #if VERBOSE_DETAILS
              std::cout << "\n-Predess: " << p << std::endl;
            #endif

            double g1 = gamma(p, i-d) + sumFeatures(p, i, b_t, i, k, sequence);
            if(gmax < g1){
              gmax = g1;
              pmax = p;
            }
          }

          if(gamma(k, i) < gmax){
            gamma(k, i) = gmax;
            psi(k, i) = pmax;
            psilen(k, i) = d;
          }
        }
      }
    }

    //backtracing
    path.resize(size);
    int L = size-1;
    int state = 0;
    double max = gamma(0, L);
    for(int k = 1; k < nstates; k++){
      if(max < gamma(k, L)){
        max = gamma(k, L);
        state = k;
      }
    }

    while(L > 0){
      int d = psilen(state, L);
      int p = psi(state, L);
      for(int i = 0; i < d; i++){
        path[L] = state;
        L--;
      }
      state = p;
    }

    #if VERBOSE
      std::cout << std::endl << "Inefficient SMCRF Viterbi algorithm:" << std::endl;
      std::cout << "Best path value: " << max << std::endl;
      std::cout << "Viterbi matrix: " << std::endl;
      print_matrix(gamma);
    #endif

    return max;
  }

  void SemiMarkovCRFModel::print_matrix(Matrix& m) const {
    std::cout << "Log values:" << std::endl;
    for(unsigned i = 0; i < m.size1(); ++i) {
      for(unsigned j = 0; j < m.size2(); ++j) {
            std::cout << m(i,j) << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << "Real values:" << std::endl;
    for(unsigned i = 0; i < m.size1(); ++i) {
      for(unsigned j = 0; j < m.size2(); ++j) {
            std::cout << exp(m(i,j)) << "\t";
        }
        std::cout << std::endl;
    }
  }


  Sequence& SemiMarkovCRFModel::chooseObservation(Sequence& h, int i, int state) const {
    // if((state >= 0) && (!getState(state)->isSilent()) )
    //   return getState(state)->emission()->chooseWithHistory(h,i,1);
    return h;
  }

  int SemiMarkovCRFModel::chooseState(int state) const {
    // return getState(state)->transitions()->choose();
    return -1;
  }

  int SemiMarkovCRFModel::chooseFirstState() const {
    // return _initial_probability->choose();
    return -1;
  }

  double SemiMarkovCRFModel::forward(const Sequence& sequence, Matrix& a) const {
    /*
    int nstates = _states.size();
    int size = sequence.size();
    Matrix alpha (nstates, size);

    for(int k = 0; k < nstates; k++)
      alpha(k,0) = _initial_probability->log_probability_of(k) + getState(k)->emission()->log_probability_of(sequence[0]);

    for(int t = 0; t < size-1; t++) {
      for(int i = 0; i < nstates; i++) {
        int j  = 0;
        if(j < nstates){
          alpha(i,t+1) =  alpha(j, t) + getState(j)->transitions()->log_probability_of(i);
          for( j = 1; j < nstates; j++)
            alpha(i, t+1) = log_sum_2(alpha(i,t+1), alpha(j, t) + getState(j)->transitions()->log_probability_of(i));
        }
        alpha(i,t+1) += getState(i)->emission()->log_probability_of(sequence[t+1]);
      }
    }
    a = alpha;
    double sum =  alpha(0, size-1);
    for(int k = 1; k < nstates; k++)
      sum = log_sum_2(sum, alpha(k, size-1));
    return sum;
    */
    return -1.0;
  }

  double SemiMarkovCRFModel::backward(const Sequence& sequence, Matrix& b) const {
    /*
    int nstates = _states.size();
    int size = sequence.size();
    Matrix beta (nstates, size);
    for(int k = 0; k < nstates; k++)
      beta(k,size-1) = 0.0;
    for(int t = size-2; t >= 0; t--) {
      for(int i = 0; i < nstates; i++) {
        int j = 0;
        if(j < nstates) {
          beta(i,t) =  getState(i)->transitions()->log_probability_of(j) + getState(j)->emission()->log_probability_of(sequence[t+1]) + beta(j, t+1);
          for(j = 1; j < nstates; j++)
            beta(i, t) = log_sum_2(beta(i, t), getState(i)->transitions()->log_probability_of(j) + getState(j)->emission()->log_probability_of(sequence[t+1]) + beta(j, t+1));
        }
      }
    }
    b = beta;
    double sum = -HUGE;
    for(int k = 0; k < nstates; k++)
      sum = log_sum_2(sum, beta(k, 0) + _initial_probability->log_probability_of(k) + getState(k)->emission()->log_probability_of(sequence[0]));
    return sum;
    */
   return -1.0;
  }

}
