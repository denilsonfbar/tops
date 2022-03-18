/*
 *       SemiMarkovCRFModel.hpp
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

#ifndef SEMI_MARKOV_CRF_MODEL_HPP
#define SEMI_MARKOV_CRF_MODEL_HPP

#include "crossplatform.hpp"

#include "ProbabilisticModel.hpp"
#include "DecodableModel.hpp"
#include "Sequence.hpp"
#include "Alphabet.hpp"
#include "ContextTree.hpp"

#include "SemiMarkovCRFModelCreator.hpp"
#include "util.hpp"
#include <cstdarg>
#include <vector>

namespace tops {

  class DLLEXPORT CRFState {

    protected:
      int _id;
      SymbolPtr _name;
      DiscreteIIDModelPtr _emission;
      DiscreteIIDModelPtr _transitions;
  
    public:
      CRFState(){}
      CRFState(int id, SymbolPtr name, DiscreteIIDModelPtr emission, DiscreteIIDModelPtr transitions) : _id(id), _name(name), _emission(emission), _transitions(transitions) {}
    
      void setId(int i) { _id = i; }
      void setName(SymbolPtr name) { _name = name; }
      void setEmissions(DiscreteIIDModelPtr e) { _emission = e; }
      void setTransition(DiscreteIIDModelPtr t) { _transitions = t; }
    
      int getId() { return _id; }
      SymbolPtr getName() const { return _name; }
      DiscreteIIDModelPtr& emission() { return _emission; }
      DiscreteIIDModelPtr& transitions() { return _transitions; }

      bool isSilent() {
        return (_emission == NULL);
      }
  };
  typedef boost::shared_ptr <CRFState> CRFStatePtr;



class DLLEXPORT FeatureFunction {
    public:
      FeatureFunction() {}
      virtual ~FeatureFunction() {}

      virtual double ff(int t, int y_t, int y_tm1, const Sequence& x) {
        std::cerr << "Not implemented: FeatureFunction.ff()" << std::endl;
        exit(-1);
        return 0.0;
      }
  };
  typedef boost::shared_ptr<FeatureFunction> FeatureFunctionPtr;

  class DLLEXPORT HMMEmissionFeatureFunction : public FeatureFunction {
    private:
      std::vector<CRFStatePtr> _states;
    public:
      HMMEmissionFeatureFunction(std::vector<CRFStatePtr> states) {
        _states = states;
      }
      virtual ~HMMEmissionFeatureFunction() {}

      virtual double ff(int t, int y_t, int y_tm1, const Sequence& x) {
        return _states[y_t]->emission()->log_probability_of(x[t]);
      }
  };
  typedef boost::shared_ptr<HMMEmissionFeatureFunction> HMMEmissionFeatureFunctionPtr;

  class DLLEXPORT HMMTransitionFeatureFunction : public FeatureFunction {
    private:
      std::vector<CRFStatePtr> _states;
      DiscreteIIDModelPtr _initial_probability;
    public:
      HMMTransitionFeatureFunction(std::vector<CRFStatePtr> states, DiscreteIIDModelPtr initial) {
        _states = states;
        _initial_probability = initial;
      }
      virtual ~HMMTransitionFeatureFunction() {}

      virtual double ff(int t, int y_t, int y_tm1, const Sequence& x) {
          if (t == 0){
            return _initial_probability->log_probability_of(y_t);
          } else {
            return _states[y_tm1]->transitions()->log_probability_of(y_t);
          }
      }
  };
  typedef boost::shared_ptr<HMMTransitionFeatureFunction> HMMTransitionFeatureFunctionPtr;

  class DLLEXPORT SemiMarkovCRFModel : public DecodableModel {
  
    private:
      std::vector <CRFStatePtr> _states;
      DiscreteIIDModelPtr _initial_probability;
      std::vector<double> _ctFactors;
      AlphabetPtr _state_names;
      std::vector<double> iterate(Sequence& obs);

      void scale(std::vector<double>& in, int t);

      std::vector<FeatureFunctionPtr> _features_functions;
      std::vector<double> _weights;
  
    public:
      SemiMarkovCRFModel() {}
      SemiMarkovCRFModel(std::vector <CRFStatePtr> states, DiscreteIIDModelPtr initial_probability, AlphabetPtr state_names, AlphabetPtr observation_symbols) : _states(states), _initial_probability(initial_probability), _state_names(state_names) {
        tops::ProbabilisticModel::setAlphabet(observation_symbols);
      }
      virtual ~SemiMarkovCRFModel() {}

      virtual ProbabilisticModelCreatorPtr getFactory() const {
        return SemiMarkovCRFModelCreatorPtr(new SemiMarkovCRFModelCreator());
      }

      virtual void initialize(const ProbabilisticModelParameters& par);
      void setStates(std::vector<CRFStatePtr> states, AlphabetPtr state_names);
      void setInitialProbability(DiscreteIIDModelPtr initial);
      void setObservationSymbols(AlphabetPtr obs) ;

      void setStates(std::vector<CRFStatePtr> states) { _states = states; }
      virtual void setState(int id, CRFStatePtr state) {
        if(_states.size() < _state_names->size())
          _states.resize(_state_names->size());
        _states[id] = state;
        state->setId(id);
      }

      virtual DecodableModel* decodable() { return this; }
      virtual std::string model_name() const { return "SemiMarkovCRFModel"; }
      virtual CRFStatePtr getState(int id) const { return _states[id]; }
      virtual std::string getStateName(int state) const;
      virtual AlphabetPtr getStateNames() const { return _state_names; }
      virtual ProbabilisticModelParameters parameters() const;
      virtual std::string str() const;

      //! Choose the observation given a state
      virtual Sequence& chooseObservation(Sequence& h, int i, int state) const;
      
      //! Choose a state
      virtual int chooseState(int state) const;
      //! Choose first state
      virtual int chooseFirstState() const;
    
      //! Forward algorithm
      virtual double forward(const Sequence& s, Matrix& alpha) const;

      //! Backward algorithm
      virtual double backward(const Sequence& s, Matrix& beta) const;

      //! Viterbi algorithm
      virtual double viterbiHMM(const Sequence& s, Sequence& path, Matrix& gamma) const;
    
      virtual void trainBaumWelch(SequenceList& training_set, int maxiterations, double diff) ;


      void configureCRF();

      //! Return the weighted sum of features functions results of a time t of a sequence x
      double sumFeatures(int t, int y_t, int y_tm1, const Sequence& x) const;

      //! Viterbi LCCRF algorithm
      virtual double viterbi(const Sequence& s, Sequence& path, Matrix& gamma) const;

      void print_matrix(Matrix& m) const;

  };
  typedef boost::shared_ptr<SemiMarkovCRFModel> SemiMarkovCRFModelPtr;

}

#endif
