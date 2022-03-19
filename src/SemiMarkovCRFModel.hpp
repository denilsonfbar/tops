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
#include "GeneralizedHiddenMarkovModel.hpp"
#include "Sequence.hpp"
#include "Alphabet.hpp"
#include "ContextTree.hpp"

#include "SemiMarkovCRFModelCreator.hpp"
#include "util.hpp"
#include <cstdarg>
#include <vector>

namespace tops {

  class DLLEXPORT FeatureFunction {
    public:
      FeatureFunction() {}
      virtual ~FeatureFunction() {}

      virtual double ff(int y_tm1, int t, int b_t, int e_t, int y_t, const Sequence& x) const = 0;
  };
  typedef boost::shared_ptr<FeatureFunction> FeatureFunctionPtr;


  class DLLEXPORT FeatureFunctionGHMMTransitions : public FeatureFunction {
    private:
      GeneralizedHiddenMarkovModelPtr _ghmm_model;
    public:
      FeatureFunctionGHMMTransitions(GeneralizedHiddenMarkovModelPtr ghmm_model) : _ghmm_model(ghmm_model) {}
      virtual ~FeatureFunctionGHMMTransitions() {}

      virtual double ff(int y_tm1, int t, int b_t, int e_t, int y_t, const Sequence& x) {
          if (t == 0)
            return _ghmm_model->getInitialProbabilities()->log_probability_of(y_t);
          else
            return _ghmm_model->getAllStates()[y_tm1]->transition()->log_probability_of(y_t);
      }
  };
  typedef boost::shared_ptr<FeatureFunctionGHMMTransitions> FeatureFunctionGHMMTransitionsPtr;


  class DLLEXPORT FeatureFunctionGHMMObservations : public FeatureFunction {
    private:
      GeneralizedHiddenMarkovModelPtr _ghmm_model;
    public:
      FeatureFunctionGHMMObservations(GeneralizedHiddenMarkovModelPtr ghmm_model) : _ghmm_model(ghmm_model) {}
      virtual ~FeatureFunctionGHMMObservations() {}

      virtual double ff(int y_tm1, int t, int b_t, int e_t, int y_t, const Sequence& x) {
        return _ghmm_model->getAllStates()[y_t]->observation()->log_probability_of(x[t]);
      }
  };
  typedef boost::shared_ptr<FeatureFunctionGHMMObservations> FeatureFunctionGHMMObservationsPtr;


  class DLLEXPORT FeatureFunctionGHMMDurations : public FeatureFunction {
    private:
      GeneralizedHiddenMarkovModelPtr _ghmm_model;
    public:
      FeatureFunctionGHMMDurations(GeneralizedHiddenMarkovModelPtr ghmm_model) : _ghmm_model(ghmm_model) {}
      virtual ~FeatureFunctionGHMMDurations() {}

      virtual double ff(int y_tm1, int t, int b_t, int e_t, int y_t, const Sequence& x) {
        int d = e_t - b_t + 1;
        return _ghmm_model->getAllStates()[y_t]->duration_probability(d);
      }
  };
  typedef boost::shared_ptr<FeatureFunctionGHMMDurations> FeatureFunctionGHMMDurationsPtr;


  class DLLEXPORT SemiMarkovCRFModel : public DecodableModel {
    private:
      std::vector<FeatureFunctionPtr> _features_functions;
      std::vector<double> _weights;
  
    public:
      SemiMarkovCRFModel() {}
      virtual ~SemiMarkovCRFModel() {}

      virtual ProbabilisticModelCreatorPtr getFactory() const {
        return SemiMarkovCRFModelCreatorPtr(new SemiMarkovCRFModelCreator());
      }
      virtual DecodableModel* decodable() { return this; }
      virtual std::string model_name() const { return "SemiMarkovCRFModel"; }

      virtual ProbabilisticModelParameters parameters() const;
      virtual std::string str() const;
      virtual std::string getStateName(int state) const;
      virtual AlphabetPtr getStateNames() const;

      virtual void initialize(const ProbabilisticModelParameters& par);

      //! Choose the observation given a state
      virtual Sequence& chooseObservation(Sequence& h, int i, int state) const;
      virtual int chooseState(int state) const;
      virtual int chooseFirstState() const;

      virtual double forward(const Sequence& s, Matrix& alpha) const;
      virtual double backward(const Sequence& s, Matrix& beta) const;
      virtual double viterbi(const Sequence& s, Sequence& path, Matrix& gamma) const;

/*
      virtual void trainBaumWelch(SequenceList& training_set, int maxiterations, double diff) ;

      void configureCRF();

      //! Return the weighted sum of features functions results of a time t of a sequence x
      double sumFeatures(int t, int y_t, int y_tm1, const Sequence& x) const;

      //! Viterbi LCCRF algorithm
      virtual double viterbi(const Sequence& s, Sequence& path, Matrix& gamma) const;

      void print_matrix(Matrix& m) const;
*/

  };
  typedef boost::shared_ptr<SemiMarkovCRFModel> SemiMarkovCRFModelPtr;

}

#endif
