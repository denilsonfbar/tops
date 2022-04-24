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
#include "ProbabilisticModelCreatorClient.hpp"
#include "util.hpp"
#include <cstdarg>
#include <vector>


#include "FeatureFunctionGHMM.hpp"


namespace tops {

  // TODO: fix this model (methods relocated to GeneralizedHiddenMarkovModel)
  class DLLEXPORT SemiMarkovCRFModel : public DecodableModel {
    private:
      int _n_states;
      AlphabetPtr _state_names;
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

      virtual void setStates(AlphabetPtr state_names);
      virtual void setObservationSymbols(AlphabetPtr obs);

      //! Return the weighted sum of features functions results of a time t of a sequence x
      double sumFeatures(int y_tm1, int t, int b_t, int e_t, int y_t, const Sequence& x) const;

      virtual double viterbi(const Sequence& s, Sequence& path, Matrix& gamma) const;

      void print_matrix(Matrix& m) const;


      //! Choose the observation given a state
      virtual Sequence& chooseObservation(Sequence& h, int i, int state) const;
      virtual int chooseState(int state) const;
      virtual int chooseFirstState() const;

      virtual double forward(const Sequence& s, Matrix& alpha) const;
      virtual double backward(const Sequence& s, Matrix& beta) const;
  };
  typedef boost::shared_ptr<SemiMarkovCRFModel> SemiMarkovCRFModelPtr;

}

#endif
