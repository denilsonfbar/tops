/*
 *       SemiMarkovCRFModelCreator.cpp
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

#include "SemiMarkovCRFModelCreator.hpp"
#include "DiscreteIIDModelCreator.hpp"
#include "DiscreteIIDModel.hpp"
#include "ConfigurationReader.hpp"
#include "SemiMarkovCRFModel.hpp"
#include "Alphabet.hpp"
#include "Symbol.hpp"
namespace tops {
  ProbabilisticModelPtr SemiMarkovCRFModelCreator::create(ProbabilisticModelParameters & parameters) const {

    ProbabilisticModelParameterValuePtr state_names = parameters.getMandatoryParameterValue("state_names");
    ProbabilisticModelParameterValuePtr observation_symbols = parameters.getMandatoryParameterValue("observation_symbols");

    if((state_names == NULL)||
       (observation_symbols == NULL)) {
        std::cerr << help() << std::endl;
      }

    ProbabilisticModelPtr model = SemiMarkovCRFModelPtr(new SemiMarkovCRFModel());
    model->initialize(parameters);
    return model;
  }
}
