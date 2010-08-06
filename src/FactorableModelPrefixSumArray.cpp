#include "FactorableModelPrefixSumArray.hpp"

namespace tops 
{

    
  //! Initialize the prefix sum array
  void FactorableModelPrefixSumArray::initialize(const Sequence & s) {
    if(_last == s)
      return;
    _last =s ;
    _alpha.resize(s.size() + 1);
    _precision.resize(s.size() + 1);
    _alpha[0] = 0;
    for(int i = 0 ; i < (int) s.size() ; i++) {
      double prob = _model->evaluate(s, i, i); 
      if(close(prob, -HUGE, 1e-1)) 
	{
	  _precision[i+1] = _precision[i]+1;
	  _alpha[i+1] = _model->evaluate(s,i,i);
	}
      else 
	{
	  _alpha[i+1] = _alpha[i] +   prob;
	  _precision[i+1] = _precision[i];
	}
    }
#if 0
    for(int i = 0; i < s.size(); i++)
      std::cerr << " "  << i << " " << m->alphabet()->getSymbol(s[i])->name() << " " << _alpha[i] << " " << _precision[i] << std::endl;
#endif 

  }
  
  //! compute the prefix sum array
  double FactorableModelPrefixSumArray::compute(int begin, int end ) const {
    if((begin >= 0) && ((end + 1) < (int) _alpha.size())) 
      {
	if((_precision[end+1] - _precision[begin]) > 0)
	  return -HUGE;
	return _alpha[end+1] - _alpha[begin];
      }
    else 
      return -HUGE;
  }
  

}