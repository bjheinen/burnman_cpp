/*
  TODO: Copyright Notice!
*/
#include "../../include/burnman/classes/mineral.hpp"


// EOS properties - in P,T form

double Mineral::compute_molar_gibbs() const {
  // calls method.gibbs_free_energy()
  // adds ._property_modifiers['G']
}

// TODO: implement _property_modifiers
//       implement set_method so know how to call EOS funcs
//         should EOS be functional to reduce overheads??


// Mineral parameter properties