/*
  TODO: Copyright Notice!
*/
#ifndef BURNMAN_UTILS_MAKE_EOS_HPP_INCLUDED
#define BURNMAN_UTILS_MAKE_EOS_HPP_INCLUDED

#include <memory>
#include <stdexcept>
#include "burnman/utils/eos.hpp"
#include "burnman/core/equation_of_state.hpp"
#include "burnman/eos/vinet.hpp"
#include "burnman/eos/birch_murnaghan.hpp"
#include "burnman/eos/mie_grueneisen_debye.hpp"

/**
 * @brief Makes a pointer to a predefined EOS class
 * @param eos_type EOSType enum specifying equation of state
 * @return std::unique_ptr to an instance of the specific EOS
 */
std::unique_ptr<EquationOfState> make_eos(EOSType eos_type) {
  switch (eos_type) {
    case EOSType::Vinet:
      return std::make_unique<Vinet>();
    case EOSType::BM2:
      return std::make_unique<BM2>();
    case EOSType::BM3:
      return std::make_unique<BM3>();
    case EOSType::MGD2:
      return std::make_unique<MGD2>();
    case EOSType::MGD3:
      return std::make_unique<MGD3>();
    default:
      throw std::invalid_argument("Unknown EOS type!");
  }
}

#endif // BURNMAN_UTILS_MAKE_EOS_HPP_INCLUDED