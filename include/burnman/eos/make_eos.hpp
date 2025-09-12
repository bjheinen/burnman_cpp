/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_EOS_MAKE_EOS_HPP_INCLUDED
#define BURNMAN_EOS_MAKE_EOS_HPP_INCLUDED

#include <memory>
#include "burnman/utils/types/simple_types.hpp"
#include "burnman/core/equation_of_state.hpp"

namespace burnman {
  namespace eos {
    /**
    * @brief Makes a pointer to a predefined EOS class
    * @param eos_type types::EOSType enum specifying equation of state
    * @return std::shared_ptr to an instance of the specific EOS
    */
    std::shared_ptr<EquationOfState> make_eos(types::EOSType eos_type);
  } // namespace eos
} // namespace burnman

#endif // BURNMAN_EOS_MAKE_EOS_HPP_INCLUDED
