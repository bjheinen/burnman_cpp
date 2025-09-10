/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include "burnman/eos/make_eos.hpp"
#include <stdexcept>
#include "burnman/eos/vinet.hpp"
#include "burnman/eos/birch_murnaghan.hpp"
#include "burnman/eos/mie_grueneisen_debye.hpp"
#include "burnman/eos/slb.hpp"

// TODO: Add HP types here

std::shared_ptr<EquationOfState> make_eos(EOSType eos_type) {
  switch (eos_type) {
    case EOSType::Vinet:
      return std::make_shared<Vinet>();
    case EOSType::BM2:
      return std::make_shared<BM2>();
    case EOSType::BM3:
      return std::make_shared<BM3>();
    case EOSType::MGD2:
      return std::make_shared<MGD2>();
    case EOSType::MGD3:
      return std::make_shared<MGD3>();
    case EOSType::SLB2:
      return std::make_shared<SLB2>();
    case EOSType::SLB3:
      return std::make_shared<SLB3>();
    case EOSType::SLB3Conductive:
      return std::make_shared<SLB3Conductive>();
    default:
      throw std::invalid_argument("Unknown EOS type!");
  }
}
