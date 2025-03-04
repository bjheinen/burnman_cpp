/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
//--> test_mineral
// Check set_method
//  Auto
//  Custom
//  All others etc.
// Check set method --> custom EOS

// Check set_state
//  Property modifiers calculated
//  Only if they exist

// compute_molar_gibbs
//  matches eos
//  matches eos + excesses.G

// compute_molar_volume_unmodified
//  matches eos.volume

// compute_molar_volume
//  matches unmodified + excesses.dGdP

// compute_molar_entropy
//  matches + dGdT

// compute_isothermal_bulk_modulus_reuss
//  d2GdP2

// compute_molar_heat_capacity_p
//  -T * d2GdT2

// compute_thermal_expansivity
//  + d2GdPdT) / volume

// compute_shear_modulus
//  matches eos

// compute_molar_mass
//  check throw

// compute_density
//  check some known value

// compute_molar_internal_energy
//  check

// compute_molar_helmholtz
//  check

// compute_molar_enthalpy
//  check

// compute_isentropic_bulk_modulus_reuss
//  check

// compute_?_compressibility_reuss
//  1 / bulk modulus

// compute velocities
//  check known values

// compute_greuenseisen_parameter
//  check throw
//  check known value

// compute_molar_heat_capacity_v
//  check

// compute_isentropic_thermal_gradient
//  check

// Check via get() functions
// --> check caching?