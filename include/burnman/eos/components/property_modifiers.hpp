/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_EOS_PROPERTY_MODIFIERS_HPP_INCLUDED
#define BURNMAN_EOS_PROPERTY_MODIFIERS_HPP_INCLUDED

#include "burnman/utils/types/excess_params.hpp"

/**
 * TODO: Add Docs to EOS module
 * @brief Modification functions for thermodynamic properties.
 *
 * Currently this includes modifications for:
 *   - second order transitions (landau, landau_slb_2022, landau_hp)
 *   - order-disorder (bragg_williams)
 *   - magnetism (magnetic_chs)
 *   - linear
 */
namespace excesses {

  /**
   * @brief Computes a tricritical Landau correction.
   *
   * Used to apply a tricritical Landau correction to the properties
   * of an endmember which undergoes a displacive phase transition.
   * Correction follows Putnis (1992), and is done relative to the
   * completely *ordered* state (at 0 K).
   * This version differs in implemenation from Stixrude & 
   * Lithgow-Bertelloni (2022) and Holland and Powell (2011), who
   * both compute properties relative to the completely disordered state
   * and standard state respectively.
   *
   * The excess entropy (and heat capacity) terms are equal to zero at 0 K.
   *
   * N.B. The excesses are for a *completely relaxed* mineral; for example,
   * seismic wave propagation is *slow compared to the rate of change in 
   * order parameter.
   *
   * @param pressure In [Pa].
   * @param temperature In [K].
   * @param params Function specific parameters.
   *
   * @return excesses:Excesses object.
   */
  Excesses compute_excesses(
    double pressure,
    double temperature,
    LandauParams params);

  /**
   * @brief Computes a tricritical Landau correction.
   *
   * Used to apply a tricritical Landau correction to the properties
   * of an endmember which undergoes a displacive phase transition.
   * Correction follows Stixrude and Lithgow-Bertelloni (2022), and is 
   * done relative to the state with order parameter Q=1.
   *
   * The order parameter of this formulation can exceed one, at odds with
   * Putnis (above), but in better agreement with atomic intuition
   * Nevertheless, this implementation is still not perfect, as the excess
   * entropy (and heat capacity) terms are not equal to zero at 0 K.
   * Q is limited to values less than or equal to 2 to avoid unrealistic
   * stabilisation at ultrahigh pressure.
   *
   * N.B. The excesses are for a *completely relaxed* mineral; for example,
   * seismic wave propagation is *slow compared to the rate of change in 
   * order parameter.
   *
   * @param pressure In [Pa].
   * @param temperature In [K].
   * @param params Function specific parameters.
   *
   * @return excesses:Excesses object.
   */
  Excesses compute_excesses(
    double pressure,
    double temperature,
    LandauSLB2022Params params);

  /**
   * @brief Computes a tricritical Landau correction.
   *
   * Used to apply a tricritical Landau correction to the properties
   * of an endmember which undergoes a displacive phase transition.
   * Correction follows Holland and Powell (2011), and is done relative
   * to the standard state.
   *
   * Includes the correction published within landaunote.pdf
   * (Holland, pers. comm), which 'corrects' the terms involving
   * the critical temperature Tc / Tc*
   *
   * This formalism predicts that the order parameter can be greate than one.
   *
   * N.B. The excesses are for a *completely relaxed* mineral; for example,
   * seismic wave propagation is *slow compared to the rate of change in 
   * order parameter.
   *
   * @param pressure In [Pa].
   * @param temperature In [K].
   * @param params Function specific parameters.
   *
   * @return excesses:Excesses object.
   */
  Excesses compute_excesses(
    double pressure,
    double temperature,
    LandauHPParams params);

  /**
   * @brief Computes a linear property correction.
   *
   * Applies a 'Darken's quadratic formalism' correction (Powell, 1987)
   * to the thermodynamic properties of a mineral endmember.
   * This correction is relative to P = 0 and T = 0 and linear in P and T
   * and therefore corresponds to a constant volume and entropy correction.
   *
   * Applying either a volume or entropy term will generally break
   * equations of state (i.e. the properties of the mineral will no longer
   * obey the equation of state defined in the params. However, this form of
   * excess is extremely useful as a first order tweak to free energies
   * (especially in solid solution calculations).
   *
   * @param pressure In [Pa].
   * @param temperature In [K].
   * @param params Function specific parameters.
   *
   * @return excesses:Excesses object.
   */
  Excesses compute_excesses(
    double pressure,
    double temperature,
    LinearParams params);

  /**
   * @brief Computes a Bragg-Williams-type correction for order-disorder.
   *
   * Applies a Bragg-Williams type correction to the thermodynamic
   * properties of a mineral endmember. 
   * Used for modelling order-disorder processes.
   * Expressions are from Holland and Powell (1996).
   *
   * N.B. The excesses are for a *completely relaxed* mineral; i.e. the 
   * seismic wave propagation is *slow* compared to the rate of reaction.
   *
   * This may not be reasonable for order-disorder, especially
   * for slow or coupled diffusers (Si-Al, for example).
   * The completely *unrelaxed* mineral (in terms of order-disorder)
   * can be calculated with a solid solution model.
   *
   * @param pressure In [Pa].
   * @param temperature In [K].
   * @param params Function specific parameters.
   *
   * @return excesses:Excesses object.
   */
  Excesses compute_excesses(
    double pressure,
    double temperature,
    BraggWilliamsParams params);

  /**
   * @brief Computes a magnetic excess contribution.
   *
   * Applies a magnetic contribution to the thermodynamic properties of
   * a mineral endmember.
   * The expression for the gibbs energy contribution is from Chin,
   * Hertzmann and Sundman (1987) and Sundman (1991).
   *
   * @param pressure In [Pa].
   * @param temperature In [K].
   * @param params Function specific parameters.
   *
   * @return excesses:Excesses object.
   */
  Excesses compute_excesses(
    double pressure,
    double temperature,
    MagneticChsParams params);

  /**
   * @brief Excess correction based on Debye model.
   *
   * Applies an excess contribution based on a Debye model. The excess
   * heat capacity tends toward a constant value at high temperature.
   *
   * @param pressure In [Pa].
   * @param temperature In [K].
   * @param params Function specific parameters.
   *
   * @return excesses:Excesses object.
   */
  Excesses compute_excesses(
    double pressure,
    double temperature,
    DebyeParams params);

  /**
   * @brief Excess correction based on thermal derivatives of Debye model.
   *
   * Applies an excess contribution based on thermal derivatives of a Debye
   * model. The excess entropy tends toward a constant value at high
   * temperature and behaves like the heat capacity of a Debye model at
   * finite temperature.
   *
   * @param pressure In [Pa].
   * @param temperature In [K].
   * @param params Function specific parameters.
   *
   * @return excesses:Excesses object.
   */
  Excesses compute_excesses(
    double pressure,
    double temperature,
    DebyeDeltaParams params);

  /**
   * @brief Excess correction based on an Einstein model.
   *
   * Applies an excess contribution based on Einstein model. The excess
   * heat capacity tends toward a constant value at high temperature.
   *
   * @param pressure In [Pa].
   * @param temperature In [K].
   * @param params Function specific parameters.
   *
   * @return excesses:Excesses object.
   */
  Excesses compute_excesses(
    double pressure,
    double temperature,
    EinsteinParams params);

  /**
   * @brief Excess correction based on thermal derivatives of Einstein model.
   *
   * Applies an excess contribution based on thermal derivatives of an
   * Einstein model. The excess entropy tends toward a constant value at
   * high temperature and behaves like the heat capacity of an Einstein
   * model at finite temperature.
   *
   * @param pressure In [Pa].
   * @param temperature In [K].
   * @param params Function specific parameters.
   *
   * @return excesses:Excesses object.
   */
  Excesses compute_excesses(
    double pressure,
    double temperature,
    EinsteinDeltaParams params);

} // End namespace excesses

#endif // BURNMAN_EOS_PROPERTY_MODIFIERS_HPP_INCLUDED
