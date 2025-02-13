/*
  TODO: Copyright Notice!
*/
#ifndef BURNMAN_EOS_BIRCH_MURNAGHAN_HPP_INCLUDED
#define BURNMAN_EOS_BIRCH_MURNAGHAN_HPP_INCLUDED

#include "burnman/core/equation_of_state.hpp"
#include "burnman/utils/eos.hpp"

/**
 * @class BM3
 * @brief Base class for the Birch Murnaghan equation of state.
 *
 * The EOS is third order in strain, and has no temperature dependence.
 * Uses third order expansion also for shear modulus. To fit shear modulus
 * to second order expansion @see BM2
 * 
 * @note All functions assume SI units for all properties.
 */
class BM3 : public EquationOfState{
 public:

  // Helper functions
  bool validate_parameters(MineralParams& params) override;

  // Specific EOS functions
  double compute_volume(
    double pressure,
    double temperature,
    const MineralParams& params) const override;

  double compute_pressure(
    double temperature,
    double volume,
    const MineralParams& params) const override;

  double compute_grueneisen_parameter(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

  double compute_isothermal_bulk_modulus_reuss(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

  double compute_isentropic_bulk_modulus_reuss(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

  double compute_shear_modulus(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

  double compute_molar_heat_capacity_v(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

  double compute_molar_heat_capacity_p(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

  double compute_thermal_expansivity(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

  double compute_gibbs_free_energy(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

  double compute_entropy(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

  double compute_molar_internal_energy(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

 private:

  /**
   * @brief Evaluate the BM EOS pressure.
   *
   * @param inv_compression V_0/V.
   * @param params Mineral parameters object of type MineralParams
   *
   * @return Pressure in [Pa].
   */
  static double compute_birch_murnaghan(
    double inv_compression,
    const MineralParams& params);

  /**
   * @brief GSL function wrapper to compute P(V) - P
   * 
   * @param x Volume to test (passed by solver)
   * @param p Generic pointer for parameter object
   * @see `ParamsGSL::SolverParams_P`
   */
  static double bm_gsl_wrapper(double x, void* p);

};

#endif // BURNMAN_EOS_BIRCH_MURNAGHAN_HPP_INCLUDED