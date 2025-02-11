/*
  TODO: Copyright Notice!
*/
#ifndef BURNMAN_EOS_VINET_HPP_INCLUDED
#define BURNMAN_EOS_VINET_HPP_INCLUDED

#include "./equation_of_state.hpp"
#include "../util/eos.hpp"


/**
 * @class Vinet
 * @brief Base class for the isothermal Vinet equation of state.
 *
 * References for this equation of state are Vinet (1986) and Vinet (1987).
 * This equation of state actually predates Vinet by 55 years (Rydberg, 1932)
 * and was investigated further by Stacey (1981).
 * Also called the Rose-Vinet, Vinet-Rydberg and Morse-Rydberg EOS.
 *
 * @note All functions assume SI units for all properties.
 */
class Vinet : public EquationOfState{
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
   * @brief Evaluate the Vinet EOS pressure.
   *
   * @param compression V/V_0.
   * @param params Mineral parameters object of type MineralParams
   *
   * @return Pressure in [Pa].
   */
  double compute_vinet(double compression, const MineralParams& params) const;

  /**
   * @brief GSL function wrapper to compute P(V) - P
   * 
   * @param x Volume to test (passed by solver)
   * @param p Generic pointer for parameter object
   * @see `ParamsGSL::SolverParams_P`
   */
  double vinet_gsl_wrapper(double x, void* p); //const void* p) const; ?
};

#endif // BURNMAN_EOS_VINET_HPP_INCLUDED