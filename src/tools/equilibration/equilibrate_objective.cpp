/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include "burnman/tools/equilibration/equilibrate_objective.hpp"
#include <cstddef>
#include <vector>
#include "burnman/tools/equilibration/equilibrate_utils.hpp"
#include "burnman/core/solution.hpp"

namespace burnman::equilibration {

Eigen::VectorXd F(
  const Eigen::VectorXd& x,
  Assemblage& assemblage,
  const ConstraintGroup& equality_constraints,
  const Eigen::VectorXd& reduced_composition_vector,
  const Eigen::MatrixXd& reduced_free_composition_vectors
) {
  // Update assemblage state
  set_composition_and_state_from_parameters(assemblage, x);
  // Retrieve updated endmember amounts
  Eigen::VectorXd new_endmember_amounts = get_endmember_amounts(assemblage);
  // Allocate F
  Eigen::Index n_eqc = static_cast<Eigen::Index>(equality_constraints.size());
  Eigen::VectorXd eqns = Eigen::VectorXd::Zero(assemblage.get_n_endmembers() + n_eqc);
  // Fill equality constraint portion of F
  for (Eigen::Index i = 0; i < n_eqc; ++i) {
    eqns(i) = equality_constraints[static_cast<std::size_t>(i)]->evaluate(x, assemblage);
  }
  // Compute reduced composition vector
  Eigen::VectorXd new_reduced_composition_vector = reduced_composition_vector;
  if (n_eqc > 2) {
    new_reduced_composition_vector +=
      x.tail(n_eqc - 2).transpose() * reduced_free_composition_vectors;
  }
  // TODO:: Assemblage::get_reaction_affinities
  Eigen::Index n_reac = assemblage.get_n_reactions();
  eqns.segment(n_eqc, n_reac) = assemblage.get_reaction_affinities();
  // TODO:: Assemblage::get_reduced_stoichiometric_matrix
  eqns.tail(eqns.size() - (n_eqc + n_reac)) = (
    assemblage.get_reduced_stoichiometric_matrix().transpose()
    * new_endmember_amounts
    ) - new_reduced_composition_vector;
  return eqns;
}

Eigen::MatrixXd J(
  const Eigen::VectorXd& x,
  Assemblage& assemblage,
  const ConstraintGroup& equality_constraints,
  const Eigen::MatrixXd& reduced_free_composition_vectors
) {
  Eigen::Index n_eqc = static_cast<Eigen::Index>(equality_constraints.size());
  Eigen::Index n_end = assemblage.get_n_endmembers();
  Eigen::Index jacobian_size = n_end + n_eqc;
  Eigen::MatrixXd jacobian = Eigen::MatrixXd::Zero(jacobian_size, jacobian_size);
  // Build constraints part of Jacobian
  for (Eigen::Index i = 0; i < n_eqc; ++i) {
    jacobian.row(i) =
      equality_constraints[static_cast<std::size_t>(i)]->derivative(
        x, assemblage, jacobian_size);
  }
  // P-T effects on each independent reaction
  // i.e. dF(i, reactions)/dx[0] and dF(i, reactions)/dx[1]
  Eigen::VectorXd partial_volumes = Eigen::VectorXd::Zero(n_end);
  Eigen::VectorXd partial_entropies = Eigen::VectorXd::Zero(n_end);
  std::vector<int> embr_per_phase = assemblage.get_endmembers_per_phase();
  Eigen::Index j = 0;
  for (std::size_t i = 0; i < embr_per_phase.size(); ++i) {
    Eigen::Index n = static_cast<Eigen::Index>(embr_per_phase[i]);
    if (n == 1) {
      // endmembers
      partial_volumes(j) = assemblage.get_phase(i)->get_molar_volume();
      partial_entropies(j) = assemblage.get_phase(i)->get_molar_entropy();
    } else {
      // solutions
      auto ph = assemblage.get_phase<Solution>(i);
      partial_volumes.segment(j, n) = ph->get_partial_volumes();
      partial_entropies.segment(j, n) = ph->get_partial_entropies();
    }
    j += n;
  }
  Eigen::VectorXd reaction_volumes = assemblage.get_reaction_basis() * partial_volumes;
  Eigen::VectorXd reaction_entropies = assemblage.get_reaction_basis() * partial_entropies;
  // dGi/dP = deltaVi; dGi/dT = -deltaSi
  jacobian.col(0).tail(n_end) = reaction_volumes;
  jacobian.col(1).tail(n_end) = -reaction_entropies;
  // Bulk composition constraints
  // P & T have no effect, i.e. dF(i, bulk)/dx[0] and dF(i, bulk)/dx[1] = 0
  // Build composition Hessian d2G/dfidfj = dmui/dfj
  // where fj is the fraction of endmember j in a phase.
  Eigen::VectorXd phase_amounts = assemblage.get_molar_fractions() * assemblage.get_n_moles();
  Eigen::MatrixXd comp_hessian = Eigen::MatrixXd::Zero(n_end, n_end);
  Eigen::MatrixXd dfi_dxj = Eigen::MatrixXd::Zero(n_end, n_end);
  Eigen::MatrixXd dpi_dxj = Eigen::MatrixXd::Zero(n_end, n_end);
  j = 0;
  for (std::size_t i = 0; i < embr_per_phase.size(); ++i) {
    Eigen::Index n = static_cast<Eigen::Index>(embr_per_phase[i]);
    if (n == 1) {
      // changing the amount (p) of a pure phase doesn't change its fraction in that phase
      dpi_dxj(j, j) = 1.0;
    } else {
      auto ph = assemblage.get_phase<Solution>(i);
      comp_hessian.block(j, j, n, n) = ph->get_gibbs_hessian();
      // x[0] = p(phase) and x[1:] = f[1:] - f[0]
      // Therefore
      // df[0]/dx[0] = 0
      // df[0]/dx[1:] = -1
      // (because changing the fraction of any endmember
      // depletes the fraction of the first endmember)
      // df[1:]/dx[1:] = 1 on diagonal, 0 otherwise
      // (because all other fractions are independent of each other)
      dfi_dxj.block(j, j, n, n) = Eigen::MatrixXd::Identity(n, n);
      dfi_dxj.row(j).segment(j, n).array() -= 1.0;
      // Total amounts of endmembers (p) are the fractions
      // multiplied by the amounts of their representative phases
      dpi_dxj.block(j, j, n, n) = dfi_dxj.block(j, j, n, n)
        * phase_amounts(static_cast<Eigen::Index>(i));
      // The derivative of the amount of each endmember with respect to
      // the amount of each phase is equal to the molar fractions of
      // the endmembers
      dpi_dxj.col(j).segment(j, n) = ph->get_molar_fractions();
    }
    j += n;
  }
  // dfi_dxj converts the endmember hessian to the parameter hessian.
  Eigen::MatrixXd reaction_hessian = assemblage.get_reaction_basis()
    * comp_hessian * dfi_dxj;
  Eigen::MatrixXd bulk_hessian = assemblage.get_reduced_stoichiometric_matrix().transpose()
    * dpi_dxj;
  if (reaction_hessian.rows() > 0) {
    jacobian.block(
      n_eqc, 2,
      reaction_hessian.rows(), reaction_hessian.cols()
    ) = reaction_hessian;
    jacobian.block(
      n_eqc + reaction_hessian.rows(),
      2, bulk_hessian.rows(), bulk_hessian.cols()
    ) = bulk_hessian;
  } else {
    jacobian.block(
      n_eqc, 2,
      bulk_hessian.rows(), bulk_hessian.cols()
    ) = bulk_hessian;
  }
  if (reduced_free_composition_vectors.rows() > 0) {
    // Swap because transposing
    Eigen::Index nrows = reduced_free_composition_vectors.cols();
    Eigen::Index ncols = reduced_free_composition_vectors.rows();
    jacobian.block(
      jacobian.rows() - nrows, 2 + reaction_hessian.cols(),
      nrows, ncols
    ) = -reduced_free_composition_vectors.transpose();
  }
  return jacobian;
}

} // namespace burnman::equilibration
