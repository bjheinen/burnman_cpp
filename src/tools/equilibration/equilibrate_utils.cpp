/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include "burnman/tools/equilibration/equilibrate_utils.hpp"
#include <cassert>
#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include "burnman/core/solution.hpp"

EquilibrationParameters get_equilibration_parameters(
  const Assemblage& assemblage,
  const FormulaMap& composition,
  const std::vector<std::unordered_map<std::string, double>>& free_compositional_vectors
) {
  // Make parameter object
  EquilibrationParameters prm;
  // Make temporary for storing phase_amount_indices
  std::vector<int> phase_amount_ind_temp;
  // Process parameter names
  prm.parameter_names.push_back("Pressure (Pa)");
  prm.parameter_names.push_back("Temperature (K)");
  int embr_start_idx = 0;
  std::vector<int> embr_per_phase = assemblage.get_endmembers_per_phase();
  std::vector<std::string> embr_names = assemblage.get_endmember_names();
  for (std::size_t i = 0; i < embr_per_phase.size(); ++i) {
    auto ph = assemblage.get_phase(i);
    phase_amount_ind_temp.push_back(static_cast<int>(prm.parameter_names.size()));
    prm.parameter_names.push_back("x(" + ph->get_name() + ")");
    int n_mbrs = embr_per_phase[i];
    // When n_mbrs > 1, add embr names (but skip 1st)
    for (int j = 1; j < n_mbrs; ++j) {
      prm.parameter_names.push_back("p(" + embr_names[static_cast<std::size_t>(embr_start_idx + j)] + ")");
    }
    embr_start_idx += n_mbrs;
  }
  // Add parameter names for any free_compositional_vectors
  std::size_t n_free_compositional_vectors = free_compositional_vectors.size();
  for (std::size_t i = 0; i < n_free_compositional_vectors; ++i) {
    prm.parameter_names.push_back("v_" + std::to_string(i));
  }
  // Get number of parameters
  prm.n_parameters = static_cast<Eigen::Index>(prm.parameter_names.size());
  // Map phase amount indices
  prm.phase_amount_indices = Eigen::Map<Eigen::ArrayXi>(
    phase_amount_ind_temp.data(), phase_amount_ind_temp.size());
  // Process bulk composition vector
  const std::vector<std::string>& elements = assemblage.get_elements();
  Eigen::Index n_elements = static_cast<Eigen::Index>(elements.size());
  prm.bulk_composition_vector.resize(n_elements);
  for (Eigen::Index i = 0; i < n_elements; ++i) {
    std::string el = elements[static_cast<std::size_t>(i)];
    auto it = composition.find(el);
    if (it == composition.end()) {
      throw std::runtime_error(
        "Element '" + el + "' not found in bulk composition!"
      );
    }
    prm.bulk_composition_vector(i) = it->second;
  }
  // Process free_compositional_vectors
  if (n_free_compositional_vectors > 0) {
    prm.free_compositional_vectors.resize(
      static_cast<Eigen::Index>(n_free_compositional_vectors),
      n_elements
    );
    for (std::size_t i = 0; i < n_free_compositional_vectors; ++i) {
      for (std::size_t j = 0; j < static_cast<std::size_t>(n_elements); ++j) {
        const std::string& el = elements[j];
        auto it = free_compositional_vectors[i].find(el);
        if (it == free_compositional_vectors[i].end()) {
          throw std::runtime_error(
            "Element '" + el + "' not found in composition!"
          );
        }
        prm.free_compositional_vectors(
          static_cast<Eigen::Index>(i),
          static_cast<Eigen::Index>(j)) = it->second;
      }
    }
  } else {
    // Resize but keep empty
    prm.free_compositional_vectors.resize(0, n_elements);
  }

  // Check bulk composition
  if (assemblage.get_compositional_null_basis().rows() != 0) {
    // In burnman python only first element is checked
    // assemblage.compositional_null_basis.dot(prm.bulk_composition_vector)[0]
    double val = (assemblage.get_compositional_null_basis() * prm.bulk_composition_vector)(0);
    if (std::abs(val) > constants::precision::abs_tolerance) {
      throw std::runtime_error(
        "The bulk composition is not within the compositional space of the assemblage."
      );
    }
    // TODO: Should we check if whole Vector is zero??
    // i.e.
    //Eigen::VectorXd v = assemblage.get_compositional_null_basis() * prm.bulk_composition_vector;
    //if (!v.isZero(constants::precision::abs_tolerance)) {
    //  throw ...
    //}
  }
  // Reduce vector to independent elements
  const auto& indep = assemblage.get_independent_element_indices();
  prm.reduced_composition_vector = prm.bulk_composition_vector(indep);
  prm.reduced_free_compositional_vectors = prm.free_compositional_vectors(Eigen::all, indep);
  // Process constraints
  auto [constraint_matrix, constraint_vector] = calculate_constraints(assemblage, static_cast<int>(n_free_compositional_vectors));
  prm.constraint_matrix = constraint_matrix;
  prm.constraint_vector = constraint_vector;
  return prm;
}

std::pair<Eigen::MatrixXd, Eigen::VectorXd> calculate_constraints(
  const Assemblage& assemblage,
  int n_free_compositional_vectors
) {
  // Use std::optional for empty bounds
  std::vector<std::optional<Eigen::ArrayXXd>> bounds;
  Eigen::Index n_constraints = 0;
  std::vector<int> embr_per_phase = assemblage.get_endmembers_per_phase();
  for (std::size_t i = 0; i < embr_per_phase.size(); ++i) {
    std::optional<Eigen::ArrayXXd> bound;
    if (embr_per_phase[i] > 1) {
      bound = assemblage.get_phase<Solution>(i)->get_endmember_occupancies();
      n_constraints += bound->cols(); // n_elements
    }
    bounds.push_back(bound);
    ++n_constraints;
  }
  // Setup of vector/matrix
  Eigen::VectorXd c_vector = Eigen::VectorXd::Zero(n_constraints + 2);
  Eigen::MatrixXd c_matrix = Eigen::MatrixXd::Zero(
    n_constraints + 2,
    assemblage.get_n_endmembers() + 2 + n_free_compositional_vectors
  );
  // Manually set P T constraints
  c_matrix(0, 0) = -1.0; // P > 0
  c_matrix(1, 1) = -1.0; // T > 0
  Eigen::Index cidx = 2; // Index of current compositional constraint (row)
  Eigen::Index pidx = 0; // Starting index of current phase
  for (std::size_t i = 0; i < embr_per_phase.size(); ++i) {
    Eigen::Index n = static_cast<Eigen::Index>(embr_per_phase[i]);
    c_matrix(cidx, pidx + 2) = -1.0; // phase prop > 0
    // The first endmember proportion is not a free variable
    // (all endmembers in solution must sum to one)
    // Re-express the constraints without the first endmember
    ++cidx;
    if (bounds[i].has_value()) {
      const Eigen::ArrayXXd& occ = *(bounds[i]);
      Eigen::Index m = occ.cols();
      c_vector.segment(cidx, m) = -occ.row(0);
      c_matrix.block(cidx, pidx + 3, m, n - 1) =
        (occ.row(0).transpose().matrix() * Eigen::VectorXd::Ones(n - 1).transpose())
        - occ.transpose().block(0, 1, m, n-1).matrix();
      cidx += m;
    }
    pidx += n;
  }
  return {c_matrix, c_vector};
}

Eigen::VectorXd get_parameter_vector(
  const Assemblage& assemblage,
  int n_free_compositional_vectors
) {
  int n_params = assemblage.get_n_endmembers() + 2 + n_free_compositional_vectors;
  Eigen::VectorXd params = Eigen::VectorXd::Zero(n_params);
  Eigen::ArrayXd n_moles_per_phase = assemblage.n_moles * assemblage.molar_fractions;
  // check pressure & temperature are set
  if (!assemblage.has_state()) {
    throw std::runtime_error("You need to set_state before getting parameters");
  }
  // Set P & T as first two params
  params(0) = assemblage.pressure;
  params(1) = assemblage.temperature;
  std::vector<int> embr_per_phase = assemblage.get_endmembers_per_phase();
  Eigen::Index j = 2;
  for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(assemblage.get_n_phases()); ++i) {
    params(j) = n_moles_per_phase(i);
    if (auto ph = assemblage.get_phase<Solution>(static_cast<std::size_t>(i))) {
      Eigen::Index n_embr = static_cast<Eigen::Index>(embr_per_phase[i] - 1); // skip first embr
      params.segment(j + 1, n_embr) = ph->get_molar_fractions().tail(n_embr);
    }
    j += embr_per_phase[static_cast<std::size_t>(i)];
  }
  return params;
}

Eigen::VectorXd get_endmember_amounts(
  const Assemblage& assemblage
) {
  Eigen::ArrayXd phase_amounts = assemblage.get_n_moles * assemblage.get_molar_fractions();
  Eigen::VectorXd abs_amounts(assemblage.get_n_endmembers());
  std::vector<int> embr_per_phase = assemblage.get_endmembers_per_phase();
  Eigen::Index j = 0;
  for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(assemblage.get_n_phases()); ++i) {
    if (auto ph = assemblage.get_phase<Solution>(static_cast<std::size_t>(i))) {
      abs_amounts.segment(j, j + static_cast<Eigen::Index>(embr_per_phase(i))) =
        phase_amounts(i) * ph->get_molar_fractions();
    } else {
      abs_amounts(j) = phase_amounts(i);
    }
    j += embr_per_phase[static_cast<std::size_t>(i)];
  }
  return abs_amounts;
}

void set_composition_and_state_from_parameters(
  Assemblage& assemblage,
  const Eigen::VectorXd& parameters
) {
  // Set P & T (first two parameters)
  assemblage.set_state(parameters(0), parameters(1));
  Eigen::Index n_phases = static_cast<Eigen::Index>(assemblage.get_n_phases());
  Eigen::ArrayXd phase_amounts = Eigen::ArrayXd::Zero(n_phases);
  Eigen::Index i = 2;
  for (Eigen::Index phase_idx = 0; phase_idx < n_phases; ++phase_idx) {
    phase_amounts(phase_idx) = parameters(i);
    // TODO: get_phase take Eigen::Index?
    if (auto ph = assemblage.get_phase<Solution>(static_cast<std::size_t>(phase_idx))) {
      Eigen::Index n_mbrs = static_cast<Eigen::Index>(ph->get_n_endmembers());
      Eigen::ArrayXd f = Eigen::ArrayXd::Zero(n_mbrs);
      f.segment(1, n_mbrs - 1) = parameters.segment(i + 1, n_mbrs - 1);
      f(0) = 1.0 - f.tail(n_mbrs - 1).sum();
      ph->set_composition(f);
      i += n_mbrs;
    } else {
      ++i;
    }
  }
  assert((phase_amounts > -1.0e-8).all());
  phase_amounts = phase_amounts.abs();
  assemblage.set_n_moles = phase_amounts.sum();
  assemblage.set_fractions(phase_amounts / assemblage.get_n_moles());
}
