/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include <algorithm>
#include <cassert>
#include <cmath>
#include <Eigen/Dense>
#include "burnman/tools/equilibration/equilibrate.hpp"

Eigen::VectorXd get_parameter_vector(
  const Assemblage& assemblage,
  int n_free_compositional_vectors = 0
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
      params.segment(j + 1, n_embr) = ph.get_molar_fractions().tail(n_embr);
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
        phase_amounts(i) * ph.get_molar_fractions();
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
      Eigen::Index n_mbrs = static_cast<Eigen::Index>(ph.get_n_endmembers());
      Eigen::ArrayXd f = Eigen::ArrayXd::Zero(n_mbrs);
      f.segment(1, n_mbrs - 1) = parameters.segment(i + 1, n_mbrs - 1);
      f(0) = 1.0 - f.tail(n_mbrs - 1).sum();
      ph.set_composition(f);
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

std::pair<double, double> lambda_bounds_func(
  const Eigen::VectorXd& dx,
  const Eigen::VectorXd& x,
  const std::vector<int>& endmembers_per_phase
) {
  Eigen::ArrayXd max_steps = Eigen::ArrayXd::Constant(x.size(), 100000.0);
  // First two constraints are P & T - use biggest reasonable P-T steps
  max_steps(0) = 20.0e9;
  max_steps(1) = 500.0;
  int j = 2;
  for (int n : endmembers_per_phase) {
    if (x(j) + dx(j) < 0.0) {
      max_steps(j) = std::max(x(j)*0.999, 0.001);
    }
    for (int k = 1; k < n; ++k) {
      max_steps(j + k) = std::max(x(j + k) * 0.99, 0.01);
    }
    j += n;
  }
  double max_lambda = 1.0;
  for (Eigen::Index i = 0; i < dx.size(); ++i) {
    double step = std::abs(dx(i));
    double ratio = (step <= max_steps(i)) ? 1.0 : max_steps(i) / step;
    max_lambda = std::min(max_lambda, ratio);
  }
  return {1.0e-8, max_lambda};
}

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
    prm.parameter_names.push_back("x(" + ph.get_name() + ")");
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
  return prm
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
      bound = assemblage.get_phase<Solution>(i).get_endmember_occupancies();
      n_constraints += bound.cols(); // n_elements
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
  c_matrix(0, 0) = -1.0 // P > 0
  c_matrix(1, 1) = -1.0 // T > 0
  Eigen::Index cidx = 2; // Index of current compositional constraint (row)
  Eigen::Index pidx = 0; // Starting index of current phase
  for (std::size_t i = 0; i < embr_per_phase.size(); ++i) {
    Eigen::Index n = static_cast<Eigen::Index>(embr_per_phase[i]);
    c_matrix(cidx, pidx + 2) = -1.0 // phase prop > 0
    // The first endmember proportion is not a free variable
    // (all endmembers in solution must sum to one)
    // Re-express the constraints without the first endmember
    ++c_idx;
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

Eigen::VectorXd F(
  const Eigen::VectorXd& x,
  Assemblage& assemblage,
  const std::vector<std::unique_ptr<EqualityConstraint>>& equality_constraints,
  const Eigen::VectorXd& reduced_composition_vector,
  const Eigen::MatrixXd& reduced_free_composition_vectors
) {
  // Update assemblage state
  set_composition_and_state_from_parameters(assemblage, x);
  // Retrieve updated endmember amounts
  Eigen::VectorXd new_endmember_amounts = get_endmember_amounts(assemblage);
  // Allocate F
  Eigen::Index n_eqc = static_cast<Eigen::Index>(equality_constraints.size());
  Eigen::VectorXd eqns = Eigen::VectorXd::Zero(
    static_cast<Eigen::Index>(assemblage.get_n_endmembers()) + n_eqc);
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
  Eigen::Index n_reac = static_cas<Eigen::Index>(assemblage.get_n_reactions());
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
  const std::vector<std::unique_ptr<EqualityConstraint>>& equality_constraints,
  const Eigen::MatrixXd& reduced_free_composition_vectors
) {
  Eigen::Index n_eqc = static_cast<Eigen::Index>(equality_constraints.size());
  Eigen::Index n_end = static_cast<Eigen::Index>(assemblage.get_n_endmembers());
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
      n_eq, 2,
      bulk_hessian.rows(), bulk_hessian.cols()
    ) = bulk_hessian;
  }
  if (reduced_free_composition_vectors.rows() > 0) {
    // Swap because transposing
    Eigen::Index nrows = reduced_free_composition_vectors.cols();
    Eigen::Index ncols = reduced_free_composition_vectors.rows();
    jacobian.block(
      jacobian.rows() - nrows, 2 + reaction_hessian.cols(),
      nrows(), ncols()
    ) = -reduced_free_composition_vectors.transpose();
  }
  return jacobian;
}

// TODO: EquilibrateResult or similar to hold output?
EquilibrateResult equilibrate(
  const FormulaMap& composition,
  Assemblage& assemblage,
  const ConstraintList& equality_constraints,
  const std::vector<FreeVectorMap>& free_compositional_vectors = {},
  double tol = 1.0e-3,
  bool store_iterates = false,
  bool store_assemblage = true,
  int max_iterations = 100,
  bool verbose = false
) {
  // Check compositions of solutions set
  // TODO:: could implement a has_composition() for convenience here
  for (std::size_t i = 0; i < static_cast<std::size_t>(assemblage.get_n_phases()); ++i) {
    if (auto& ph = assemblage.get_phase<Solution>(i)) {
      if (!ph->get_molar_fractions().size()) {
        throw std::runtime_error(
          "Set composition for solution " + ph->get_name() + " before running equilibrate!"
        );
      }
    }
  }
  // Check number constraints OK
  std::size_t n_equality_constraints = equality_constraints.size();
  std::size_t n_free_compositional_vectors = free_compositional_vectors.size();
  if (n_equality_constraints != n_free_compositional_vectors + 2) {
    throw std::runtime_error(
      "The number of equality constraints (" + std::to_string(n_equality_constraints)
      + ") must be two more than the number of free_compositional_vectors ("
      + std::to_string(n_free_compositional_vectors) + ")."
    );
  }
  // Check free_compositional_vectors values sum to zero
  for (auto const& vec : free_compositional_vectors) {
    double sum = 0;
    for (const auto& pair : vec) {
      sum += pair.second;
    }
    if (std::abs(sum) > constants::precision::abs_tolerance) {
      throw std::runtime_error(
        "The amounts of each free_compositional_vector must sum to zero");
    }
  }


  // Set default assemblage molar_fractions if none
  if (!assemblage.get_molar_fractions().size()) {
    int n_phases = assemblage->get_n_phases();
    Eigen::ArrayXd f = Eigen::ArrayXd::Constant(n_phases, 1.0 / n_phases);
    assemblage->set_fractions(f);
  }
  // Set n_moles
  double comp_sum = 0
  for (const auto& pair : composition) {
    comp_sum += pair.second;
  }
  double form_sum = 0;
  for (const auto& pair : assemblage.get_formula()) {
    form_sum += pair.second;
  }
  assemblage.set_n_moles(comp_sum / form_sum);

  // Make parameters
  EquilibrationParameters prm = get_equilibration_parameters(
    assemblage, composition, free_compositional_vectors);

  // Set default state if none
  double initial_pressure = 5.0e9;
  double initial_temperature = 1200.0;
  // overwrite if state already set
  if (assemblage.has_state()) {
    initial_pressure = assemblage.get_pressure();
    initial_temperature = assemblage.get_temperature();
  }
  // Overwrite with state from constraints if applicable
  for (const ConstraintGroup& c_group : equality_constraints) {
    const auto& constraint = c_group[0];
    if (auto P_c = dynamic_cast<PressureConstraint*>(constraint.get())) {
      initial_pressure = P_c->value;
    } else if (auto T_c = dynamic_cast<TemperatureConstraint*>(constraint.get())) {
      initial_temperature = T_c->value;
    } else if (auto PTE_c = dynamic_cast<PTEllipseConstraint*>(constraint.get())) {
      Eigen::Vector2d value = PTE_c->scaling;
      initial_pressure = value(0);
      initial_temperature = value(1);
    }
  }
  // Set the state
  assemblage.set_state(initial_pressure, initial_temperature);

  // Get parameter vector (x)
  Eigen::VectorXd parameter_vector = get_parameter_vector(assemblage, n_free_compositional_vectors);

  // Set up solves constraint indices from ConstraintList
  std::vector<std::size_t> grid_shape;
  grid_shape.reserve(n_equality_constraints);
  for (const ConstraintGroup& c_group : equality_constraints) {
    grid_shape.push_back(c_group.size());
  }
  // Store strides for mapping back to grid_index
  std::vector<std::size_t> strides = utils::compute_strides(grid_shape);
  // Set up grid of solver results
  NDArray<DampedNewtonResult> sol_array(grid_shape);

  // Make solver settings
  DampedNewtonSettings solver_settings;
  solver_settings.store_iterates = store_iterates;
  solver_settings.max_iterations = max_iterations;
  solver_settings.tol = tol;
  std::vector<int> embr_per_phase = assemblage.get_endmembers_per_phase();
  solver_settings.lambda_bounds_func =
    [embr_per_phase](const Eigen::VectorXd& dx, const Eigen::VectorXd& x) {
      return lambda_bounds_func(dx, x, embr_per_phase);
    };
  // Make solver object
  DampedNewtonSolver solver(solver_settings);

  // Loop over all the problems
  std::size_t n_problems = sol_array.size();
  for (std::size_t problem_idx = 0; problem_idx < n_problems; ++problem_idx) {
    // Get constraints idx and build current constraints
    std::vector<std::size_t> constraint_indices = utils::map_index(problem_idx, strides);
    ConstraintGroup problem_constraints;
    problem_constraints.reserve(n_equality_constraints);
    for (std::size_t i = 0; i < n_equality_constraints; ++i) {
      // Clone constraints to have independent instances
      // Could also use .get() and use raw pointers (need to modify F, J arg types)
      problem_constraints.emplace_back(
        equality_constraints[i][constraint_indices[i]]->clone());
    }

    // Note - mutates assemblage object (might need to clone before)
    DampedNewtonResult sol = solver.solve(
      parameter_vector,
      [&](const Eigen::VectorXd& x) {
        return F(
          x,
          assemblage,
          problem_constraints,
          prm.reduced_composition_vector,
          prm.reduced_free_composition_vectors);
      },
      [&](const Eigen::VectorXd& x) {
        return J(
          x,
          assemblage,
          problem_constraints,
          prm.reduced_free_composition_vectors);
      },
      {prm.constraint_matrix, prm.constraint_vector}
    );

    double maxres = 0.0;
    if (sol.success && assemblage.get_reaction_affinities().size() > 0) {
      // TODO: assemblage.get_reaction_affinities()
      maxres = assemblage.get_reaction_affinities().cwiseAbs().maxCoeff() + 1.0e-5;
      // TODO: assemblage.set_equilibrium_tolerance()
      assemblage.set_equilibrium_tolerance(maxres);
    }
    // TODO: can't store assemblage as not separate object for each problem in loop
    // Need to implement clone (might work) - see below
    if (store_assemblage) {
      // TODO!
      // Clone with? sol.assemblage = assemblage.clone();
      // Or better, make a new clone at start of each loop, need to be careful with
      // lifetime as clone should return a pointer.
      ;
      // TODO: Maybe set equilibrium tolerance only here and not above?
    }

    // TODO: include iostream
    if (verbose) {
      std::cout << sol.message << std::endl;
    }

    // Store result!
    sol_array(problem_idx) = std::move(sol);

    if (problem_idx < n_problems - 1) {
      // Get next set of constraints
      std::vector<std::size_t> next_c_idx = utils::map_index(
        problem_idx + 1, strides);
      ConstraintGroup next_constraints;
      next_constraints.reserve(n_equality_constraints);
      for (std::size_t i = 0; i < n_equality_constraints; ++i) {
        next_constraints.emplace_back(
          equality_constraints[i][next_c_idx[i]]->clone());
      }
      // Get neighbouring indices in grid
      std::vector<std::size_t> prev_indices;
      for (std::size_t i = 0; i < n_equality_constraints; ++i) {
        if (next_c_idx[i] != 0) {
          std::vector<std::size_t> neighbor_idx = next_c_idx;
          neighbor_idx[i] -= 1;
          prev_indices.push_back(utils::flatten_index(neighbor_idx, strides));
        }
      }
      bool updated_params = false;
      for (std::size_t idx : prev_indices) {
        const DampedNewtonResult& s = sol_array(idx);
        if (s.success && !updated_params) {
          // Next guess based on a Newton step using
          // the old solution vector and Jacobian with new constraints
          Eigen::VectorXd dF = F(
            s.x,
            assemblage,
            next_constraints,
            prm.reduced_composition_vector,
            prm.reduced_free_composition_vectors
          );
          Eigen::PartialPivLU<Eigen::MatrixXd> luJ(s.J);
          Eigen::VectorXd new_parameters = s.x + luJ.solve(-dF);
          // Check constraints are satisfied
          Eigen::VectorXd c = prm.constraint_matrix * new_parameters + prm.constraint_vector;
          if ((c.array() <= 0.0).all()) {
            // Constraints satisfied - accept new guess
            parameters = new_parameters;
          } else {
            // Fallback to last parameter vector
            parameters = s.x;
            if (verbose) {
              std::vector<std::string> exhausted_phases;
              Eigen::ArrayXd phase_amounts = new_parameters(prm.phase_amount_indices);
              for (Eigen::Index i = 0; i < phase_amounts.size(); ++i) {
                if (phase_amounts(i) < 0.0) {
                  // TODO: update to Eigen::Index
                  exhausted_phases.push_back(
                    assemblage.get_phase(i)->get_name()
                  );
                }
              }
              if (!exhausted_phases.empty()) {
                std::string ex_p_message =
                  "A phase might be exhausted before the next step: [";
                for (const std::string& name : exhausted_phases) {
                  ex_p_message += name + " ";
                }
                ex_p_message += "]";
                std::cout << ex_p_message << std::endl;
              }
            } // exhausted phase check if verbose
          } // parameter fallback
          updated_params = true;
        } // update param check
      } // neighbouring constraints loop
    } // not last problem check
  } // loop over all problems

  EquilibrateResult result;
  result.sol_array = std::move(sol_array);
  result.prm = std::move(prm);
  return result;
}

// TODO: if we want to store assemblages like in the python then they need to be copyable:
//       maybe possible with clones?
//       Material:
//         virtual std::shared_ptr<Material> clone() const = 0;
//       Mineral:
//         std::shared_ptr<Material> clone() const override {
//           return std::make_shared<Mineral>(*this);
//         }
//       SolutionModel:
//         virtual std::shared_ptr<SolutionModel> clone() const = 0;
//       IdealSolution etc.:
//         std::shared_ptr<SolutionModel> clone() const override {
//           return std::make_shared<IdealSolutionModel>(*this);
//         }
//       Solution: (implement to copy solution models)
//         std::shared_ptr<Material> clone() const override {
//           auto copy = std::make_shared<Solution>(*this);
//           if (solution_model) {
//             copy->solution_model = solution_model->clone();
//           }
//           return copy;
//          }
//       Assemblage (implement to copy phases)
//         std::shared_ptr<Material> clone() const override {
//           auto copy = std::make_shared<Assemblage>(*this);
//           copy->phases.clear();
//           for (const auto& phase : phases) {
//             copy->phases.push_back(phase->clone());
//           }
//           return copy;
//         }
//       Then: std::shared_ptr<Assemblage> copy = std::dynamic_pointer_cast<Assemblage>(assemblage.clone());
