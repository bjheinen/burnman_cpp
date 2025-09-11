/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include "burnman/tools/equilibration/equilibrate_impl.hpp"
#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include "burnman/utils/types/ndarray.hpp"
#include "burnman/core/solution.hpp"
#include "burnman/optim/roots/damped_newton_types.hpp"
#include "burnman/optim/roots/damped_newton_solver.hpp"
#include "burnman/tools/equilibration/equality_constraint_variants.hpp"
#include "burnman/tools/equilibration/equilibrate_lambda_bounds.hpp"
#include "burnman/tools/equilibration/equilibrate_utils.hpp"
#include "burnman/tools/equilibration/equilibrate_objective.hpp"

EquilibrateResult equilibrate(
  const FormulaMap& composition,
  Assemblage& assemblage,
  const ConstraintList& equality_constraints,
  const std::vector<FreeVectorMap>& free_compositional_vectors,
  double tol,
  bool store_iterates,
  bool store_assemblage ,
  int max_iterations,
  bool verbose
) {
  // Check compositions of solutions set
  // TODO:: could implement a has_composition() for convenience here
  for (std::size_t i = 0; i < static_cast<std::size_t>(assemblage.get_n_phases()); ++i) {
    if (const auto& ph = assemblage.get_phase<Solution>(i)) {
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
    int n_phases = assemblage.get_n_phases();
    Eigen::ArrayXd f = Eigen::ArrayXd::Constant(n_phases, 1.0 / n_phases);
    assemblage.set_fractions(f);
  }
  // Set n_moles
  double comp_sum = 0;
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
  NDArray<optim::roots::DampedNewtonResult> sol_array(grid_shape);

  // Make solver settings
  optim::roots::DampedNewtonSettings solver_settings;
  solver_settings.store_iterates = store_iterates;
  solver_settings.max_iterations = max_iterations;
  solver_settings.tol = tol;
  std::vector<int> embr_per_phase = assemblage.get_endmembers_per_phase();
  solver_settings.lambda_bounds_func =
    [embr_per_phase](const Eigen::VectorXd& dx, const Eigen::VectorXd& x) {
      return lambda_bounds_func(dx, x, embr_per_phase);
    };
  // Make solver object
  optim::roots::DampedNewtonSolver solver(solver_settings);

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
    optim::roots::DampedNewtonResult sol = solver.solve(
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
        const optim::roots::DampedNewtonResult& s = sol_array(idx);
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
            parameter_vector = new_parameters;
          } else {
            // Fallback to last parameter vector
            parameter_vector = s.x;
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
