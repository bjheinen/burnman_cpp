/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_UTILS_MATRIX_UTILS_INCLUDED
#define BURNMAN_UTILS_MATRIX_UTILS_INCLUDED

#include <cmath>
#include <algorithm>
#include <cstddef>
#include <vector>
#include <Eigen/Dense>

#include <iostream>

namespace burnman {
namespace utils {

  /**
  * @struct RREFResult
  * @brief Holds the result of a Reduced Row Echelon Form (RREF) computation.
  *
  * Contains the RREF matrix, the indices of pivot columns, and the rank of the matrix.
  */
  struct RREFResult {
    Eigen::MatrixXd rref_matrix;
    Eigen::VectorXi pivot_columns;
    int rank;
  };

  /**
  * @brief Computes the Reduced Row Echelon Form (RREF) of a given matrix.
  *
  * Computes RREF by manual gaussian elimination. Keeps track of indices of pivot columns.
  *
  * @param M The input matrix..
  * @param tol A tolerance value to determine when an element is zero. Default is 1.0e-12.
  *
  * @return RREFResult Struct containing the RREF matrix, pivot column indices, and matrix rank.
  */
  inline RREFResult compute_rref(Eigen::MatrixXd M, double tol = 1.0e-12) {
    // Get number of rows/cols
    const Eigen::Index rows = M.rows();
    const Eigen::Index cols = M.cols();
    // Keep track of pivot row index
    Eigen::Index pivot_row = 0;
    // Keep track of pivot columns
    std::vector<int> pivot_columns_list;
    // Loop through and construct RREF
    for (Eigen::Index col = 0; col < cols && pivot_row < rows; ++col) {
      // Find pivot in the column - the max value
      Eigen::Index pivot;
      M.col(col).tail(rows - pivot_row).cwiseAbs().maxCoeff(&pivot);
      pivot += pivot_row;
      // Skip if no pivot in column
      if (std::abs(M(pivot, col)) < tol) {
        continue;
      }
      // Swap pivot row into position
      M.row(pivot_row).swap(M.row(pivot));
      // Normalise pivot row
      M.row(pivot_row) /= M(pivot_row, col);
      // Eliminate other rows in the column
      for (Eigen::Index i = 0; i < rows; ++i) {
        if (i != pivot_row) {
          M.row(i) -= M(i, col) * M.row(pivot_row);
        }
      }
      // Store pivot column index
      pivot_columns_list.push_back(static_cast<int>(col));
      ++pivot_row;
    }
    // Cleanup - zero small values
    M = (M.cwiseAbs().array() < tol).select(0.0, M);
    // Convert pivot columns to VectorXi
    Eigen::VectorXi pivot_columns = Eigen::Map<Eigen::VectorXi>(
      pivot_columns_list.data(), static_cast<Eigen::Index>(pivot_columns_list.size())
    );
    return {M, pivot_columns, static_cast<int>(pivot_columns.size())};
  }

  /**
  * @brief Computes the canonical nullspace basis of a given matrix.
  *
  * Uses the Reduced Row Echelon Form (RREF) of the matrix to find the
  * canonical basis for its nullspace. Basis vectors are returned as rows
  * in the output matrix.
  *
  * @param M Input matrix.
  * @param tol A tolerance value to determine when an element is considered zero. Default is 1.0e-12.
  *
  * @return Nullspace basis matrix.
  */
  inline Eigen::MatrixXd nullspace(const Eigen::MatrixXd& M, double tol = 1.0e-12) {
    // Get RREF
    RREFResult rref = compute_rref(M, tol);
    const Eigen::MatrixXd& R = rref.rref_matrix;
    const Eigen::Index cols = R.cols();
    const Eigen::VectorXi& pivots = rref.pivot_columns;
    const Eigen::Index num_pivots = pivots.size();
    // Make a boolean array for pivot locations
    Eigen::Array<bool, Eigen::Dynamic, 1> is_pivot = Eigen::Array<bool, Eigen::Dynamic, 1>::Constant(cols, false);
    is_pivot(pivots) = true;
    // Pre-allocate the nullspace matrix
    const Eigen::Index nullity = cols - num_pivots;
    Eigen::MatrixXd basis(cols, nullity);
    // Loop through and fill
    Eigen::Index basis_col = 0;
    for (Eigen::Index free_index = 0; free_index < cols; ++free_index) {
      // Skip if pivot (not free)
      if (is_pivot(free_index)) {
        continue;
      }
      Eigen::VectorXd vec = Eigen::VectorXd::Zero(cols);
      vec(free_index) = 1.0;
      for (Eigen::Index i = 0; i < num_pivots; ++i) {
        vec(pivots(i)) = -R(i, free_index);
      }
      // Store column and increment
      basis.col(basis_col++) = vec;
    }
    // Transpose to match expected form
    return basis.transpose();
  }

  /**
   * @brief Convert jagged upper triangle to Eigen::Matrix
   *
   * e.g.:
   *   std::vector<std::vector<double>> v = {
   *     {0.0, 24.74e3, 26.0e3, 24.3e3},
   *     {24.74e3, 0.0, 0.0e3},
   *     {60.53136e3, 0.0},
   *     {10.0e3}
   *   };
   *   int n = 5;
   *   jagged2square(x, n) -->
   *     0, 0, 24740, 26000, 24300,
   *     0, 0, 24740, 0, 0,
   *     0, 0, 0, 60531.36, 0,
   *     0, 0, 0, 0, 10000,
   *     0, 0, 0, 0, 0;
   *
   * @param v jagged array to convert.
   * @param n size of square matrix.
   * @return Square matrix.
   */
  inline Eigen::MatrixXd jagged2square(const std::vector<std::vector<double>>& v, Eigen::Index n) {
    Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(n, n);
    for (std::size_t i = 0; i < v.size(); i++) {
      std::ptrdiff_t col = static_cast<std::ptrdiff_t>(n) - 1;
      for (std::ptrdiff_t j = static_cast<std::ptrdiff_t>(v[i].size()) - 1; j >= 0; --j) {
        mat(static_cast<Eigen::Index>(i), col) = v[i][j];
        --col;
      }
    }
    // Ensure matrix is upper triangle only
    mat.triangularView<Eigen::Lower>().setZero();
    return mat;
  }

  inline Eigen::MatrixXd populate_interaction_matrix(
    const Eigen::MatrixXd& interaction,
    const Eigen::ArrayXd& alphas,
    Eigen::Index n
  ) {
    Eigen::MatrixXd alphas_outer_sum = alphas.replicate(1, n) + alphas.transpose().replicate(n, 1);
    // Compute 2.0 / alpha_sum element-wise
    Eigen::MatrixXd weights = 2.0 / alphas_outer_sum.array();
    // Set lower triangles to zero
    weights.triangularView<Eigen::Lower>().setZero();
    // Element-wise product (only upper triangle is non-zero)
    Eigen::MatrixXd interaction_matrix = weights.cwiseProduct(interaction);
    return interaction_matrix;
  }

  /**
   * @brief Retrieves indices of linearly independent columns
   *
   * Helper function to return the list of linearly indepedent columns
   * of a matrix (i.e. the pivot column indices). This functions is
   * left here for convenience. It currently maps the Eigen::VectorXi
   * to a std::vector.
   *
   * @param mat Matrix
   * @returns Sorted list of indices of linearly independent cols
   */
  inline std::vector<Eigen::Index> get_independent_col_indices(
    const Eigen::MatrixXd& mat
  ) {
    // TODO: Use Eigen::ArrayXi for indices?
    // Use manual RREF for pivots
    // RREF currently guarantees increasing index order. If
    // the implementation is changed we should sort here for safety.
    RREFResult rref = compute_rref(mat);
    const Eigen::VectorXi& pivots = rref.pivot_columns;
    std::vector<Eigen::Index> indices(pivots.begin(), pivots.end());
    //std::sort(indices.begin(), indices.end());
    return indices;
  }

  // TODO: Make move to separate cpp files...
  /**
   * @brief Creates a full basis by filling remaining rows with
   * selected rows of the identity matrix that have indices not in
   * the column pivot list of the basis RREF.
   *
   * @param basis Partial basis
   * @returns full_basis Complete basis.
   */
  inline Eigen::MatrixXd complete_basis(
    const Eigen::MatrixXd& basis
  ) {
    const Eigen::Index n = basis.rows();
    const Eigen::Index m = basis.cols();
    if (n >= m) {
      return basis;
    }
    // Flip basis
    Eigen::MatrixXd basis_flipped = basis.rowwise().reverse();
    // Get pivot columns from manual RREF
    RREFResult rref = compute_rref(basis_flipped);
    const Eigen::VectorXi& pivots_flipped = rref.pivot_columns;
    // Flip pivot column indices back
    Eigen::VectorXi pivots = (static_cast<int>(m) - 1) - pivots_flipped.array();
    // Get excluded indices
    Eigen::Array<bool, Eigen::Dynamic, 1> excluded = Eigen::Array<bool, Eigen::Dynamic, 1>::Constant(m, false);
    excluded(pivots) = true;
    std::vector<Eigen::Index> basis_indices;
    for (Eigen::Index i = 0; i < m; ++i) {
      if (!excluded(i)) {
        basis_indices.push_back(i);
      }
    }
    // Get identity matrix and concatenate with basis
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(m, m);
    Eigen::MatrixXd full_basis(n + static_cast<Eigen::Index>(basis_indices.size()), m);
    full_basis <<
      basis,
      I(basis_indices, Eigen::all);
    return full_basis;
  }

} // namespace utils
} // namespace burnman

#endif // BURNMAN_UTILS_MATRIX_UTILS_INCLUDED
