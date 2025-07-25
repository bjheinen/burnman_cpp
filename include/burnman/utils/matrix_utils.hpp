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

#include <algorithm>
#include <vector>
#include <Eigen/Dense>

namespace utils {

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
  inline Eigen::MatrixXd jagged2square(const std::vector<std::vector<double>>& v, int n) {
    Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(n, n);
    for (std::size_t i = 0; i < v.size(); i++) {
      std::ptrdiff_t col = n - 1;
      for (std::ptrdiff_t j = static_cast<std::ptrdiff_t>(v[i].size()) - 1; j >= 0; --j) {
        mat(i, col) = v[i][j];
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
    int n
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
   * @brief Retrieves indices of linearly independent rows
   *
   * Uses QR decomposition to get independent columns. A transpose
   * is applied at the start so that independent rows are retrieved.
   *
   * @param mat Matrix
   * @returns Sorted list of indices of linearly independent rows
   */
  inline std::vector<int> get_independent_row_indices(
    const Eigen::MatrixXd& mat
  ) {
    // Transpose to get rows instead of cols
    Eigen::MatrixXd A = mat.transpose();
    // QR decomp with column pivoting
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(A);
    // Use rank as num independent cols
    Eigen::Index rank = qr.rank();
    // colsPermutation() retrieves P
    // head(rank) takes needed block
    Eigen::VectorXi indices_vector = qr.colsPermutation().indices().head(rank);
    std::vector<int> indices(indices_vector.begin(), indices_vector.begin() + rank);
    std::sort(indices.begin(), indices.end());
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
    Eigen::Index n = basis.rows();
    Eigen::Index m = basis.cols();
    if (n >= m) {
      return basis;
    }
    // Flip basis
    Eigen::MatrixXd basis_flipped = basis.rowwise().reverse();
    // FullPivLU decomposition to get RREF
    Eigen::FullPivLU<Eigen::MatrixXd> lu_basis(basis_flipped);
    // Get the rank of the matrix
    Eigen::Index rank = lu_basis.rank();
    Eigen::MatrixXd U = lu_basis.matrixLU().triangularView<Eigen::Upper>();
    // Get the indices of the pivot columns for each row -
    // first non-zero element in each row.
    std::vector<Eigen::Index> pivot_columns;
    for (Eigen::Index i = 0; i < rank; ++i) {
      for (Eigen::Index j = i; j < m; ++j) {
        if (U(i, j) != 0) {
          pivot_columns.push_back(m - 1 - j);
          break;
        }
      }
    }
    // Get excluded indices
    std::vector<bool> excluded(static_cast<std::size_t>(m), false);
    for (Eigen::Index val : pivot_columns) {
      if (val >= 0 && val < m)
        excluded[static_cast<std::size_t>(val)] = true;
    }
    // Get indices of identity matrix to use
    std::vector<Eigen::Index> basis_indices;
    basis_indices.reserve(static_cast<std::size_t>(m) - pivot_columns.size());
    for (Eigen::Index i = 0; i < m; ++i) {
      if (!excluded[static_cast<std::size_t>(i)])
        basis_indices.push_back(i);
    }
    // Get identity matrix and concatenate with basis
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(m, m);
    Eigen::MatrixXd full_basis(n + basis_indices.size(), m);
    full_basis <<
      basis,
      I(basis_indices, Eigen::all);
    return full_basis;
  }

}

#endif // BURNMAN_UTILS_MATRIX_UTILS_INCLUDED
