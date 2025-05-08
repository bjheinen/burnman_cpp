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
   * jagged2mat(x) --> 
   *   0, 0, 24740, 26000, 24300,
   *   0, 0, 24740, 0, 0,
   *   0, 0, 0, 60531.36, 0,
   *   0, 0, 0, 0, 10000,
   *   0, 0, 0, 0, 0;
   * 
   * @param v jagged array to convert.
   * @param n size of square matrix.
   * @return Square matrix.
   */
  Eigen::MatrixXd jagged2square(std::vector<std::vector<double>>& v, int n) {
    Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(n, n);
    for (int i = 0; i < v.size(); i++) {
      int col = n - 1;
      for (int j = v[i].size() - 1; j >= 0; --j) {
        mat(i, col) = v[i][j];
        --col;
      }
    }
    // Ensure matrix is upper triangle only
    mat.triangularView<Eigen::Lower>().setZero();
  }

  Eigen::MatrixXd populate_interaction_matrix(
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

}

#endif // BURNMAN_UTILS_MATRIX_UTILS_INCLUDED
