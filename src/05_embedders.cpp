#include <RcppEigen.h>

#include <algorithm>
#include <cmath>
#include <vector>

//[[Rcpp::depends(RcppEigen)]]

inline double stable_sigmoid(double value) {
  if (value >= 0.0) {
    const double exponential = std::exp(-value);
    return 1.0 / (1.0 + exponential);
  }
  const double exponential = std::exp(value);
  return exponential / (1.0 + exponential);
}

inline double stable_softplus(double value) {
  return std::max(value, 0.0) + std::log1p(std::exp(-std::abs(value)));
}

double lsm_log_likelihood_cpp(const Eigen::MatrixXd& A,
                              const Eigen::MatrixXd& theta) {
  const Eigen::Index n = A.rows();
  double objective = 0.0;
  for (Eigen::Index row = 0; row < n; ++row) {
    for (Eigen::Index column = row + 1; column < n; ++column) {
      objective += A(row, column) * theta(row, column) -
        stable_softplus(theta(row, column));
    }
  }
  return objective;
}

void validate_lsm_pgd_cpp_inputs(const Eigen::MatrixXd& A,
                                 const Eigen::MatrixXd& Z,
                                 const Eigen::VectorXd& alpha,
                                 double step_size_Z,
                                 double step_size_alpha,
                                 int niter) {
  const Eigen::Index n = A.rows();
  if (A.cols() != n) Rcpp::stop("A must be square");
  if (Z.rows() != n || Z.cols() < 1) {
    Rcpp::stop("Z must have nrow(A) rows and at least one column");
  }
  if (alpha.size() != n) Rcpp::stop("alpha must have length nrow(A)");
  if (!A.allFinite() || !Z.allFinite() || !alpha.allFinite()) {
    Rcpp::stop("A, Z, and alpha must contain only finite values");
  }
  if (!std::isfinite(step_size_Z) || step_size_Z <= 0.0 ||
      !std::isfinite(step_size_alpha) || step_size_alpha <= 0.0) {
    Rcpp::stop("step sizes must be finite and positive");
  }
  if (niter < 1) Rcpp::stop("niter must be positive");
  const double symmetry_tolerance = 1e-10;
  for (Eigen::Index row = 0; row < n; ++row) {
    if (A(row, row) != 0.0) Rcpp::stop("A must have a zero diagonal");
    for (Eigen::Index column = 0; column < n; ++column) {
      const double value = A(row, column);
      if (value != 0.0 && value != 1.0) {
        Rcpp::stop("A must contain only binary values 0 and 1");
      }
      if (std::abs(value - A(column, row)) > symmetry_tolerance) {
        Rcpp::stop("A must be symmetric");
      }
    }
  }
}

// Projected gradient ascent for an undirected logistic latent-space model.
// Diagonal dyads are excluded from both the objective and gradient.
// [[Rcpp::export]]
Rcpp::List lsm_pgd_cpp(const Eigen::Map<Eigen::MatrixXd>& A_map,
                       Eigen::MatrixXd Z,
                       Eigen::VectorXd alpha,
                       double step_size_Z,
                       double step_size_alpha,
                       int niter,
                       bool trace = false) {
  const Eigen::MatrixXd A = A_map;
  validate_lsm_pgd_cpp_inputs(
    A,
    Z,
    alpha,
    step_size_Z,
    step_size_alpha,
    niter
  );
  const Eigen::Index n = A.rows();
  const Eigen::VectorXd ones = Eigen::VectorXd::Ones(n);
  std::vector<double> objective;
  objective.reserve(static_cast<std::size_t>(niter + 1));

  for (int iteration = 0; iteration < niter; ++iteration) {
    Eigen::MatrixXd theta = alpha * ones.transpose();
    theta.noalias() += ones * alpha.transpose();
    theta.noalias() += Z * Z.transpose();
    objective.push_back(lsm_log_likelihood_cpp(A, theta));
    if (trace) {
      Rcpp::Rcout << "Iteration " << iteration + 1
                  << " objective: " << objective.back() << "\n";
    }

    Eigen::MatrixXd P_hat = theta.unaryExpr(
      [](double value) { return stable_sigmoid(value); }
    );
    Eigen::MatrixXd residual = A - P_hat;
    residual.diagonal().setZero();
    const Eigen::MatrixXd update_Z = residual * Z;
    Z += step_size_Z * update_Z;
    alpha.noalias() += step_size_alpha * residual.rowwise().sum();
    Z.rowwise() -= Z.colwise().mean();

    if ((iteration + 1) % 10 == 0) Rcpp::checkUserInterrupt();
  }

  Eigen::MatrixXd theta = alpha * ones.transpose();
  theta.noalias() += ones * alpha.transpose();
  theta.noalias() += Z * Z.transpose();
  objective.push_back(lsm_log_likelihood_cpp(A, theta));
  Eigen::MatrixXd P_hat = theta.unaryExpr(
    [](double value) { return stable_sigmoid(value); }
  );
  P_hat.diagonal().setZero();

  return Rcpp::List::create(
    Rcpp::Named("Z_hat") = Z,
    Rcpp::Named("alpha_hat") = alpha,
    Rcpp::Named("P_hat") = P_hat,
    Rcpp::Named("objective") = objective
  );
}
