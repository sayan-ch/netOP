#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <queue>
#include <utility>
#include <vector>

struct Edge {
  int node;
  double weight;
};

using adjacency_list = std::vector<std::vector<Edge>>;

inline bool is_nonzero(double value) {
  return !R_IsNA(value) && !R_IsNaN(value) && value != 0.0;
}

inline double undirected_weight(double outgoing, double incoming) {
  if (!is_nonzero(outgoing)) return incoming;
  if (!is_nonzero(incoming)) return outgoing;
  return std::min(outgoing, incoming);
}

adjacency_list build_dense_adjacency(const Rcpp::NumericMatrix& A,
                                     bool directed,
                                     bool weighted,
                                     bool include_self_loops) {
  const int n = A.nrow();
  if (A.ncol() != n) Rcpp::stop("A must be square");
  for (R_xlen_t index = 0; index < A.size(); ++index) {
    if (!std::isfinite(A[index])) {
      Rcpp::stop("A must contain only finite values");
    }
    if (weighted && A[index] < 0.0) {
      Rcpp::stop("weighted shortest paths require nonnegative edge weights");
    }
  }
  adjacency_list neighbors(static_cast<std::size_t>(n));
  for (int node = 0; node < n; ++node) {
    neighbors[node].reserve(static_cast<std::size_t>(n));
    for (int neighbor = 0; neighbor < n; ++neighbor) {
      if (!include_self_loops && node == neighbor) continue;
      const double outgoing = A(node, neighbor);
      const double incoming = A(neighbor, node);
      const bool connected = is_nonzero(outgoing) ||
        (!directed && is_nonzero(incoming));
      if (!connected) continue;
      const double edge_weight = weighted
        ? (directed ? outgoing : undirected_weight(outgoing, incoming))
        : 1.0;
      neighbors[node].push_back({neighbor, edge_weight});
    }
  }
  return neighbors;
}

adjacency_list build_sparse_adjacency(SEXP A_sexp,
                                      bool directed,
                                      bool weighted,
                                      bool include_self_loops) {
  Rcpp::S4 A(A_sexp);
  const Rcpp::IntegerVector dimensions = A.slot("Dim");
  const int n = dimensions[0];
  if (dimensions[1] != n) Rcpp::stop("A must be square");
  const Rcpp::IntegerVector column_pointers = A.slot("p");
  const Rcpp::IntegerVector row_indices = A.slot("i");
  const Rcpp::NumericVector values = A.slot("x");
  adjacency_list neighbors(static_cast<std::size_t>(n));

  for (R_xlen_t index = 0; index < values.size(); ++index) {
    if (!std::isfinite(values[index])) {
      Rcpp::stop("A must contain only finite values");
    }
    if (weighted && values[index] < 0.0) {
      Rcpp::stop("weighted shortest paths require nonnegative edge weights");
    }
  }

  for (int column = 0; column < n; ++column) {
    for (int index = column_pointers[column];
         index < column_pointers[column + 1]; ++index) {
      const double value = values[index];
      if (!is_nonzero(value)) continue;
      const int row = row_indices[index];
      if (!include_self_loops && row == column) continue;
      const double edge_weight = weighted ? value : 1.0;
      neighbors[row].push_back({column, edge_weight});
      if (!directed && row != column) {
        neighbors[column].push_back({row, edge_weight});
      }
    }
  }

  for (std::vector<Edge>& node_neighbors : neighbors) {
    std::sort(
      node_neighbors.begin(), node_neighbors.end(),
      [](const Edge& left, const Edge& right) {
        return left.node < right.node;
      }
    );
    std::vector<Edge> unique_neighbors;
    unique_neighbors.reserve(node_neighbors.size());
    for (const Edge& edge : node_neighbors) {
      if (unique_neighbors.empty() ||
          unique_neighbors.back().node != edge.node) {
        unique_neighbors.push_back(edge);
      } else {
        unique_neighbors.back().weight = std::min(
          unique_neighbors.back().weight, edge.weight
        );
      }
    }
    node_neighbors.swap(unique_neighbors);
  }
  return neighbors;
}

Rcpp::List connected_components_from_adjacency(
    const adjacency_list& neighbors) {
  const int n = static_cast<int>(neighbors.size());
  std::vector<unsigned char> visited(static_cast<std::size_t>(n), 0);
  std::vector<std::vector<int>> components;
  std::vector<int> queue(static_cast<std::size_t>(n));

  for (int start_node = 0; start_node < n; ++start_node) {
    if (visited[start_node]) continue;
    std::size_t head = 0;
    std::size_t tail = 0;
    queue[tail++] = start_node;
    visited[start_node] = 1;
    std::vector<int> component;
    while (head < tail) {
      const int node = queue[head++];
      component.push_back(node + 1);
      for (const Edge& edge : neighbors[node]) {
        if (!visited[edge.node]) {
          visited[edge.node] = 1;
          queue[tail++] = edge.node;
        }
      }
    }
    components.push_back(component);
  }

  Rcpp::List result(components.size());
  for (std::size_t index = 0; index < components.size(); ++index) {
    result[index] = Rcpp::wrap(components[index]);
  }
  return result;
}

Rcpp::NumericMatrix unweighted_distances_from_adjacency(
    const adjacency_list& neighbors) {
  const int n = static_cast<int>(neighbors.size());
  Rcpp::NumericMatrix distances = Rcpp::no_init_matrix(n, n);
  std::fill(distances.begin(), distances.end(),
            std::numeric_limits<double>::infinity());
  std::vector<int> queue(static_cast<std::size_t>(n));

  for (int source = 0; source < n; ++source) {
    std::size_t head = 0;
    std::size_t tail = 0;
    queue[tail++] = source;
    distances(source, source) = 0.0;
    while (head < tail) {
      const int node = queue[head++];
      const double next_distance = distances(source, node) + 1.0;
      for (const Edge& edge : neighbors[node]) {
        if (!std::isfinite(distances(source, edge.node))) {
          distances(source, edge.node) = next_distance;
          queue[tail++] = edge.node;
        }
      }
    }
  }
  return distances;
}

Rcpp::NumericMatrix weighted_distances_from_adjacency(
    const adjacency_list& neighbors) {
  const int n = static_cast<int>(neighbors.size());
  Rcpp::NumericMatrix distances = Rcpp::no_init_matrix(n, n);
  std::fill(distances.begin(), distances.end(),
            std::numeric_limits<double>::infinity());
  using queue_entry = std::pair<double, int>;

  for (int source = 0; source < n; ++source) {
    std::priority_queue<queue_entry, std::vector<queue_entry>,
                        std::greater<queue_entry>> queue;
    distances(source, source) = 0.0;
    queue.push({0.0, source});
    while (!queue.empty()) {
      const double node_distance = queue.top().first;
      const int node = queue.top().second;
      queue.pop();
      if (node_distance > distances(source, node)) continue;
      for (const Edge& edge : neighbors[node]) {
        const double candidate_distance = node_distance + edge.weight;
        if (candidate_distance < distances(source, edge.node)) {
          distances(source, edge.node) = candidate_distance;
          queue.push({candidate_distance, edge.node});
        }
      }
    }
  }
  return distances;
}

Rcpp::NumericMatrix shortest_path_distances_from_adjacency(
    const adjacency_list& neighbors, bool weighted) {
  return weighted
    ? weighted_distances_from_adjacency(neighbors)
    : unweighted_distances_from_adjacency(neighbors);
}

// [[Rcpp::export]]
Rcpp::List connected_components_dense_cpp(
    Rcpp::NumericMatrix A, bool include_self_loops = false) {
  return connected_components_from_adjacency(
    build_dense_adjacency(A, false, false, include_self_loops)
  );
}

// [[Rcpp::export]]
Rcpp::List connected_components_sparse_cpp(
    SEXP A, bool include_self_loops = false) {
  return connected_components_from_adjacency(
    build_sparse_adjacency(A, false, false, include_self_loops)
  );
}

// [[Rcpp::export]]
Rcpp::NumericMatrix shortest_path_distances_dense_cpp(
    Rcpp::NumericMatrix A,
    bool directed = false,
    bool weighted = false,
    bool include_self_loops = false) {
  return shortest_path_distances_from_adjacency(
    build_dense_adjacency(A, directed, weighted, include_self_loops), weighted
  );
}

// [[Rcpp::export]]
Rcpp::NumericMatrix shortest_path_distances_sparse_cpp(
    SEXP A,
    bool directed = false,
    bool weighted = false,
    bool include_self_loops = false) {
  return shortest_path_distances_from_adjacency(
    build_sparse_adjacency(A, directed, weighted, include_self_loops), weighted
  );
}
