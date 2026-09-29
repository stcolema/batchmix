// minVIGreedy.cpp
// =============================================================================
// Greedy item-relocation refinement of a point-estimate clustering under the
// Jensen lower bound to the posterior expected Variation of Information
// (Wade & Ghahramani, 2018, Bayesian Analysis 13(2), Sec. 3.2; the "greedy"
// method of mcclust.ext::minVI). Local, dependency-free replacement for the
// stochastic search in salso (Dahl, Johnson & Müller, 2022), so that no Rust
// toolchain is needed to install the package.
# include <RcppArmadillo.h>
# include <vector>

using namespace Rcpp ;
using namespace arma ;

//' @title Greedy refinement of a clustering under the VI lower bound
//' @description Starting from \code{init}, repeatedly moves single items to
//' the cluster (or to a new singleton) that most reduces
//' \deqn{f(c) = \frac{1}{n} \sum_i \left[ \log_2 n_{c_i} + \log_2 \sum_j
//' \psi_{ij} - 2 \log_2 \sum_{j : c_j = c_i} \psi_{ij} \right],}
//' the lower bound on the posterior expected Variation of Information
//' obtained from Jensen's inequality, until no move improves it (or
//' \code{max_passes} sweeps have run). Each candidate move costs time linear
//' in the sizes of the two clusters involved, so a sweep is \eqn{O(n^2 G)}
//' for \eqn{G} clusters. Memory is \eqn{O(nG)} beyond the posterior
//' similarity matrix itself.
//' @param psm N x N posterior similarity matrix.
//' @param init N-vector of initial cluster labels (any integer labelling).
//' @param max_passes Maximum number of full sweeps over the items.
//' @return A list with \code{labels} (1-indexed, consecutive) and
//' \code{value}, the lower-bound loss of the returned clustering.
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::List minVIGreedyRefine(arma::mat psm, arma::uvec init, arma::uword max_passes = 50) {

  const uword n = psm.n_rows;
  if(psm.n_cols != n || init.n_elem != n) {
    Rcpp::stop("psm must be square and init must have one label per row of psm.");
  }

  // Relabel init to 0..G-1.
  std::vector<uword> cl(n);
  {
    uvec u = unique(init);
    for(uword i = 0; i < n; i++) {
      cl[i] = (uword) arma::as_scalar(find(u == init(i), 1));
    }
  }
  uword G = 0;
  for(uword i = 0; i < n; i++) G = std::max(G, cl[i] + 1);

  // M[g](j) = sum_{l in g} psm(j, l); size[g] = |g|.
  std::vector<vec> M(G, zeros<vec>(n));
  std::vector<double> size(G, 0.0);
  for(uword i = 0; i < n; i++) {
    M[cl[i]] += psm.col(i);
    size[cl[i]] += 1.0;
  }

  auto members_of = [&](uword g) {
    std::vector<uword> out;
    for(uword j = 0; j < n; j++) if(cl[j] == g) out.push_back(j);
    return out;
  };

  const double inv_ln2 = 1.0 / std::log(2.0);
  auto lg2 = [&](double x) { return std::log(x) * inv_ln2; };

  for(uword pass = 0; pass < max_passes; pass++) {
    bool moved = false;

    for(uword i = 0; i < n; i++) {
      const uword a = cl[i];

      // Remove i from its cluster (tentatively).
      M[a] -= psm.col(i);
      size[a] -= 1.0;
      cl[i] = n; // sentinel: belongs to no cluster while being placed

      // Cost, relative to "i placed nowhere", of adding i to each candidate.
      uword best_g = G; // G denotes a new singleton
      double best_delta = 0.0; // singleton: log2(1) - 2 log2(1) = 0

      for(uword g = 0; g < G; g++) {
        if(size[g] < 0.5) continue;
        std::vector<uword> mem = members_of(g);
        double delta = 0.0;
        for(uword j : mem) {
          delta += lg2(size[g] + 1.0) - 2.0 * lg2(M[g](j) + psm(j, i))
                 - lg2(size[g]) + 2.0 * lg2(M[g](j));
        }
        delta += lg2(size[g] + 1.0) - 2.0 * lg2(M[g](i) + psm(i, i));
        if(delta < best_delta - 1e-12) {
          best_delta = delta;
          best_g = g;
        }
      }

      // Prefer staying put on ties: compare with the original placement.
      // Returning to a (which may now be empty, in which case it equals the
      // singleton value of 0).
      uword target = best_g;
      double delta_a = 0.0;
      if(size[a] >= 0.5) {
        std::vector<uword> mem = members_of(a);
        for(uword j : mem) {
          delta_a += lg2(size[a] + 1.0) - 2.0 * lg2(M[a](j) + psm(j, i))
                   - lg2(size[a]) + 2.0 * lg2(M[a](j));
        }
        delta_a += lg2(size[a] + 1.0) - 2.0 * lg2(M[a](i) + psm(i, i));
      }
      if(delta_a <= best_delta + 1e-12) {
        target = a;
      }

      if(target == G) {
        // New singleton, reusing an emptied column when there is one.
        uword slot = G;
        for(uword g = 0; g < G; g++) {
          if(size[g] < 0.5) { slot = g; break; }
        }
        if(slot == G) {
          M.push_back(zeros<vec>(n));
          size.push_back(0.0);
          G++;
        }
        target = slot;
      }

      M[target] += psm.col(i);
      size[target] += 1.0;
      cl[i] = target;
      if(target != a) moved = true;
    }

    if(!moved) break;
  }

  // Relabel to consecutive 1-indexed labels and evaluate the loss.
  std::vector<uword> map(G + 1, n);
  uword next = 0;
  uvec labels_out(n);
  for(uword i = 0; i < n; i++) {
    if(map[cl[i]] == n) map[cl[i]] = next++;
    labels_out(i) = map[cl[i]] + 1;
  }

  vec s = sum(psm, 1);
  double loss = 0.0;
  for(uword i = 0; i < n; i++) {
    double within = 0.0, sz = 0.0;
    for(uword j = 0; j < n; j++) {
      if(cl[j] == cl[i]) { within += psm(i, j); sz += 1.0; }
    }
    loss += (lg2(sz) + lg2(s(i)) - 2.0 * lg2(within)) / (double) n;
  }

  return Rcpp::List::create(
    Rcpp::Named("labels") = labels_out,
    Rcpp::Named("value") = loss
  );
}
