#include "beachmat/numeric_matrix.h"
#include "pcg_random.hpp"
#include "convert_seed.h"
#include "boost/random.hpp"
#include "boost/range/algorithm.hpp"

//' @importFrom Rcpp sourceCpp
//' @useDynLib PCAtools
// [[Rcpp::export(rng=false)]]
SEXP shuffle_matrix(SEXP incoming, SEXP seed, int stream) {
    auto in=beachmat::create_matrix<beachmat::numeric_matrix>(incoming);
    const size_t NR=in->get_nrow(), NC=in->get_ncol();
    auto out=beachmat::create_output<beachmat::numeric_output>(NR, NC, beachmat::output_param(in.get()));
    typename Rcpp::NumericVector tmp(NR);

    auto gen=pcg32(dqrng::convert_seed<uint64_t>(seed), stream);
    for (size_t c=0; c<NC; ++c) {
        in->get_col(c, tmp.begin());
        boost::range::random_shuffle(tmp, gen);
        out->set_col(c, tmp.begin());
    }

    return out->yield();
}
