/**
 * @file mtl_test.cpp
 * Test script for interfacing with MTL4 and it's linear solvers.
 */

#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>

#include "IdentityMatrix.hpp"

/** Traits that MTL uses to determine properties of our IdentityMatrix. */
namespace mtl {
namespace ashape {

/** Define IdentityMatrix to be a non-scalar type. */
template <>
struct ashape_aux <IdentityMatrix> {
typedef nonscal type;
};
} // end namespace ashape

/** IdentityMatrix implements the Collection concept with value_type and size_type */
template <>
struct Collection <IdentityMatrix> {
	typedef double value_type;
	typedef unsigned size_type;
};
} // end namespace mtl

using namespace mtl;
using namespace itl;

int main()
{
	const int size = 50, N = size * size;
	IdentityMatrix I = IdentityMatrix(N);
	dense_vector<int> b(N);
	dense_vector<int> x(N, 0);
	b = 15;

	// Preconditionner
	pc::identity<IdentityMatrix> L(I);

    // Termination criterion: r < 1e-6 * b or N iterations
    noisy_iteration<double> iter(b, 500, 1.e-6);
	
	// Use CG solver
	cg(I, x, b, L, iter);
}
