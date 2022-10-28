#[macro_use]
extern crate lazy_static;

use ark_poly::DenseMVPolynomial;
use ark_poly::polynomial::multivariate::{SparsePolynomial, SparseTerm, Term};
use rstest::rstest;
use thaler::small_fields::{F251};
use thaler::sumcheck::{self, MultiPoly};
use thaler::triangles::{TriangleMLE, TriangleGraph, MLEAlgorithm};

fn convert_vec (m: &[i32]) -> Vec<Vec<F251>> {
	let len = (m.len() as f64).sqrt() as usize;

	let mut matrix = Vec::new();
	for (i, e) in m.iter().enumerate() {
		if i % len == 0 {
			let row = Vec::new();
			matrix.push(row);
		}
		let v = F251::from(*e);
		let row = matrix.last_mut().unwrap();
		
		row.push(v);
	}

	matrix
}

lazy_static! {
	static ref G_0: MultiPoly = SparsePolynomial::from_coefficients_vec(
		3,
		vec![
			(2u32.into(), SparseTerm::new(vec![(0, 3)])),
			(1u32.into(), SparseTerm::new(vec![(0, 1), (2, 1)])),
			(1u32.into(), SparseTerm::new(vec![(1, 1), (2, 1)])),
		],
	);
	static ref G_0_SUM: F251 = sumcheck::Prover::<MultiPoly>::new(&G_0).slow_sum_g();
	// Test with a larger g
	static ref G_1: sumcheck::MultiPoly = SparsePolynomial::from_coefficients_vec(
		4,
		vec![
			(2u32.into(), SparseTerm::new(vec![(0, 3)])),
			(1u32.into(), SparseTerm::new(vec![(0, 1), (2, 1)])),
			(1u32.into(), SparseTerm::new(vec![(1, 1), (2, 1)])),
			(1u32.into(), SparseTerm::new(vec![(3, 1), (2, 1)])),
		],
	);
	static ref G_1_SUM: F251 = sumcheck::Prover::<MultiPoly>::new(&G_1).slow_sum_g();

	static ref M: TriangleGraph<F251> = TriangleGraph::new(thaler::utils::gen_matrix(4));
	static ref G_2: TriangleMLE<F251> = M.derive_mle(MLEAlgorithm::DynamicMLE);
	static ref G_2_SUM: F251 = sumcheck::Prover::<TriangleMLE<F251>>::new(&G_2).slow_sum_g();

}

#[rstest]
#[case(&G_0, &G_0_SUM)]
#[case(&G_1, &G_1_SUM)]
fn sumcheck_multi_poly_test(#[case] p: &sumcheck::MultiPoly, #[case] c: &F251) {
	assert!(sumcheck::verify::<MultiPoly>(p, *c));
}

#[rstest]
#[case(&G_2, &G_2_SUM)]
fn sumcheck_triangles_test(#[case] p: &TriangleMLE<F251>, #[case] c: &F251) {
	assert!(sumcheck::verify::<TriangleMLE<F251>>(p, *c));
}

#[rstest]
fn sumcheck_triangles_2_test() {
	let m = convert_vec(&[0, 1, 1, 0]);
    let matrix = TriangleGraph::new(m.to_vec());
	let g = matrix.derive_mle(MLEAlgorithm::SlowMLE);
	let sum: F251 = sumcheck::Prover::<TriangleMLE<F251>>::new(&g).slow_sum_g();

	assert!(sumcheck::verify::<TriangleMLE<F251>>(&g, sum));
}

