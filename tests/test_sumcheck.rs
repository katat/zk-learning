#[macro_use]
extern crate lazy_static;

use ark_poly::DenseMVPolynomial;
use ark_poly::polynomial::multivariate::{SparsePolynomial, SparseTerm, Term};
use rstest::rstest;
use thaler::mles::dynamic_mle::DynamicMultilinearExtension;
use thaler::mles::poly_mle::PolyMultilinearExtension;
use thaler::small_fields::{F251};
use thaler::sumcheck::{self, MultiPoly};
use thaler::triangles::{TriangleMLE, TriangleGraph};

type TestField = F251;
type SumCheckPoly = MultiPoly<TestField>;

fn convert_vec (m: &[i32]) -> Vec<Vec<TestField>> {
	let len = (m.len() as f64).sqrt() as usize;

	let mut matrix = Vec::new();
	for (i, e) in m.iter().enumerate() {
		if i % len == 0 {
			let row = Vec::new();
			matrix.push(row);
		}
		let v = TestField::from(*e);
		let row = matrix.last_mut().unwrap();
		
		row.push(v);
	}

	matrix
}

lazy_static! {
	static ref G_0: SumCheckPoly = SparsePolynomial::from_coefficients_vec(
		3,
		vec![
			(2u32.into(), SparseTerm::new(vec![(0, 3)])),
			(1u32.into(), SparseTerm::new(vec![(0, 1), (2, 1)])),
			(1u32.into(), SparseTerm::new(vec![(1, 1), (2, 1)])),
		],
	);
	static ref G_0_SUM: TestField = sumcheck::Prover::<TestField, SumCheckPoly>::new(&G_0).slow_sum_g();
	// Test with a larger g
	static ref G_1: SumCheckPoly = SparsePolynomial::from_coefficients_vec(
		4,
		vec![
			(2u32.into(), SparseTerm::new(vec![(0, 3)])),
			(1u32.into(), SparseTerm::new(vec![(0, 1), (2, 1)])),
			(1u32.into(), SparseTerm::new(vec![(1, 1), (2, 1)])),
			(1u32.into(), SparseTerm::new(vec![(3, 1), (2, 1)])),
		],
	);
	static ref G_1_SUM: TestField = sumcheck::Prover::<TestField, SumCheckPoly>::new(&G_1).slow_sum_g();

	static ref M: TriangleGraph<TestField> = TriangleGraph::new(thaler::utils::gen_matrix(4));

	static ref G_2: TriangleMLE<TestField, DynamicMultilinearExtension<TestField>> = M.derive_mle();
	static ref G_2_SUM: TestField = sumcheck::Prover::<TestField, TriangleMLE<TestField, DynamicMultilinearExtension<TestField>>>::new(&G_2).slow_sum_g();

	static ref G_3: TriangleMLE<TestField, PolyMultilinearExtension<TestField>> = M.derive_mle();
	static ref G_3_SUM: TestField = sumcheck::Prover::<TestField, TriangleMLE<TestField, PolyMultilinearExtension<TestField>>>::new(&G_3).slow_sum_g();

}

#[rstest]
#[case(&G_0, &G_0_SUM)]
#[case(&G_1, &G_1_SUM)]
fn sumcheck_multi_poly_test(#[case] p: &SumCheckPoly, #[case] c: &TestField) {
	assert!(sumcheck::verify::<TestField, SumCheckPoly>(p, *c));
}

#[rstest]
#[case(&G_2, &G_2_SUM)]
fn sumcheck_triangles_test(#[case] p: &TriangleMLE<TestField, DynamicMultilinearExtension<TestField>>, #[case] c: &TestField) {
	assert!(sumcheck::verify::<_, TriangleMLE<_, DynamicMultilinearExtension<_>>>(p, *c));
}

#[rstest]
#[case(&G_3, &G_3_SUM)]
fn sumcheck_triangles_poly_test(#[case] p: &TriangleMLE<TestField, PolyMultilinearExtension<TestField>>, #[case] c: &TestField) {
	assert!(sumcheck::verify::<_, TriangleMLE<_, PolyMultilinearExtension<_>>>(p, *c));
}

#[rstest]
fn sumcheck_triangles_2_test() {
	let m = convert_vec(&[0, 1, 1, 0]);
    let matrix = TriangleGraph::new(m.to_vec());
	let g: TriangleMLE<TestField, DynamicMultilinearExtension<TestField>> = matrix.derive_mle();
	let sum: TestField = sumcheck::Prover::<TestField, TriangleMLE<TestField, DynamicMultilinearExtension<TestField>>>::new(&g).slow_sum_g();

	assert!(sumcheck::verify::<TestField, TriangleMLE<TestField, DynamicMultilinearExtension<TestField>>>(&g, sum));
}

