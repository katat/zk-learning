#[macro_use]
extern crate lazy_static;

use rstest::rstest;
use thaler::mles::{
	value_mle::{
		ValueBasedMultilinearExtension,
		methods::DynamicEvaluationMethod
	},
	PolyMultilinearExtension
};
use thaler::small_fields::{F251};
use thaler::sumcheck;
use thaler::triangles::{TriangleMLE, TriangleGraph};

type TestField = F251;

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
	static ref M: TriangleGraph<TestField> = TriangleGraph::new(thaler::utils::gen_matrix(4));

	static ref G_2: TriangleMLE<TestField, ValueBasedMultilinearExtension<TestField, DynamicEvaluationMethod>> = M.derive_mle();
	static ref G_2_SUM: TestField = sumcheck::Prover::<TestField, TriangleMLE<TestField, ValueBasedMultilinearExtension<TestField, DynamicEvaluationMethod>>>::new(&G_2).slow_sum_g();

	static ref G_3: TriangleMLE<TestField, PolyMultilinearExtension<TestField>> = M.derive_mle();
	static ref G_3_SUM: TestField = sumcheck::Prover::<TestField, TriangleMLE<TestField, PolyMultilinearExtension<TestField>>>::new(&G_3).slow_sum_g();

}

#[rstest]
#[case(&G_2, &G_2_SUM)]
fn sumcheck_triangles_test(#[case] p: &TriangleMLE<TestField, ValueBasedMultilinearExtension<TestField, DynamicEvaluationMethod>>, #[case] c: &TestField) {
	assert!(sumcheck::verify::<_, TriangleMLE<_, ValueBasedMultilinearExtension<_, DynamicEvaluationMethod>>>(p, *c));
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
	let g: TriangleMLE<TestField, ValueBasedMultilinearExtension<TestField, DynamicEvaluationMethod>> = matrix.derive_mle();
	let sum: TestField = sumcheck::Prover::<TestField, TriangleMLE<TestField, ValueBasedMultilinearExtension<TestField, DynamicEvaluationMethod>>>::new(&g).slow_sum_g();

	assert!(sumcheck::verify::<TestField, TriangleMLE<TestField, ValueBasedMultilinearExtension<TestField, DynamicEvaluationMethod>>>(&g, sum));
}

