#[macro_use]
extern crate lazy_static;

use ark_ff::Field;
use rstest::rstest;
use thaler::{mles::{
	value_mle::{
		ValueBasedMultilinearExtension,
		methods::DynamicEvaluationMethod
	},
	PolyMultilinearExtension
}, sumcheck::{SumCheckPolynomial, Prover, Verifier, Round}, lagrange::UniPoly};
use thaler::small_fields::{F251};
use thaler::triangles::{TriangleMLE, TriangleGraph};

type TestField = F251;

type DynamicMLE<F> = ValueBasedMultilinearExtension<F, DynamicEvaluationMethod>;
type TriangleDynamicMLE<F> = TriangleMLE<F, DynamicMLE<F>>;

type PolyMLE<F> = PolyMultilinearExtension<F>;
type TrianglePolyMLE<F> = TriangleMLE<F, PolyMLE<F>>;

lazy_static! {
	static ref M: TriangleGraph<TestField> = TriangleGraph::new(thaler::utils::gen_matrix(4));

	static ref G_1: TriangleDynamicMLE<TestField> = M.derive_mle();
	static ref G_1_SUM: TestField = G_1.hypercube_sum();

	static ref G_2: TrianglePolyMLE<TestField> = M.derive_mle();
	static ref G_2_SUM: TestField = G_2.hypercube_sum();
}

fn run_protocol<F: Field, P: SumCheckPolynomial<F>>(g: &P, c: F) -> bool where P: Clone{
	let mut p = Prover::new(g);

	let mut v: Verifier<F, UniPoly<F>, P> = Verifier::new(c, g.to_owned());
	while v.current_round != Some(Round::Final()) {
		let r = v.r_vec.last();
		let gi = match r {
			None => p.gen_uni_polynomial(None),
			_default => p.gen_uni_polynomial(Some(*r.unwrap()))
		};

		v.verify(Some(gi));
	}

	// final round
	v.verify(None);
	true
}

#[rstest]
#[case(&G_1, &G_1_SUM)]
fn triangle_graph_count_value_mle_test(#[case] p: &TriangleDynamicMLE<TestField>, #[case] c: &TestField) {
	assert_eq!(p.matrix.count() * TestField::from(6), *c);
}

#[rstest]
#[case(&G_2, &G_2_SUM)]
fn triangle_graph_count_poly_mle_test(#[case] p: &TrianglePolyMLE<TestField>, #[case] c: &TestField) {
	assert_eq!(p.matrix.count() * TestField::from(6), *c);
}


#[rstest]
#[case(&G_1, &G_1_SUM)]
fn sumcheck_triangle_dynamic_mle_test(#[case] p: &TriangleDynamicMLE<TestField>, #[case] c: &TestField) {
	assert!(run_protocol::<_, TriangleDynamicMLE<_>>(p, *c));
}

#[rstest]
#[case(&G_2, &G_2_SUM)]
fn sumcheck_triangle_poly_mle_test(#[case] p: &TrianglePolyMLE<TestField>, #[case] c: &TestField) {
	assert!(run_protocol::<_, TrianglePolyMLE<_>>(p, *c));
}
