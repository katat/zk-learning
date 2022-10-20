#[macro_use]
extern crate lazy_static;
use ark_ff::{One, Field};
use ark_poly::{multivariate::{SparsePolynomial, SparseTerm, Term}, Polynomial, DenseMVPolynomial};

use rstest::rstest;
use thaler::lagrange::{self, count_triangles, eval_chi_step, eval_poly_chi_step};
use thaler::small_fields::{F5};

type TestField = F5;

fn convert_field<F: Field>(vec: &[u32]) -> Vec<F> {
	vec.iter().map(|e| F::from(*e)).collect()
}

lazy_static! {
	static ref F_2: Vec<TestField> = convert_field(&[1, 2, 1, 4]);
	static ref R_0: Vec<TestField> = convert_field(&[0, 0]);
	static ref R_1: Vec<TestField> = convert_field(&[0, 2]);
	static ref R_2: Vec<TestField> = convert_field(&[3, 4]);
	static ref R_3: Vec<TestField> = convert_field(&[4, 1]);
	// Larger v tests
	static ref R_4: Vec<TestField> = convert_field(&[0, 0, 0]);
	static ref R_4_CHI: Vec<TestField> = convert_field(&[1, 0, 0, 0, 0, 0, 0, 0]);
	static ref R_5: Vec<TestField> = convert_field(&[1, 2, 3]);
	static ref R_5_CHI: Vec<TestField> = convert_field(&[0, 0, 0, 0, 2, 7, 6, 6]);
	static ref R_6: Vec<TestField> = convert_field(&[5, 5, 5]);
	static ref R_6_CHI: Vec<TestField> = convert_field(&[66, 80, 80, 100, 80, 100, 100, 125]);

	static ref F_W: Vec<TestField> = convert_field(&[
		0, 0, 0, 0,
		0, 1, 1, 0,
		0, 1, 0, 0,
		0, 0, 0, 0
	]);
	static ref F_W_1: Vec<TestField> = convert_field(&[
		0, 1, 1, 0,
		1, 0, 1, 1,
		1, 1, 0, 1,
		0, 1, 1, 0
	]);
}

#[rstest]
#[case(&F_2, &R_0, TestField::from(1))]
#[case(&F_2, &R_1, TestField::from(3))]
#[case(&F_2, &R_2, TestField::from(4))]
#[case(&F_2, &R_3, TestField::from(0))]
fn slow_lagrange_test(
	#[case] fw: &Vec<TestField>,
	#[case] r: &Vec<TestField>,
	#[case] expected: TestField,
) {
	assert_eq!(lagrange::slow_mle(fw, r), expected);
}

#[rstest]
#[case(&F_W_1, TestField::from(2))]
fn count_triangles_test(
	#[case] fw: &Vec<TestField>,
	#[case] expected: TestField
) {
	assert_eq!(count_triangles(fw), expected);

	let e = eval_poly_chi_step(
		true, 
		&SparsePolynomial::from_coefficients_vec(
			1, 
			vec![(TestField::one(), SparseTerm::new(vec![(0, 1)]))]
		)
	);
	assert_eq!(
		eval_chi_step(true, TestField::from(3)),
		e.evaluate(&vec![TestField::from(3)])
	);
}

#[rstest]
#[case(&F_2, &R_0, TestField::from(1))]
#[case(&F_2, &R_1, TestField::from(3))]
#[case(&F_2, &R_2, TestField::from(4))]
#[case(&F_2, &R_3, TestField::from(0))]
fn stream_lagrange_test(
	#[case] fw: &Vec<TestField>,
	#[case] r: &Vec<TestField>,
	#[case] expected: TestField,
) {
	assert_eq!(lagrange::stream_mle(fw, r), expected);
}

#[rstest]
#[case(&F_2, &R_0, TestField::from(1))]
#[case(&F_2, &R_1, TestField::from(3))]
#[case(&F_2, &R_2, TestField::from(4))]
#[case(&F_2, &R_3, TestField::from(0))]
fn dynamic_mle_test(
	#[case] fw: &Vec<TestField>,
	#[case] r: &Vec<TestField>,
	#[case] expected: TestField,
) {
	assert_eq!(lagrange::dynamic_mle(fw, r), expected);
}

#[rstest]
#[case(&R_4, &R_4_CHI)]
#[case(&R_5,&R_5_CHI)]
#[case(&R_6, &R_6_CHI)]
fn memoize_test(#[case] r: &Vec<TestField>, #[case] expected: &Vec<TestField>) {
	assert_eq!(lagrange::memoize(r, r.len()), *expected);
}
