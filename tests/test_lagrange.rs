#[macro_use]
extern crate lazy_static;

use ark_ff::{Zero, One};
use ark_poly::polynomial::Polynomial;
use rstest::rstest;
use thaler::{
	lagrange::{self, MultilinearExtension}, 
	mles::{
		slow_mle::SlowMultilinearExtension, 
		dynamic_mle::DynamicMultilinearExtension, 
		stream_mle::StreamMultilinearExtension,
	},
	utils::{convert_field, n_to_vec}, 
	sumcheck::UniPoly
};
use thaler::small_fields::{F5};

type TestField = F5;

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
	let mle: SlowMultilinearExtension<TestField> = lagrange::MultilinearExtension::new(fw.clone());
	assert_eq!(mle.evaluate(r), expected);
}

#[rstest]
fn t() {
	let evals = convert_field(&[2, 4, 3, 2]);
	let mle: DynamicMultilinearExtension<TestField> = lagrange::MultilinearExtension::new(evals);
	let result: TestField = mle.evaluate(&convert_field(&[0, 0]));
	assert_eq!(result, TestField::from(2));
	let result: TestField = mle.evaluate(&convert_field(&[0, 1]));
	assert_eq!(result, TestField::from(4));
	
	let result: TestField = mle.evaluate(&convert_field(&[1, 0]));
	assert_eq!(result, TestField::from(3));
	let result: TestField = mle.evaluate(&convert_field(&[1, 1]));
	assert_eq!(result, TestField::from(2));


	let mut sum: TestField = TestField::from(0);
	for i in 0..4 {
		let vars = n_to_vec(i, 2);
		println!("{:?}", vars);
		let eval = mle.evaluate(&vars.iter().map(|e| TestField::from(*e)).collect());
		println!("eval {:?}", eval);
		sum += eval;
	}

	assert_eq!(sum, TestField::from(11));

	// 2 * 4 * 3 * 2 = 48
	// [6,8]
}

#[rstest]
fn test_fix_vars() {
	let evals: Vec<TestField> = convert_field(&[2, 4, 3, 2]);

	// full point
	let mut mle0: DynamicMultilinearExtension<TestField> = MultilinearExtension::new(evals.clone());
	mle0.fix_vars(&[], convert_field(&[0, 1]));
	assert_eq!(mle0.to_evals(), convert_field(&[4]));

	// x1
	let mut mle1: DynamicMultilinearExtension<TestField> = MultilinearExtension::new(evals.clone());
	mle1.fix_vars(&[0], [TestField::one()].to_vec());
	assert_eq!(mle1.to_evals(), convert_field(&[4, 2]));

	// x1 replace full point
	let mut mle1: DynamicMultilinearExtension<TestField> = MultilinearExtension::new(evals.clone());
	mle1.fix_vars(&[0], [TestField::zero(), TestField::one()].to_vec());
	assert_eq!(mle1.to_evals(), convert_field(&[4, 2]));

	// interpolate for non-constant polynomial
	let uni: UniPoly<TestField> = mle1.interpolate();
	assert_eq!(uni.evaluate(&TestField::zero()), TestField::from(4));

	let uni: UniPoly<TestField> = mle0.interpolate();
	assert_eq!(uni.degree(), 0);
	assert_eq!(mle0.to_evals().len(), 1);
	assert_eq!(uni.evaluate(&TestField::one()), mle0.to_evals()[0]);

	// x0
	let mut mle2: DynamicMultilinearExtension<TestField> = MultilinearExtension::new(evals.clone());
	mle2.fix_vars(&[0], [TestField::zero()].to_vec());
	assert_eq!(mle2.to_evals(), convert_field(&[2, 3]));

	// 2 * 4 * 3 * 2
	// let mle3 = mle1.clone().mul(mle2.clone());
	// let result = mle3.evaluate(&convert_field(&[0])) * mle3.evaluate(&convert_field(&[1]));
	// assert_eq!(result, evals[0] * evals[1] * evals[2] * evals[3]);

	// let mle3 = mle1.clone().add(mle2.clone());
	// let result = mle3.evaluate(&convert_field(&[0])) + mle3.evaluate(&convert_field(&[1]));
	// assert_eq!(result, evals[0] + evals[1] + evals[2] + evals[3]);

	// support 
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
	let mle: StreamMultilinearExtension<TestField> = MultilinearExtension::new(fw.clone());
	assert_eq!(mle.stream_eval(r), expected);
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
	let mle = DynamicMultilinearExtension::new(fw.clone());
	assert_eq!(mle.dynamic_eval(r), expected);
}

#[rstest]
#[case(&R_4, &R_4_CHI)]
#[case(&R_5,&R_5_CHI)]
#[case(&R_6, &R_6_CHI)]
fn memoize_test(#[case] r: &Vec<TestField>, #[case] expected: &Vec<TestField>) {
	assert_eq!(DynamicMultilinearExtension::memoize(r, r.len()), *expected);
}
