#[macro_use]
extern crate lazy_static;

use ark_poly::polynomial::Polynomial;
use rstest::rstest;
use thaler::{
	lagrange::{self, MultilinearExtension, UniPoly}, 
	mles::{
		value_mle::{
			ValueBasedMultilinearExtension, 
			methods::{
				SlowEvaluationMethod, 
				DynamicEvaluationMethod, 
				StreamEvaluationMethod
			}
		},
	},
	utils::{convert_field, n_to_vec, bools_to_ints}, 
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
	let indexes: Option<Vec<usize>> = Some(vec![0, 1]);
	let mle: ValueBasedMultilinearExtension<TestField, SlowEvaluationMethod> = lagrange::MultilinearExtension::new(fw.clone(), indexes);
	assert_eq!(mle.evaluate(r), expected);
}

#[rstest]
fn test_evaluate_points() {
	let evals = convert_field(&[1, 2, 1, 4]);
	let indexes: Option<Vec<usize>> = Some(vec![0, 1]);
	let mle: ValueBasedMultilinearExtension<TestField, DynamicEvaluationMethod> = lagrange::MultilinearExtension::new(evals, indexes);
	let result: TestField = mle.evaluate(&convert_field(&[0, 0]));
	assert_eq!(result, TestField::from(1));
	let result: TestField = mle.evaluate(&convert_field(&[0, 1]));
	assert_eq!(result, TestField::from(2));
	
	let result: TestField = mle.evaluate(&convert_field(&[1, 0]));
	assert_eq!(result, TestField::from(1));
	let result: TestField = mle.evaluate(&convert_field(&[1, 1]));
	assert_eq!(result, TestField::from(4));


	let mut sum: TestField = TestField::from(0);
	for i in 0..4 {
		let vars = n_to_vec(i, 2);
		let eval = mle.evaluate(&vars.iter().map(|e| TestField::from(*e)).collect());
		sum += eval;
	}

	assert_eq!(sum, TestField::from(8));
}

#[rstest]
#[case(&F_2, 0)]
#[case(&F_2, 1)]
#[case(&F_2, 2)]
#[case(&F_2, 3)]
fn test_fixed_vars_full_point(
	#[case] f: &Vec<TestField>,
	#[case] var_index: usize,
) {
	let num_vars = (f.len() as f64).log2() as usize;
	let indexes: Option<Vec<usize>> = Some(vec![0, 1]);

	// full point
	let mut mle0: ValueBasedMultilinearExtension<TestField, DynamicEvaluationMethod> = MultilinearExtension::new(f.to_vec(), indexes);
	let point: Vec<u32> = bools_to_ints(n_to_vec(var_index, num_vars));
	
	mle0.fix_vars(&[], convert_field(&point));
	assert_eq!(mle0.to_evals(), vec![f[var_index]]);

	// constant polynomial
	let uni: UniPoly<TestField> = mle0.interpolate();
	assert_eq!(uni.degree(), 0);
	assert_eq!(mle0.to_evals().len(), 1);

	for i in 0..F_2.len() {
		assert_eq!(uni.evaluate(&TestField::from(i as u32)), f[var_index]);
	}
}

#[rstest]
#[case(&F_2, &vec![0, 0], 0, &vec![1, 1])] // (x, 0) 
#[case(&F_2, &vec![1, 0], 0, &vec![1, 1])] // (x, 0) 
#[case(&F_2, &vec![0, 1], 0, &vec![2, 4])] // (x, 1) 
#[case(&F_2, &vec![1, 1], 0, &vec![2, 4])] // (x, 1) 
#[case(&F_2, &vec![0, 0], 1, &vec![1, 2])] // (0, y) 
#[case(&F_2, &vec![0, 1], 1, &vec![1, 2])] // (0, y) 
#[case(&F_2, &vec![1, 0], 1, &vec![1, 4])] // (1, y) 
#[case(&F_2, &vec![1, 1], 1, &vec![1, 4])] // (1, y) 
fn test_fixed_vars(
	#[case] f: &Vec<TestField>,
	#[case] xy: &Vec<u32>,
	#[case] var_index: usize,
	#[case] fixed_evals: &Vec<u32>,
) {
	let indexes: Option<Vec<usize>> = Some(vec![0, 1]);

	let mut mle: ValueBasedMultilinearExtension<TestField, DynamicEvaluationMethod> = MultilinearExtension::new(f.to_vec(), indexes.clone());
	mle.fix_vars(&[var_index], convert_field(xy));
	assert_eq!(mle.to_evals(), convert_field(fixed_evals));

	// interpolation
	for i in 0..fixed_evals.len() {
		assert_eq!(
			mle.interpolate().evaluate(&TestField::from(i as u32)), 
			TestField::from(fixed_evals[i])
		);
	}
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
	let indexes: Option<Vec<usize>> = Some(vec![0, 1]);
	let mle: ValueBasedMultilinearExtension<TestField, StreamEvaluationMethod> = MultilinearExtension::new(fw.clone(), indexes);
	assert_eq!(mle.evaluate(r), expected);
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
	let indexes: Option<Vec<usize>> = Some(vec![0, 1]);
	let mle: ValueBasedMultilinearExtension<TestField, DynamicEvaluationMethod> = MultilinearExtension::new(fw.clone(), indexes);
	assert_eq!(mle.evaluate(r), expected);
}

