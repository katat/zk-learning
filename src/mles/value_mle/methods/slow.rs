use ark_ff::Field;

use crate::{lagrange::EvaluationMethod, utils::n_to_vec, mles::value_mle::mle::chi_w};

#[derive(Debug, Clone)]
pub struct SlowEvaluationMethod{}

impl <F: Field> EvaluationMethod<F> for SlowEvaluationMethod {
    fn run(evals: &Vec<F>, var_indexes: &Vec<usize>, point: &Vec<F>) -> F{
		let p: Vec<F> = var_indexes.iter().map(|i| point[*i]).collect();

		assert_eq!(p.len() as f64, (evals.len() as f64).log2());
		let sum: F = evals
			.iter()
			.enumerate()
			.map(|(i, val)| *val * chi_w(&n_to_vec(i, p.len()), &p))
			.sum();
		sum
    }
}