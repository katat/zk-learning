use ark_ff::Field;

use crate::{lagrange::{EvaluationMethod}, mles::value_mle::mle::chi_step};

#[derive(Debug, Clone)]
pub struct DynamicEvaluationMethod{}

impl DynamicEvaluationMethod {
	fn memoize<F: Field>(r: &Vec<F>, v: usize) -> Vec<F> {
		match v {
			1 => {
				vec![
					chi_step(false, r[v - 1]), 
					chi_step(true, r[v - 1])
				]
			}
			_ => Self::memoize(r, v - 1)
				.iter()
				.flat_map(|val| {
					[
						*val * chi_step(false, r[v - 1]),
						*val * chi_step(true, r[v - 1]),
					]
				})
				.collect(),
		}
	}
}

impl <F: Field> EvaluationMethod<F> for DynamicEvaluationMethod {
    fn run(evals: &Vec<F>, var_indexes: &Vec<usize>, point: &Vec<F>) -> F{
		let p: Vec<F> = var_indexes.iter().map(|i| point[*i]).collect();

		let chi_lookup = Self::memoize(&p, p.len());
		let result: F = evals
			.iter()
			.zip(chi_lookup.iter())
			.map(|(left, right)| *left * right)
			.sum();
		result
    }
}