use ark_ff::Field;

use crate::{lagrange::{EvaluationMethod}, utils::n_to_vec, mles::value_mle::mle::chi_w};

#[derive(Debug, Clone)]
pub struct StreamEvaluationMethod{}

impl StreamEvaluationMethod {
	fn recurse<F: Field>(fw: &Vec<F>, r: &Vec<F>, n: usize) -> F {
		match n {
			0 => F::zero(),
			_ => Self::recurse(fw, r, n - 1) + fw[n - 1] * chi_w(&n_to_vec(n - 1, r.len()), r),
		}
	}
}

impl <F: Field> EvaluationMethod<F> for StreamEvaluationMethod {
    fn run(evals: &Vec<F>, var_indexes: &Vec<usize>, point: &Vec<F>) -> F{
		let p: Vec<F> = var_indexes.iter().map(|i| point[*i]).collect();

		Self::recurse(evals, &p, 2usize.pow(p.len() as u32))
    }
}