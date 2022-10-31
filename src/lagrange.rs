use ark_ff::{Field};
use crate::{sumcheck::UniPoly};


pub trait MultilinearExtension<F>: Clone {
	fn new(evals: Vec<F>, var_indexes: Option<Vec<usize>>) -> Self;
	fn num_vars(&self) -> usize;
	fn fix_vars(&mut self, fixed_vars: &[usize], partial_point: Vec<F>);
	fn evaluate(&self, point: &Vec<F>) -> F;
	fn to_evals(&self) -> Vec<F>;
	fn interpolate(&self) -> UniPoly<F> where F: Field;
}

