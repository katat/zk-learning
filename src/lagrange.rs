use ark_ff::{Field};
use ark_poly::multivariate::{SparsePolynomial, SparseTerm};
use ark_poly::univariate::{SparsePolynomial as UniSparsePolynomial};

pub type MultiPoly<F> = SparsePolynomial<F, SparseTerm>;
pub type UniPoly<F> = UniSparsePolynomial<F>;

pub trait MultilinearExtension<F>: Clone {
	fn new(evals: Vec<F>, var_indexes: Option<Vec<usize>>) -> Self;
	fn num_vars(&self) -> usize;
	fn fix_vars(&mut self, fixed_vars: &[usize], partial_point: Vec<F>);
	fn evaluate(&self, point: &Vec<F>) -> F;
	fn to_evals(&self) -> Vec<F>;
	fn interpolate(&self) -> UniPoly<F> where F: Field;
}

pub trait EvaluationMethod<F>: Clone {
	fn run(evals: &Vec<F>, var_indexes: &Vec<usize>, point: &Vec<F>) -> F;
}

