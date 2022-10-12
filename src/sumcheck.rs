use ark_ff::Field;
use ark_ff::{One};
use ark_poly::polynomial::multivariate::{SparsePolynomial, SparseTerm, Term};
use ark_poly::polynomial::univariate::SparsePolynomial as UniSparsePolynomial;
use ark_poly::polynomial::{Polynomial};
use ark_std::cfg_into_iter;

use crate::small_fields::F251;

pub type MultiPoly = SparsePolynomial<F251, SparseTerm>;
pub type UniPoly = UniSparsePolynomial<F251>;

// Converts i into an index in {0,1}^v
pub fn n_to_vec(i: usize, n: usize) -> Vec<F251> {
	format!("{:0>width$}", format!("{:b}", i), width = n)
		.chars()
		.map(|x| if x == '1' { 1.into() } else { 0.into() })
		.collect()
}

pub trait SumCheckPolynomial<F> {
    fn terms(&self) -> Vec<(F, SparseTerm)>;
	fn var_fixed_evaluate<C>(&self, cb: C) -> UniPoly where C: FnMut((F251, SparseTerm)) -> UniPoly;
    fn num_vars(&self) -> usize;
    fn evaluate(&self, point: &Vec<F>) -> F;
}

impl SumCheckPolynomial<F251> for SparsePolynomial<F251, SparseTerm> {
    fn terms(&self) -> Vec<(F251, SparseTerm)> {
		self.terms.to_vec()
    }

    fn var_fixed_evaluate<C>(&self, mut cb: C) -> UniPoly where C: FnMut((F251, SparseTerm)) -> UniPoly {
        self.terms().clone().into_iter().fold(
			UniPoly::from_coefficients_vec(vec![]),
			|sum, term| {
				let curr = cb(term);
				sum + curr
			}
		)

    }

    fn num_vars(&self) -> usize {
		self.num_vars
    }

    fn evaluate(&self, point: &Vec<F251>) -> F251 {
        Polynomial::evaluate(self, point)
    }
}

// Simulates memory of a single prover instance
#[derive(Debug, Clone)]
pub struct Prover<P: SumCheckPolynomial<F251>> {
	pub g: P,
	pub r_vec: Vec<F251>,
}

impl <P: SumCheckPolynomial<F251>> Prover<P> where P: Clone {
	pub fn new(g: &P) -> Self {
		Prover {
			g: g.clone(),
			r_vec: vec![],
		}
	}

	// Given polynomial g, fix Xj, evaluate over xj+1
	pub fn gen_uni_polynomial(&mut self, r: Option<F251>) -> UniPoly {
		if r.is_some() {
			self.r_vec.push(r.unwrap());
		}
		let v = self.g.num_vars() - self.r_vec.len();
		(0..(2u32.pow(v as u32 - 1))).fold(
			UniPoly::from_coefficients_vec(vec![(0, 0u32.into())]),
			|sum, n| {
				let gj = self.evaluate_gj(n_to_vec(n as usize, v));
				sum + gj
			},
		)
	}
	// Evaluates gj over a vector permutation of points, folding all evaluated terms together into one univariate polynomial
	pub fn evaluate_gj(&self, points: Vec<F251>) -> UniPoly {
		self.g.var_fixed_evaluate(|(coeff, term)| {
			let (coeff_eval, fixed_term) = self.evaluate_term(&term, &points);
			match fixed_term {
				None => UniPoly::from_coefficients_vec(vec![(0, coeff * coeff_eval)]),
				_ => UniPoly::from_coefficients_vec(vec![(
					fixed_term.unwrap().degree(),
					coeff * coeff_eval,
				)]),
			}
		})
	}

	// Evaluates a term with a fixed univar, returning (new coefficent, fixed term)
	pub fn evaluate_term(
		&self,
		term: &SparseTerm,
		point: &Vec<F251>,
	) -> (F251, Option<SparseTerm>) {
		let mut fixed_term: Option<SparseTerm> = None;
		let coeff: F251 =
			cfg_into_iter!(term).fold(1u32.into(), |product, (var, power)| match *var {
				j if j == self.r_vec.len() => {
					fixed_term = Some(SparseTerm::new(vec![(j, *power)]));
					product
				}
				j if j < self.r_vec.len() => self.r_vec[j].pow(&[*power as u64]) * product,
				_ => point[*var - self.r_vec.len()].pow(&[*power as u64]) * product,
			});
		(coeff, fixed_term)
	}

	// Sum all evaluations of polynomial `g` over boolean hypercube
	pub fn slow_sum_g(&self) -> F251 {
		let v = self.g.num_vars();
		let n = 2u32.pow(v as u32);
		(0..n)
			.map(|n| self.g.evaluate(&n_to_vec(n as usize, v)))
			.sum()
	}
}

// Verifier procedures
pub fn get_r() -> Option<F251> {
	// let mut rng = rand::thread_rng();
	// let r: F251 = rng.gen();
	let r: F251 = F251::one();
	Some(r)
}

// A degree look up table for all variables in g
pub fn max_degrees<P: SumCheckPolynomial<F251>>(g: &P) -> Vec<usize> {
	let mut lookup: Vec<usize> = vec![0; g.num_vars()];
	cfg_into_iter!(g.terms().clone()).for_each(|(_, term)| {
		cfg_into_iter!(term).for_each(|(var, power)| {
			if *power > lookup[*var] {
				lookup[*var] = *power
			}
		});
	});
	lookup
}

// Verify prover's claim c_1
// Presented pedantically:
pub fn verify<P: SumCheckPolynomial<F251>>(g: &P, c_1: F251) -> bool where P: Clone{
	// 1st round
	let mut p = Prover::new(g);
	let mut gi = p.gen_uni_polynomial(None);
	let mut expected_c = gi.evaluate(&0u32.into()) + gi.evaluate(&1u32.into());
	assert_eq!(c_1, expected_c);
	let lookup_degree = max_degrees(g);
	assert!(gi.degree() <= lookup_degree[0]);

	// middle rounds
	for j in 1..p.g.num_vars() {
		let r = get_r();
		expected_c = gi.evaluate(&r.unwrap());
		gi = p.gen_uni_polynomial(r);
		let new_c = gi.evaluate(&0u32.into()) + gi.evaluate(&1u32.into());
		assert_eq!(expected_c, new_c);
		assert!(gi.degree() <= lookup_degree[j]);
	}
	// final round
	let r = get_r();
	expected_c = gi.evaluate(&r.unwrap());
	p.r_vec.push(r.unwrap());
	let new_c = p.g.evaluate(&p.r_vec);
	assert_eq!(expected_c, new_c);
	true
}

pub fn slow_verify<P: SumCheckPolynomial<F251>>(g: &P, c_1: F251) -> bool where P: Copy{
	let p = Prover::new(g);
	let manual_sum = p.slow_sum_g();
	manual_sum == c_1
}
