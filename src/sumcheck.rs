use ark_ff::{Field};
use ark_poly::polynomial::{Polynomial};
use rand::Rng;

use crate::lagrange::UniPoly;

#[derive(Debug, Clone, Eq, PartialEq)]
pub enum Round {
	Middle(usize),
	Final(),
}

//todo use the one from utils instead
// Converts i into an index in {0,1}^v
pub fn n_to_vec<F: Field>(i: usize, n: usize) -> Vec<F> {
	format!("{:0>width$}", format!("{:b}", i), width = n)
		.chars()
		.map(|x| if x == '1' { F::one() } else { F::zero() })
		.collect()
}

pub trait SumCheckPolynomial<F> where F: Field {
	fn var_fixed_evaluate(&self, var: usize, point: Vec<F>) -> UniPoly<F>;
    fn num_vars(&self) -> usize;
    fn evaluate(&self, point: &Vec<F>) -> F;
}


// Simulates memory of a single prover instance
#[derive(Debug, Clone)]
pub struct Prover<F, P> where F: Field, P: SumCheckPolynomial<F> {
	pub g: P,
	pub r_vec: Vec<F>,
}

#[derive(Debug, Clone)]
pub struct Verifier<F: Field, P: Polynomial<F>, G: SumCheckPolynomial<F>> {
	pub c_1: F,
	pub g: G,
	pub rounds: usize,
	pub prev_g_i: Option<P>,
	pub r_vec: Vec<F>,
	pub current_round: Option<Round>,
	//todo this is not closure, find a better way to achieve closure.
	pub random_func: Option<fn() -> F>
}

impl <F: Field, P: Polynomial<F, Point = F>, G: SumCheckPolynomial<F>> Verifier<F, P, G> {
	pub fn new(c_1: F, g: G) -> Self {
		let rounds = g.num_vars();
		Verifier {
			c_1,
			g,
			rounds,
			prev_g_i: None,
			r_vec: vec![],
			current_round: None,
			random_func: None,
		}
	}

	pub fn random_func(&mut self, f: fn() -> F) {
		self.random_func = Some(f);
	}

	fn generate_random(&mut self) -> F {
		let mut rng = rand::thread_rng();
		let mut r = F::from(rng.gen::<u32>());

		if self.random_func != None {
			r = self.random_func.unwrap()();
		}

		self.r_vec.push(r);
		r
	}

	fn advance(&mut self) -> Round {
		self.generate_random();
		let cur_round = self.r_vec.len();

		if cur_round < self.rounds {
			Round::Middle(cur_round)
		}
		else {
			Round::Final()
		}
	}

	pub fn verify(&mut self, g_i: Option<P>) {
		match &self.current_round {
			None => {
				let g_i = g_i.unwrap();
				let new_c = g_i.evaluate(&0u32.into()) + g_i.evaluate(&1u32.into());
				assert_eq!(self.c_1, new_c);
				self.prev_g_i = Some(g_i);
			},
			Some(Round::Middle(_)) => {
				let g_i = g_i.unwrap();
				let last_r = self.r_vec.last().unwrap();
				let expected_c = self.prev_g_i.as_ref().unwrap().evaluate(last_r);
				let new_c = g_i.evaluate(&0u32.into()) + g_i.evaluate(&1u32.into());
				assert_eq!(expected_c, new_c);
				self.prev_g_i = Some(g_i);
			},
			Some(Round::Final()) => {
				let last_r = self.r_vec.last().unwrap();
				let expected_c = self.prev_g_i.as_ref().unwrap().evaluate(last_r);
				let new_c = self.g.evaluate(&self.r_vec);
				assert_eq!(expected_c, new_c);
				return
			}
		};

		self.current_round = Some(self.advance());
	}
}

impl <F: Field, P: SumCheckPolynomial<F>> Prover<F, P> where P: Clone {
	pub fn new(g: &P) -> Self {
		Prover {
			g: g.clone(),
			r_vec: vec![],
		}
	}

	// Given polynomial g, fix Xj, evaluate over xj+1
	pub fn gen_uni_polynomial(&mut self, r: Option<F>) -> UniPoly<F> {
		if r.is_some() {
			self.r_vec.push(r.unwrap());
		}
		let v = self.g.num_vars() - self.r_vec.len();
		(0..(2u32.pow(v as u32 - 1))).fold(
			UniPoly::from_coefficients_vec(vec![(0, 0u32.into())]),
			|sum, n| {
				let point = n_to_vec(n as usize, v);
				let gj = self.evaluate_gj(point);
				sum + gj
			},
		)
	}

	// Evaluates gj over a vector permutation of points, folding all evaluated terms together into one univariate polynomial
	pub fn evaluate_gj(&self, points: Vec<F>) -> UniPoly<F> {
		self.g.var_fixed_evaluate(self.r_vec.len(), [self.r_vec.to_vec(), points].concat())
	}
}


