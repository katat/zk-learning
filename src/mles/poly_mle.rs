use ark_ff::{Field, Zero};
use ark_poly::{multivariate::{SparsePolynomial, SparseTerm, Term}, DenseMVPolynomial, Polynomial};

use crate::{lagrange::{MultilinearExtension, MultiPoly, UniPoly}, utils::n_to_vec};

#[derive(Debug, Clone)]
pub struct PolyMultilinearExtension<F: Field> {
	evals: Vec<F>,
	p: MultiPoly<F>,
	u: Option<UniPoly<F>>,
}

impl <F: Field> MultilinearExtension<F> for PolyMultilinearExtension<F> {
	fn new(evals: Vec<F>, indexes: Option<Vec<usize>>) -> Self {
		// println!("indexes {:?}", indexes);
		PolyMultilinearExtension {
			evals: evals.clone(),
			p: Self::poly_slow_mle(&evals, &indexes.unwrap()),
			u: None,
		}
	}

	fn fix_vars(&mut self, fixed_vars: &[usize], partial_point: Vec<F>) {
		let point = partial_point;

		if fixed_vars.is_empty() {
			let e = self.p.evaluate(&point);
			self.u = Some(UniPoly::from_coefficients_vec(vec![(0, e)]));
			return
		}

		let p = self.p.terms.clone().into_iter().fold(
			UniPoly::from_coefficients_vec(vec![]),
			|sum, (coeff, term)| {
				let curr = {
					let (coeff_eval, fixed_term) = {
						let mut fixed_term: Option<SparseTerm> = None;
						let var_index = fixed_vars[0];
						let coeff: F =
							term.iter().fold(1u32.into(), |product, (var, power)| match *var {
								j if j == var_index => {
									fixed_term = Some(SparseTerm::new(vec![(j, *power)]));
									product
								}
								j if j < var_index => point[j].pow(&[*power as u64]) * product,
								// _ => point[*var - var_index].pow(&[*power as u64]) * product,
								_ => point[*var].pow(&[*power as u64]) * product,
							});
						(coeff, fixed_term)
					};

					// println!("fixed term {:?} {:?}", fixed_term, coeff_eval);

					match fixed_term {
						None => UniPoly::from_coefficients_vec(vec![(0, coeff * coeff_eval)]),
						_ => UniPoly::from_coefficients_vec(vec![(
							fixed_term.unwrap().degree(),
							coeff * coeff_eval,
						)]),
					}
				};

				sum + curr
			}
		);

		self.u = Some(p);
	}

	fn num_vars(&self) -> usize {
		(self.evals.len() as f64).log2() as usize
	}

	fn evaluate(&self, point: &Vec<F>) -> F {
		self.p.clone().evaluate(point)
	}

	fn to_evals(&self) -> Vec<F> {
		self.evals.clone()
	}

	fn interpolate(&self) -> UniPoly<F> {
		self.u.clone().unwrap()
	}
}

impl <F: Field> PolyMultilinearExtension<F> {
	pub fn poly_constant(c: F) -> SparsePolynomial<F, SparseTerm> {
		SparsePolynomial::from_coefficients_vec(
			0, 
			vec![(c, SparseTerm::new(vec![]))]
		)
	}
	// One step in chi
	pub fn poly_chi_step(w: bool, v: usize) -> SparsePolynomial<F, SparseTerm> {
		let one_minus_w = F::from(1 - (w as u32));
		let w_minus_one = F::from(w as u32) - F::one();
	
		let f: SparsePolynomial<F, SparseTerm> = SparsePolynomial::from_coefficients_vec(
			v + 1,
			vec![
				//x * w
				(w.into(), SparseTerm::new(vec![(v, 1)])),
				// 1 - w
				(one_minus_w, SparseTerm::new(vec![(0, 0)])),
				// (-1 + w) * x
				(w_minus_one, SparseTerm::new(vec![(v, 1)])),
			],
		);
	
		f
	}
	
	/// Perform a naive n^2 multiplication of `self` by `other`.
	pub fn naive_mul(
		cur: &SparsePolynomial<F, SparseTerm>,
		other: &SparsePolynomial<F, SparseTerm>,
	) -> SparsePolynomial<F, SparseTerm> {
		if cur.is_zero() || other.is_zero() {
			// //println!("zero");
			SparsePolynomial::zero()
		} else {
			let mut result_terms = Vec::new();
			for (cur_coeff, cur_term) in cur.terms.iter() {
				for (other_coeff, other_term) in other.terms.iter() {
					let mut term: Vec<(usize, usize)> = cur_term.vars().iter().zip(cur_term.powers()).map(|(v, p)| (*v, p)).collect();
					term.extend(
						other_term.vars().iter().zip(other_term.powers()).map(|(v, p)| (*v, p))
					);
					let coeff = *cur_coeff * *other_coeff;
					result_terms.push((coeff, SparseTerm::new(term)));
				}
			}
			let maxvars = core::cmp::max(cur.num_vars, other.num_vars);
			SparsePolynomial::from_coefficients_slice(maxvars, result_terms.as_slice())
		}
	}
	
	// Computes Chi_w(r) for all w, O(log n) operations
	pub fn poly_chi_w(w: &[bool], vars: &Vec<usize>) -> SparsePolynomial<F, SparseTerm> {
		let product: SparsePolynomial<F, SparseTerm> = w
			.iter()
			.enumerate()
			.zip(vars)
			.fold(
				SparsePolynomial::from_coefficients_slice(0, &[(F::one(), SparseTerm::new(Vec::new()))]), 
				|p, ((_, w), v)| {
					let step = Self::poly_chi_step(*w, *v);
					Self::naive_mul(&p, &step)
				}
			);
	
		product
	}
	
	// Calculating the slow way, for benchmarking
	pub fn poly_slow_mle(evals: &Vec<F>, vars: &Vec<usize>) -> SparsePolynomial<F, SparseTerm> {
		let sum: SparsePolynomial<F, SparseTerm> = evals
			.iter()
			.enumerate()
			.fold(
				SparsePolynomial::from_coefficients_slice(0, &[]),
				|sum, (i, val)| {
					let factor = SparsePolynomial::from_coefficients_slice(0, &[(*val, SparseTerm::new(Vec::new()))]);
					sum + Self::naive_mul(&factor, &Self::poly_chi_w(&n_to_vec(i, vars.len()), vars))
				}
			);
	
		sum
	}
	
	pub fn convert_bin(x: usize, y: usize, n: usize) -> Vec<u32> {
		let xbin = format!("{:0>width$}", format!("{:b}", x), width = n);
		let ybin = format!("{:0>width$}", format!("{:b}", y), width = n);
		let bin = format!("{}{}", xbin, ybin);
		let x: Vec<u32> = bin.chars().map(|x| x.to_digit(10).unwrap())
			.collect();
		x
	}
	
	pub fn convert_bin_z(x: usize, y: usize, z: usize, n: usize) -> Vec<u32> {
		let xbin = format!("{:0>width$}", format!("{:b}", x), width = n);
		let ybin = format!("{:0>width$}", format!("{:b}", y), width = n);
		let zbin = format!("{:0>width$}", format!("{:b}", z), width = n);
		let bin = format!("{}{}{}", xbin, ybin, zbin);
		let x: Vec<u32> = bin.chars().map(|x| x.to_digit(10).unwrap())
			.collect();
		x
	}
	
	pub fn gen_var_indexes (start_index: usize, var_num: usize) -> Vec<usize> {
		let arr: Vec<usize> = (0..var_num).map(|x| x + start_index).collect();
		arr
	}
	
	// One step in chi
	pub fn chi_step(w: bool, x: F) -> F {
		x * F::from(w) + (F::one() - x) * (F::one() - F::from(w))
	}
	
	// Computes Chi_w(r) for all w, O(log n) operations
	pub fn chi_w(w: &Vec<bool>, r: &Vec<F>) -> F {
		assert_eq!(w.len(), r.len());
		let product: F = w
			.iter()
			.zip(r.iter())
			.map(|(&w, &r)| Self::chi_step(w, r))
			.product();
		product
	}
	
	pub fn dynamic_eval(&self, r: &Vec<F>) -> F {
		let chi_lookup = Self::memoize(r, r.len());
		let result: F = self.evals
			.iter()
			.zip(chi_lookup.iter())
			.map(|(left, right)| *left * right)
			.sum();
		result
	}
	
	pub fn memoize(r: &Vec<F>, v: usize) -> Vec<F> {
		match v {
			1 => {
				vec![Self::chi_step(false, r[v - 1]), Self::chi_step(true, r[v - 1])]
			}
			_ => Self::memoize(r, v - 1)
				.iter()
				.flat_map(|val| {
					[
						*val * Self::chi_step(false, r[v - 1]),
						*val * Self::chi_step(true, r[v - 1]),
					]
				})
				.collect(),
		}
	}
}