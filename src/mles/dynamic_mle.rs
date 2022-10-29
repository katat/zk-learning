use ark_ff::{Field, Zero};
use ark_poly::{multivariate::{SparsePolynomial, SparseTerm, Term}, DenseMVPolynomial};

use crate::{lagrange::MultilinearExtension, sumcheck::UniPoly, utils::n_to_vec};

#[derive(Debug, Clone)]
pub struct DynamicMultilinearExtension<F: Field> {
	evals: Vec<F>,
}

impl <F: Field> MultilinearExtension<F> for DynamicMultilinearExtension<F> {
	fn new(evals: Vec<F>) -> Self {
		DynamicMultilinearExtension {
			evals,
		}
	}

	fn fix_vars(&mut self, fixed_vars: &[usize], partial_point: Vec<F>){
		let mut points = Vec::<Vec<F>>::new();
		let mut evals = Vec::<F>::new();

		match fixed_vars.len() {
			0 => {
				points.push(partial_point);
			}
			_ => {
				for var in fixed_vars {
					for b in [F::zero(), F::one()] {
						let mut point = partial_point.clone();
						let loc = *var;
						
						if (self.to_evals().len() as f64).log2() != partial_point.len() as f64 {
							point.splice(loc..loc, vec![b].iter().cloned());
						}
						else {
							point.splice(loc..loc+1, vec![b].iter().cloned());
						}		
						points.push(point);
					}
				};
			}
		}

		// println!("points {:?}", points);

		for point in points {
			evals.push(self.evaluate(&point));
		}

		// println!("evals {:?}", evals);

		self.evals = evals;
	}

	fn num_vars(&self) -> usize {
		(self.evals.len() as f64).log2() as usize
	}

	fn evaluate(&self, point: &Vec<F>) -> F {
		self.dynamic_eval(point)
	}

	fn to_evals(&self) -> Vec<F> {
		self.evals.clone()
	}

	fn interpolate(&self) -> UniPoly<F> {
		let mut sum = UniPoly::from_coefficients_vec(vec![(0, F::zero())]);

		for i in 0..self.evals.len() {
			let e = self.evals[i];
			let eval = UniPoly::from_coefficients_vec(vec![(0, e)]);
			
			let mut cum_prod = UniPoly::from_coefficients_vec(vec![(0, F::one())]);
			for k in 0..self.evals.len() {
				
				if i == k {
					continue;
				}
				let k = F::from(k as u32);
				// (x - k) / (i - k) = x/(i - k) - k/(i - k)

				// i - k
				let ik = F::from(i as u32) - k;

				// x/(i - k)
				let x_ik = (1, F::one() / ik);

				// k/(i - k)
				let k_ik = k / ik;

				// (x - k) / (i - k)
				let prod = UniPoly::from_coefficients_vec(vec![x_ik, (0, k_ik.neg())]);
				cum_prod = cum_prod.mul(&prod);
			}
			// mul eval
			sum += &eval.mul(&cum_prod);
		}
		sum
	}

}

impl <F: Field> DynamicMultilinearExtension<F> {
	// pub fn add(&self, other: SlowMultilinearExtension<F>) -> SlowMultilinearExtension<F> {
	// 	assert_eq!(self.to_evals().len(), other.to_evals().len());
	// 	// single row addition
	// 	// todo support matrix 
	// 	let added_evals = self.to_evals().iter().zip(other.to_evals().iter()).map(|(a, b)| {
	// 		a.add(b)
	// 	}).collect();

	// 	Self {
	// 		evals: added_evals,
	// 		algo: self.algo.clone(),
	// 	}
	// }

	// pub fn mul(&self, other: SlowMultilinearExtension<F>) -> SlowMultilinearExtension<F> {
	// 	assert_eq!(self.to_evals().len(), other.to_evals().len());
	// 	// single row multiplication
	// 	// todo support matrix
	// 	let multiplicated_evals = self.to_evals().iter().zip(other.to_evals().iter()).map(|(a, b)| {
	// 		a.mul(b)
	// 	}).collect();

	// 	Self {
	// 		evals: multiplicated_evals,
	// 		algo: self.algo.clone(),
	// 	}
	// }

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
	pub fn poly_slow_mle(&self, vars: &Vec<usize>) -> SparsePolynomial<F, SparseTerm> {
		let sum: SparsePolynomial<F, SparseTerm> = self.evals
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
	
	pub fn poly_constant(c: F) -> SparsePolynomial<F, SparseTerm> {
		SparsePolynomial::from_coefficients_vec(
			0, 
			vec![(c, SparseTerm::new(vec![]))]
		)
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