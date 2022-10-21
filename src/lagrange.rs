use ark_ff::{Zero, Field};
use ark_poly::{multivariate::{SparsePolynomial, SparseTerm, Term}, DenseMVPolynomial};


// One step in chi
pub fn poly_chi_step<F: Field>(w: bool, v: usize) -> SparsePolynomial<F, SparseTerm> {
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

// Converts i into an index in {0,1}^v, used to retrieves f evaluations
pub fn n_to_vec(i: usize, n: usize) -> Vec<bool> {
	let x = format!("{:0>width$}", format!("{:b}", i), width = n);
	let x: Vec<bool> = x.chars().map(|x| x == '1')
		.collect();
	x
}

/// Perform a naive n^2 multiplication of `self` by `other`.
pub fn naive_mul<F: Field>(
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
pub fn poly_chi_w<F: Field>(w: &[bool], vars: &Vec<usize>) -> SparsePolynomial<F, SparseTerm> {
	// //println!("chi w...");
	let product: SparsePolynomial<F, SparseTerm> = w
		.iter()
		.enumerate()
		.zip(vars)
		.fold(
			SparsePolynomial::from_coefficients_slice(0, &[(F::one(), SparseTerm::new(Vec::new()))]), 
			|p, ((_, w), v)| {
				let step = poly_chi_step(*w, *v);
				naive_mul(&p, &step)
			}
		);

	product
}

// Calculating the slow way, for benchmarking
pub fn poly_slow_mle<F: Field>(fw: &[F], vars: &Vec<usize>) -> SparsePolynomial<F, SparseTerm> {
	let sum: SparsePolynomial<F, SparseTerm> = fw
		.iter()
		.enumerate()
		.fold(
			SparsePolynomial::from_coefficients_slice(0, &[]),
			|sum, (i, val)| {
				let factor = SparsePolynomial::from_coefficients_slice(0, &[(*val, SparseTerm::new(Vec::new()))]);
				sum + naive_mul(&factor, &poly_chi_w(&n_to_vec(i, vars.len()), vars))
			}
		);

	sum
}

pub fn poly_constant<F: Field>(c: F) -> SparsePolynomial<F, SparseTerm> {
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
pub fn chi_step<F: Field>(w: bool, x: F) -> F {
	x * F::from(w) + (F::one() - x) * (F::one() - F::from(w))
}

// Computes Chi_w(r) for all w, O(log n) operations
pub fn chi_w<F: Field>(w: &Vec<bool>, r: &Vec<F>) -> F {
	assert_eq!(w.len(), r.len());
	let product: F = w
		.iter()
		.zip(r.iter())
		.map(|(&w, &r)| chi_step(w, r))
		.product();
	product
}

// Calculating the slow way, for benchmarking
pub fn slow_mle<F: Field>(fw: &Vec<F>, r: &Vec<F>) -> F {
	assert_eq!(r.len() as f64, (fw.len() as f64).log2());
	let sum: F = fw
		.iter()
		.enumerate()
		.map(|(i, val)| *val * chi_w(&n_to_vec(i, r.len()), r))
		.sum();
	sum
}

// Lemma 3.7
pub fn stream_mle<F: Field>(fw: &Vec<F>, r: &Vec<F>) -> F {
	recurse(fw, r, 2usize.pow(r.len() as u32))
}

pub fn recurse<F: Field>(fw: &Vec<F>, r: &Vec<F>, n: usize) -> F {
	match n {
		0 => F::zero(),
		_ => recurse(fw, r, n - 1) + fw[n - 1] * chi_w(&n_to_vec(n - 1, r.len()), r),
	}
}

// Lemm 3.8
pub fn dynamic_mle<F: Field>(fw: &[F], r: &Vec<F>) -> F {
	let chi_lookup = memoize(r, r.len());
	let result: F = fw
		.iter()
		.zip(chi_lookup.iter())
		.map(|(left, right)| *left * right)
		.sum();
	result
}

pub fn memoize<F: Field>(r: &Vec<F>, v: usize) -> Vec<F> {
	match v {
		1 => {
			vec![chi_step(false, r[v - 1]), chi_step(true, r[v - 1])]
		}
		_ => memoize(r, v - 1)
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
