use ark_ff::{Zero, PrimeField};
use ark_poly::{multivariate::{SparsePolynomial, SparseTerm, Term}, MVPolynomial, Polynomial};
use crate::small_fields::F5;


// One step in chi
pub fn chi_step(w: bool, v: usize) -> SparsePolynomial<F5, SparseTerm> {
	// x * i128::from(w) + (1 - x) * (1 - i128::from(w));
	// println!("step v+1:{}", v+1);

	let one_minus_w = (1 - (w as i128)).into();
	// println!("1-w: {}", one_minus_w);
	let w_minus_one = w as i128 - 1;
	// println!("w-1: {}", w_minus_one);

	let f: SparsePolynomial<F5, SparseTerm> = SparsePolynomial::from_coefficients_vec(
		v + 1,
		vec![
			//x * i128::from(w)
			(w.into(), SparseTerm::new(vec![(v, 1)])),
			// 1 - w
			(one_minus_w, SparseTerm::new(vec![(0, 0)])),
			// (-1 + w) * x
			(w_minus_one.into(), SparseTerm::new(vec![(v, 1)])),
		],
	);

	// println!("2nd term coeff: {:?}", f.terms()[1]);

	f
}

// Converts i into an index in {0,1}^v, used to retrieves f evaluations
pub fn n_to_vec(i: usize, n: usize) -> Vec<bool> {
	let x = format!("{:0>width$}", format!("{:b}", i), width = n);
	// println!("x: {:?}", x);
	let x: Vec<bool> = x.chars().map(|x| x == '1')
		.collect();
	// println!("x: {:?}, i: {:?}, n: {:?}", x, i, n);
	x
}

/// Perform a naive n^2 multiplication of `self` by `other`.
fn naive_mul(
	cur: &SparsePolynomial<F5, SparseTerm>,
	other: &SparsePolynomial<F5, SparseTerm>,
) -> SparsePolynomial<F5, SparseTerm> {
	if cur.is_zero() || other.is_zero() {
		// println!("zero");
		SparsePolynomial::zero()
	} else {
		let mut result_terms = Vec::new();
		for (cur_coeff, cur_term) in cur.terms.iter() {
			for (other_coeff, other_term) in other.terms.iter() {
				let mut term: Vec<(usize, usize)> = cur_term.vars().iter().zip(cur_term.powers()).map(|(v, p)| (*v, p)).collect();
				term.extend(
					other_term.vars().iter().zip(other_term.powers()).map(|(v, p)| (*v, p))
				);
				// let mut term = cur_term.0.clone();
				// term.extend(other_term.0.clone());
				let coeff = *cur_coeff * *other_coeff;
				// println!("naive mul coeff: {}, cur_coeff: {}, other_coeff: {}", coeff, cur_coeff, other_coeff);
				result_terms.push((coeff, SparseTerm::new(term)));
			}
		}
		let maxvars = core::cmp::max(cur.num_vars, other.num_vars);
		SparsePolynomial::from_coefficients_slice(maxvars, result_terms.as_slice())
	}
}

// Computes Chi_w(r) for all w, O(log n) operations
pub fn chi_w(w: &Vec<bool>, r: &Vec<i128>) -> i128 {
	// println!("chi w...");
	assert_eq!(w.len(), r.len());
	let product: SparsePolynomial<F5, SparseTerm> = w
		.iter()
		.enumerate()
		.fold(
			SparsePolynomial::from_coefficients_slice(0, &[(F5::from(1), SparseTerm::new(Vec::new()))]), 
			|p, (i, w)| {
				let step = chi_step(*w, i);
				let m = naive_mul(&p, &step);
				m
			}
		);

	let r: Vec<F5> = r.iter().map(|i| F5::from(*i)).collect();

	let s = &product.evaluate(&r);
	// println!("{:?}", product.terms()[0].0.into_repr());
	// println!("{:?}", product.terms()[1].0.into_repr());
	// println!("{:?}", product.terms()[2].0.into_repr());
	// println!("{:?}", product.terms()[3].0.into_repr());
	let result = i128::from_str_radix(&s.into_repr().to_string(), 16).unwrap();
	// println!("chi w result:{:?}", result);
	result
}

// Calculating the slow way, for benchmarking
pub fn slow_mle(fw: &Vec<i128>, r: &Vec<i128>, p: i128) -> i128 {
	assert_eq!(r.len() as f64, (fw.len() as f64).sqrt());
	let sum: i128 = fw
		.iter()
		.enumerate()
		.map(|(i, val)| val * chi_w(&n_to_vec(i, r.len()), r))
		.sum();
	
	// println!("slow mle sum {}", sum);
	// println!("sum%p {}", sum % p);
	sum % p
}

// Lemma 3.7
pub fn stream_mle(fw: &Vec<i128>, r: &Vec<i128>, p: i128) -> i128 {
	recurse(fw, r, 2usize.pow(r.len() as u32)) % p
}

pub fn recurse(fw: &Vec<i128>, r: &Vec<i128>, n: usize) -> i128 {
	match n {
		0 => 0,
		_ => recurse(fw, r, n - 1) + fw[n - 1] * chi_w(&n_to_vec(n - 1, r.len()), r),
	}
}

// // Lemm 3.8
// pub fn dynamic_mle(fw: &Vec<i128>, r: &Vec<i128>, p: i128) -> i128 {
// 	let chi_lookup = memoize(r, r.len());
// 	let result: i128 = fw
// 		.iter()
// 		.zip(chi_lookup.iter())
// 		.map(|(left, right)| left * right)
// 		.sum();
// 	result % p
// }

// pub fn memoize(r: &Vec<i128>, v: usize) -> Vec<i128> {
// 	match v {
// 		1 => {
// 			vec![chi_step(false, r[v - 1]), chi_step(true, r[v - 1])]
// 		}
// 		_ => memoize(r, v - 1)
// 			.iter()
// 			.flat_map(|val| {
// 				[
// 					val * chi_step(false, r[v - 1]),
// 					val * chi_step(true, r[v - 1]),
// 				]
// 			})
// 			.collect(),
// 	}
// }

pub fn convert_bin(x: usize, y: usize, n: usize) -> Vec<u32> {
	let xbin = format!("{:0>width$}", format!("{:b}", x), width = n);
	let ybin = format!("{:0>width$}", format!("{:b}", y), width = n);
	let bin = format!("{}{}", xbin, ybin);
	println!("{}", bin);
	// // println!("x: {:?}", x);
	let x: Vec<u32> = bin.chars().map(|x| x.to_digit(10).unwrap())
		.collect();
	x
}

pub fn convert_bin_vec (bin: Vec<u32>) -> Vec<i128> {
	bin.iter()
		.map(|i| i128::from_str_radix(&i.to_string(), 10).unwrap())
		.collect()
}

pub fn count_triangles(matrix: &Vec<i128>) -> u32 {
	let a = matrix.clone();
	// let b = matrix.clone();
	// let c = matrix.clone();

	let len: usize = (matrix.len() as f32).sqrt() as usize;

	let modulus = 5;
	let mut total_triangles = 0;

	for x in 0..len {
		for y in 0..len {
			let bin_xy = convert_bin(x, y, len/2);
			let exist_xy = slow_mle(&a, &convert_bin_vec(bin_xy), modulus);
			for z in 0..len {
				let bin_yz = convert_bin(y, z, len/2);
				let exist_yz = slow_mle(&a, &convert_bin_vec(bin_yz), modulus);

				let bin_xz = convert_bin(x, z, len/2);
				let exist_xz = slow_mle(&a, &convert_bin_vec(bin_xz), modulus);

				let exist = exist_xy * exist_yz * exist_xz;
				if exist != 0 {
					println!("exist {} at x: {}, y: {}, z: {}", exist, x, y, z);
				}
				total_triangles += exist;
				// if z != x && z != y {
				// }
			}
		}
	}

	total_triangles as u32 / 6

}