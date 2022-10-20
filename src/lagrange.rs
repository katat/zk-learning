use ark_ff::{Zero, One, PrimeField, Field};
use ark_poly::{multivariate::{SparsePolynomial, SparseTerm, Term}, Polynomial, DenseMVPolynomial};
use crate::{small_fields::F251};


// One step in chi
pub fn poly_chi_step<F: Field>(w: bool, v: usize) -> SparsePolynomial<F, SparseTerm> {
	let one_minus_w = F::from(1 - (w as u32));
	let w_minus_one = F::from(w as u32) - F::one();

	let f: SparsePolynomial<F, SparseTerm> = SparsePolynomial::from_coefficients_vec(
		v + 1,
		vec![
			//x * i128::from(w)
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
// todo try to have this as * operator for SparsePolynomial
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
				// let mut term = cur_term.0.clone();
				// term.extend(other_term.0.clone());
				let coeff = *cur_coeff * *other_coeff;
				// //println!("naive mul coeff: {}, cur_coeff: {}, other_coeff: {}", coeff, cur_coeff, other_coeff);
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
				let factor = SparsePolynomial::from_coefficients_slice(0, &[(F::from(*val), SparseTerm::new(Vec::new()))]);
				sum + naive_mul(&factor, &poly_chi_w(&n_to_vec(i, vars.len()), vars))
			}
		);

	sum
}

pub fn eval_chi_step<F: Field>(w: bool, x: F) -> F {
	x * F::from(w) + (F::one() - x) * (F::one() - F::from(w))
}

pub fn eval_chi_w<F: Field>(w: &[bool], point: &Vec<F>) -> F {
	w.iter()
	.enumerate()
	.zip(point)
	.fold(
		F::one(),
		|prev, ((_, w), v)| {
			let step = eval_chi_step(*w, *v);
			prev * step
		}
	)
}

pub fn eval_slow_mle(fw: &[F251], point: &Vec<F251>) -> F251 {
	let sum: F251 = fw
		.iter()
		.enumerate()
		.fold(
			F251::from(0),
			|sum, (i, val)| {
				let factor = *val;
				sum + factor * eval_chi_w(&n_to_vec(i, point.len()), point)
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

pub fn eval_poly_chi_step<F: Field>(w: bool, x: &SparsePolynomial<F, SparseTerm>) -> SparsePolynomial<F, SparseTerm> {
	let one_minus_w = F::one() - F::from(w as u32);
	let w_minus_one = F::from(w as u32) - F::one();

	naive_mul(
		&poly_constant(F::from(w as u32)),
		x
	) +
	poly_constant(one_minus_w) +
	naive_mul(
		&poly_constant(w_minus_one),
		x
	)
}

pub fn eval_memoize(r: &Vec<F251>, v: usize) -> Vec<F251> {
	match v {
		1 => {
			vec![eval_chi_step(false, r[v - 1]), eval_chi_step(true, r[v - 1])]
		}
		_ => eval_memoize(r, v - 1)
			.iter()
			.flat_map(|val| {
				[
					*val * eval_chi_step(false, r[v - 1]),
					*val * eval_chi_step(true, r[v - 1]),
				]
			})
			.collect(),
	}
}

pub fn eval_dynamic_mle(fw: &Vec<F251>, r: &Vec<F251>) -> F251 {
	let chi_lookup = eval_memoize(r, r.len());
	fw.iter()
		.zip(chi_lookup.iter())
		.map(|(left, right)| *left * *right)
		.sum()
}

// Lemma 3.7
// pub fn stream_mle(fw: &Vec<i128>, r: &Vec<i128>, p: i128) -> i128 {
// 	recurse(fw, r, 2usize.pow(r.len() as u32)) % p
// }

// pub fn recurse(fw: &Vec<i128>, r: &Vec<i128>, n: usize) -> i128 {
// 	match n {
// 		0 => 0,
// 		_ => recurse(fw, r, n - 1) + fw[n - 1] * chi_w(&n_to_vec(n - 1, r.len()), r),
// 	}
// }

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
	//println!("{}", bin);
	// // //println!("x: {:?}", x);
	let x: Vec<u32> = bin.chars().map(|x| x.to_digit(10).unwrap())
		.collect();
	x
}

pub fn convert_bin_z(x: usize, y: usize, z: usize, n: usize) -> Vec<u32> {
	let xbin = format!("{:0>width$}", format!("{:b}", x), width = n);
	let ybin = format!("{:0>width$}", format!("{:b}", y), width = n);
	let zbin = format!("{:0>width$}", format!("{:b}", z), width = n);
	let bin = format!("{}{}{}", xbin, ybin, zbin);
	//println!("{}", bin);
	// // //println!("x: {:?}", x);
	let x: Vec<u32> = bin.chars().map(|x| x.to_digit(10).unwrap())
		.collect();
	x
}

pub fn convert_bin_vec (bin: Vec<u32>) -> Vec<u32> {
	bin.iter()
		.map(|i| u32::from_str_radix(&i.to_string(), 10).unwrap())
		.collect()
}

pub fn gen_var_indexes (start_index: usize, var_num: usize) -> Vec<usize> {
	let arr: Vec<usize> = (0..var_num).map(|x| x + start_index).collect();
	arr
}

pub fn poly_count_triangles<F: Field>(matrix: &Vec<F>) -> SparsePolynomial<F, SparseTerm> {
	let a = matrix.clone();

	let len: usize = (matrix.len() as f32).sqrt() as usize;


	// todo first get polynomial resulted from multiplied matrix MLEs, then evaluate over the coordinates
	// /todo in every sum loop, new variables such as x, should be added, such as adding x3 to [x1, x2]
	let var_num = (len as f32).log2() as usize;

	let x_indexes = gen_var_indexes(0, var_num);
	let y_indexes = gen_var_indexes(x_indexes.last().unwrap() + 1, var_num);
	//println!("x indexes {:?}", x_indexes);
	//println!("y indexes {:?}", y_indexes);
	let mut xy_indexes: Vec<usize> = x_indexes.clone();
	xy_indexes.append(&mut y_indexes.clone());
	//println!("xy indexes {:?}", xy_indexes);

	let z_indexes = gen_var_indexes(y_indexes.last().unwrap() + 1, var_num);
	//println!("z indexes {:?}", z_indexes);

	let mut yz_indexes: Vec<usize> = y_indexes.clone();
	yz_indexes.append(&mut z_indexes.clone());
	//println!("yz indexes {:?}", yz_indexes);
	
	let mut xz_indexes: Vec<usize> = x_indexes.clone();
	xz_indexes.append(&mut z_indexes.clone());
	//println!("xz indexes {:?}", xz_indexes);
	let poly_exist_xy = poly_slow_mle(&a, &xy_indexes);
	let poly_exist_yz = poly_slow_mle(&a, &yz_indexes);
	let poly_exist_xz = poly_slow_mle(&a, &xz_indexes);
	let poly_exist = naive_mul(&naive_mul(&poly_exist_xy, &poly_exist_yz), &poly_exist_xz);
	poly_exist
}

pub fn count_triangles<F: Field>(matrix: &Vec<F>) -> F {
	let poly_exist = poly_count_triangles(matrix);
	let len: usize = (matrix.len() as f32).sqrt() as usize;
	let var_num = (len as f32).log2() as usize;
	let mut total_triangles = F::zero();

	for x in 0..len {
		for y in 0..len {
			for z in 0..len {
				let xyz_bin = convert_bin_vec(convert_bin_z(x, y, z, var_num));
				let r: Vec<F> = xyz_bin.iter().map(|i| F::from(*i)).collect();

				let result: F = poly_exist.evaluate(&r);

				// if exist != 0 {
				// 	println!("exist {} at x: {}, y: {}, z: {}", exist, x, y, z);
				// }
				total_triangles += result;
			}
		}
	}

	// todo get polynomial multiplied before evaluation
	total_triangles.div(F::from(6u32))

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
	assert_eq!(r.len() as f64, (fw.len() as f64).sqrt());
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
pub fn dynamic_mle<F: Field>(fw: &Vec<F>, r: &Vec<F>) -> F {
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
