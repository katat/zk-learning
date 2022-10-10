use ark_ff::{Zero, PrimeField};
use ark_poly::{multivariate::{SparsePolynomial, SparseTerm, Term}, Polynomial, DenseMVPolynomial};
use crate::small_fields::F251;


// One step in chi
pub fn chi_step(w: bool, v: usize) -> SparsePolynomial<F251, SparseTerm> {
	// x * i128::from(w) + (1 - x) * (1 - i128::from(w));
	// //println!("step v+1:{}", v+1);

	let one_minus_w = (1 - (w as i128)).into();
	// //println!("1-w: {}", one_minus_w);
	let w_minus_one = w as i128 - 1;
	// //println!("w-1: {}", w_minus_one);

	let f: SparsePolynomial<F251, SparseTerm> = SparsePolynomial::from_coefficients_vec(
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

	// //println!("2nd term coeff: {:?}", f.terms()[1]);

	f
}

// Converts i into an index in {0,1}^v, used to retrieves f evaluations
pub fn n_to_vec(i: usize, n: usize) -> Vec<bool> {
	let x = format!("{:0>width$}", format!("{:b}", i), width = n);
	// //println!("x: {:?}", x);
	let x: Vec<bool> = x.chars().map(|x| x == '1')
		.collect();
	// //println!("x: {:?}, i: {:?}, n: {:?}", x, i, n);
	x
}

/// Perform a naive n^2 multiplication of `self` by `other`.
// todo try to have this as * operator for SparsePolynomial
pub fn naive_mul(
	cur: &SparsePolynomial<F251, SparseTerm>,
	other: &SparsePolynomial<F251, SparseTerm>,
) -> SparsePolynomial<F251, SparseTerm> {
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
pub fn chi_w(w: &[bool], vars: &Vec<usize>) -> SparsePolynomial<F251, SparseTerm> {
	// //println!("chi w...");
	let product: SparsePolynomial<F251, SparseTerm> = w
		.iter()
		.enumerate()
		.zip(vars)
		.fold(
			SparsePolynomial::from_coefficients_slice(0, &[(F251::from(1), SparseTerm::new(Vec::new()))]), 
			|p, ((_, w), v)| {
				let step = chi_step(*w, *v);
				naive_mul(&p, &step)
			}
		);

	product

	// let r: Vec<F251> = r.iter().map(|i| F251::from(*i)).collect();

	// lets = &product.evaluate(&r);
	// let result = i128::from_str_radix(&s.into_repr().to_string(), 16).unwrap();
	// result
}

// Calculating the slow way, for benchmarking
pub fn poly_slow_mle(fw: &[i128], vars: &Vec<usize>) -> SparsePolynomial<F251, SparseTerm> {
	let sum: SparsePolynomial<F251, SparseTerm> = fw
		.iter()
		.enumerate()
		.fold(
			SparsePolynomial::from_coefficients_slice(0, &[]),
			|sum, (i, val)| {
				let factor = SparsePolynomial::from_coefficients_slice(0, &[(F251::from(*val), SparseTerm::new(Vec::new()))]);
				// [x1, x2, y1, y2]
				sum + naive_mul(&factor, &chi_w(&n_to_vec(i, vars.len()), vars))
			}
		);

	sum
}

pub fn slow_mle(fw: &Vec<i128>, r: &Vec<i128>, vars: &Vec<usize>) -> F251 {
	let fw_len = (fw.len() as f64).sqrt();
	assert_eq!(r.len() as f64, fw_len);
	let sum: SparsePolynomial<F251, SparseTerm> = poly_slow_mle(fw, vars);
	
	let mut r: Vec<F251> = r.iter().map(|i| F251::from(*i)).collect();

	let diff = sum.num_vars() - r.len();
	r.splice::<_, _>(0..0, std::iter::repeat(F251::from(0i32)).take(diff));
	let s: F251 = sum.evaluate(&r);
	s
	// let result = i128::from_str_radix(&s.into_repr().to_string(), 16).unwrap();
	// result 
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

pub fn convert_bin_vec (bin: Vec<u32>) -> Vec<i128> {
	bin.iter()
		.map(|i| i128::from_str_radix(&i.to_string(), 10).unwrap())
		.collect()
}

pub fn gen_var_indexes (start_index: usize, var_num: usize) -> Vec<usize> {
	let arr: Vec<usize> = (0..var_num).map(|x| x + start_index).collect();
	arr
}

pub fn poly_count_triangles(matrix: &Vec<i128>) -> SparsePolynomial<F251, SparseTerm> {
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

pub fn count_triangles(matrix: &Vec<i128>) -> u32 {
	let poly_exist = poly_count_triangles(matrix);
	let len: usize = (matrix.len() as f32).sqrt() as usize;
	let var_num = (len as f32).log2() as usize;
	let mut total_triangles = 0;

	for x in 0..len {
		for y in 0..len {
			for z in 0..len {
				let xyz_bin = convert_bin_vec(convert_bin_z(x, y, z, var_num));
				let r: Vec<F251> = xyz_bin.iter().map(|i| F251::from(*i)).collect();

				let result: F251 = poly_exist.evaluate(&r);
				
				let exist = result.into_bigint().as_ref()[0];

				if exist != 0 {
					println!("exist {} at x: {}, y: {}, z: {}", exist, x, y, z);
				}
				total_triangles += exist;
			}
		}
	}

	// todo get polynomial multiplied before evaluation
	total_triangles as u32 / 6

}