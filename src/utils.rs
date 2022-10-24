use ark_ff::{Zero, One, Field};
use crate::small_fields::F251;

pub fn gen_matrix (len: usize) -> Vec<Vec<F251>> {
	let mut matrix = Vec::new();
    for i in 0..len {
		let mut row: Vec<F251> = Vec::new();
		for j in 0..len {
			if i == j {
				row.push(F251::zero());
				continue;
			}
			row.push(F251::one());
		}
		matrix.push(row)
    }

	matrix[0][len - 1] = F251::zero();
	matrix[len - 1][0] = F251::zero();
	matrix
}

pub fn convert_field<F: Field>(vec: &[u32]) -> Vec<F> {
	vec.iter().map(|e| F::from(*e)).collect()
}

// Converts i into an index in {0,1}^v, used to retrieves f evaluations
pub fn n_to_vec(i: usize, n: usize) -> Vec<bool> {
	let x = format!("{:0>width$}", format!("{:b}", i), width = n);
	let x: Vec<bool> = x.chars().map(|x| x == '1')
		.collect();
	x
}