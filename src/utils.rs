use ark_ff::{Zero, One};
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