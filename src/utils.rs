use ark_ff::{Field};

pub fn gen_matrix <F: Field> (len: usize) -> Vec<Vec<F>> {
	let mut matrix = Vec::new();
    for i in 0..len {
		let mut row: Vec<F> = Vec::new();
		for j in 0..len {
			if i == j {
				row.push(F::zero());
				continue;
			}
			row.push(F::one());
		}
		matrix.push(row)
    }

	matrix[0][len - 1] = F::zero();
	matrix[len - 1][0] = F::zero();
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