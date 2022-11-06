use std::marker::PhantomData;
use ark_ff::{Field};
use crate::{lagrange::{MultilinearExtension, EvaluationMethod, UniPoly}};

#[derive(Debug, Clone)]
pub struct ValueBasedMultilinearExtension<F: Field, E: EvaluationMethod<F>> {
	evals: Vec<F>,
    indexes: Vec<usize>,
	phantom: PhantomData<E>
}

impl <F: Field, E: EvaluationMethod<F>> MultilinearExtension<F> for ValueBasedMultilinearExtension<F, E> {
	fn new(evals: Vec<F>, indexes: Option<Vec<usize>>) -> Self {
		ValueBasedMultilinearExtension {
			evals,
			indexes: indexes.unwrap(),
			phantom: PhantomData,
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
						
						point.splice(loc..loc+1, vec![b].iter().cloned());
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
        E::run(&self.evals, &self.indexes, point)
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