use std::{iter};
use ark_ff::{Field};
use crate::{lagrange::{MultilinearExtension, UniPoly}, sumcheck::{SumCheckPolynomial, n_to_vec}};

#[derive(Debug, Clone)]
pub struct TriangleGraph <F: Field>  {
    vec: Vec<Vec<F>>
}

impl <F: Field> TriangleGraph<F> {
    pub fn new(m: Vec<Vec<F>>) -> Self {
        TriangleGraph {
            vec: m
        }
    }
    pub fn flatten(&self) -> Vec<F>{
        self.vec.iter().flatten().cloned().collect()
    }

    pub fn size(&self) -> usize{
        (self.flatten().len() as f32).sqrt() as usize
    }

    pub fn one_dimension_size(&self) -> usize {
        (self.size() as f32).log2() as usize
    }

    pub fn get(&self, x: usize, y: usize) -> F {
        self.vec[x][y]
    }

    pub fn multiply(&self, matrix: TriangleGraph<F>) -> TriangleGraph<F> {
        let size = self.size();
        let row: Vec<F> = iter::repeat(F::zero()).take(size).collect();
        let mut result_matrix: Vec<Vec<F>> = iter::repeat(row).take(size).collect();
        for i in 0..size {
            for j in 0..size {
                let mut elm = F::zero();
                for k in 0..size {
                    elm += self.get(i, k) * matrix.get(k, j);
                }
                result_matrix[i][j] = elm;
            }
        }

        TriangleGraph {
            vec: result_matrix
        }
    }

    pub fn count(&self) -> F {
        let a2 = self.multiply(self.clone());
        let a3 = self.multiply(a2);
        let mut number = F::zero();

        for i in 0..a3.size() {
            number += a3.get(i, i);
        }

        number.div(F::from(6u32))
    }

    pub fn derive_mle <E: MultilinearExtension<F>>(&self) -> TriangleMLE<F, E> {
        TriangleMLE::new(self.clone())
    }
}

#[derive(Debug, Clone)]
pub struct TriangleMLE <F: Field, E: MultilinearExtension<F>> {
	pub matrix: TriangleGraph<F>,
    f_xy: E,
    f_yz: E,
    f_xz: E,
}

impl <F: Field, E: MultilinearExtension<F>> TriangleMLE <F, E> {
    pub fn new(matrix: TriangleGraph<F>) -> Self {
        let a: Vec<F> = matrix.flatten();
        let var_num = matrix.one_dimension_size();
    
        let x_start_index = 0;
        let y_start_index = var_num;
        let z_start_index = var_num * 2;
        
        let x_indexes = Self::gen_var_indexes(x_start_index, var_num);
        let y_indexes = Self::gen_var_indexes(y_start_index, var_num);
        let mut xy_indexes: Vec<usize> = x_indexes.clone();
        xy_indexes.append(&mut y_indexes.clone());
    
        let mut z_indexes = Self::gen_var_indexes(z_start_index, var_num);
    
        let mut yz_indexes: Vec<usize> = y_indexes;
        yz_indexes.append(&mut z_indexes.clone());
        
        let mut xz_indexes: Vec<usize> = x_indexes;
        xz_indexes.append(&mut z_indexes);

        let xy_mle = E::new(a.clone(), Option::Some(xy_indexes.clone()));
        let yz_mle = E::new(a.clone(), Option::Some(yz_indexes.clone()));
        let xz_mle = E::new(a, Option::Some(xz_indexes));

        TriangleMLE {
            matrix,
            f_xy: xy_mle,
            f_yz: yz_mle,
            f_xz: xz_mle,
        }
    }

    fn gen_var_indexes (start_index: usize, var_num: usize) -> Vec<usize> {
        let arr: Vec<usize> = (0..var_num).map(|x| x + start_index).collect();
        arr
    }

    // Sum all evaluations of polynomial `g` over boolean hypercube
	pub fn hypercube_sum(&self) -> F {
		let v = self.num_vars();
		let n = 2u32.pow(v as u32);
		(0..n)
			.map(|n| self.evaluate(&n_to_vec(n as usize, v)))
			.sum()
	}
}

impl <F: Field, E: MultilinearExtension<F>> SumCheckPolynomial<F> for TriangleMLE<F, E> {
    fn var_fixed_evaluate(&self, var: usize, point: Vec<F>) -> UniPoly<F> {
        let mut f_xy = self.f_xy.clone();
        let mut f_yz = self.f_yz.clone();
        let mut f_xz = self.f_xz.clone();

        f_xy.fix_vars(&[var], point.to_vec());
        f_yz.fix_vars(&[var], point.to_vec());
        f_xz.fix_vars(&[var], point.to_vec());

        let xy_poly = f_xy.interpolate();
        let yz_poly = f_yz.interpolate();
        let xz_poly = f_xz.interpolate();

        xy_poly.mul(&yz_poly).mul(&xz_poly)
    }

    fn num_vars(&self) -> usize {
        (self.f_xy.num_vars() + self.f_yz.num_vars() + self.f_xz.num_vars()) / 2
    }

    fn evaluate(&self, point: &Vec<F>) -> F {
        let xy_eval = self.f_xy.evaluate(point);
        let yz_eval = self.f_yz.evaluate(point);
        let xz_eval = self.f_xz.evaluate(point);

        xy_eval * yz_eval * xz_eval
    }
}
