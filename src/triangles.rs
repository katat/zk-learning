use std::{iter};

use ark_ff::{Field};
use ark_poly::{multivariate::{SparseTerm, Term}};

use crate::{lagrange::{MultilinearExtension}, sumcheck::{SumCheckPolynomial, UniPoly}};

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

pub fn gen_var_indexes (start_index: usize, var_num: usize) -> Vec<usize> {
    let arr: Vec<usize> = (0..var_num).map(|x| x + start_index).collect();
    arr
}

pub fn convert_bin(x: usize, y: usize, n: usize) -> Vec<u32> {
    let xbin = format!("{:0>width$}", format!("{:b}", x), width = n);
    let ybin = format!("{:0>width$}", format!("{:b}", y), width = n);
    let bin = format!("{}{}", xbin, ybin);
    let x: Vec<u32> = bin.chars().map(|x| x.to_digit(10).unwrap())
        .collect();
    x
}

//TODO refactor for more efficient
pub fn convert_bin_z(x: usize, y: usize, z: usize, n: usize) -> Vec<u32> {
    let xbin = format!("{:0>width$}", format!("{:b}", x), width = n);
    let ybin = format!("{:0>width$}", format!("{:b}", y), width = n);
    let zbin = format!("{:0>width$}", format!("{:b}", z), width = n);
    let bin = format!("{}{}{}", xbin, ybin, zbin);
    let x: Vec<u32> = bin.chars().map(|x| x.to_digit(10).unwrap())
        .collect();
    x
}

impl <F: Field, E: MultilinearExtension<F>> TriangleMLE <F, E> {
    pub fn new(matrix: TriangleGraph<F>) -> Self {
        let a: Vec<F> = matrix.flatten();
        let var_num = matrix.one_dimension_size();
    
        let x_start_index = 0;
        let y_start_index = var_num;
        let z_start_index = var_num * 2;
        
        let x_indexes = gen_var_indexes(x_start_index, var_num);
        let y_indexes = gen_var_indexes(y_start_index, var_num);
        let mut xy_indexes: Vec<usize> = x_indexes.clone();
        xy_indexes.append(&mut y_indexes.clone());
    
        let mut z_indexes = gen_var_indexes(z_start_index, var_num);
    
        let mut yz_indexes: Vec<usize> = y_indexes;
        yz_indexes.append(&mut z_indexes.clone());
        
        let mut xz_indexes: Vec<usize> = x_indexes;
        xz_indexes.append(&mut z_indexes);

        let mle = E::new(a);

        TriangleMLE {
            matrix,
            f_xy: mle.clone(),
            f_yz: mle.clone(),
            f_xz: mle.clone(),
        }
    }
}

fn mul<F: Field>(cur: Vec<(F, SparseTerm)>, other: Vec<(F, SparseTerm)>) -> Vec<(F, SparseTerm)> {
    let mut result_terms = Vec::new();
    for (cur_coeff, cur_term) in cur.iter() {
        for (other_coeff, other_term) in other.iter() {
            let mut term: Vec<(usize, usize)> = cur_term.vars().iter().zip(cur_term.powers()).map(|(v, p)| (*v, p)).collect();
            term.extend(
                other_term.vars().iter().zip(other_term.powers()).map(|(v, p)| (*v, p))
            );
            let coeff = *cur_coeff * *other_coeff;
            result_terms.push((coeff, SparseTerm::new(term)));
        }
    }
    result_terms
}

impl <F: Field, E: MultilinearExtension<F>> SumCheckPolynomial<F> for TriangleMLE<F, E> {
    fn var_fixed_evaluate(&self, var: usize, point: Vec<F>) -> UniPoly<F> {
        // println!("\x1b[93mpoint {:?}\x1b[0m", point);

        let dim = self.matrix.one_dimension_size();
        let y_start_index = dim;
        let z_start_index = dim * 2;
        // prepend r_vec
        let point_xy = &point[0..(z_start_index)];
        let point_yz = &point[y_start_index..(z_start_index + self.matrix.one_dimension_size())];
        let point_xz = [&point[0..(y_start_index)], &point[z_start_index..]].concat();

        // println!("xy point {:?}", point_xy);
        // println!("yz point {:?}", point_yz);
        // println!("xz point {:?}", point_xz);

        let mut xy_mle = E::new(self.matrix.flatten());
        let mut yz_mle = E::new(self.matrix.flatten());
        let mut xz_mle = E::new(self.matrix.flatten());

        let mut xy_var= vec![];
        let mut yz_var= vec![];
        let mut xz_var= vec![];
        if var < z_start_index {
            xy_var.push(var);
            if var < y_start_index {
                xz_var.push(var);
            }
            else {
                yz_var.push(var - dim);
            }
        }
        else {
            yz_var.push(var - dim);
            xz_var.push(var - dim);
        }

        // println!("xy var {:?}", xy_var);
        // println!("yz var {:?}", yz_var);
        // println!("xz var {:?}", xz_var);

        xy_mle.fix_vars(&xy_var, point_xy.to_vec());
        yz_mle.fix_vars(&yz_var, point_yz.to_vec());
        xz_mle.fix_vars(&xz_var, point_xz.to_vec());

        // println!("xy mle fixed vars, evals {:?}", xy_mle.to_evals());
        // println!("yz mle fixed vars, evals {:?}", yz_mle.to_evals());
        // println!("xz mle fixed vars, evals {:?}", xz_mle.to_evals());

        let xy_poly = xy_mle.interpolate();
        let yz_poly = yz_mle.interpolate();
        let xz_poly = xz_mle.interpolate();

        // println!("xy poly {:?}", xy_poly);
        // println!("yz poly {:?}", yz_poly);
        // println!("xz poly {:?}", xz_poly);

        xy_poly.mul(&yz_poly).mul(&xz_poly)

        // let f_xy_unipoly: UniPoly = self.f_xy.terms.clone().into_iter().fold(
		// 	UniPoly::from_coefficients_vec(vec![]),
		// 	|sum, term| {
		// 		let curr = cb(term);
		// 		sum + curr
		// 	}
		// );
        // let f_yz_unipoly: UniPoly = self.f_yz.terms.clone().into_iter().fold(
		// 	UniPoly::from_coefficients_vec(vec![]),
		// 	|sum, term| {
		// 		let curr = cb(term);
		// 		sum + curr
		// 	}
		// );
        // let f_xz_unipoly: UniPoly = self.f_xz.terms.clone().into_iter().fold(
		// 	UniPoly::from_coefficients_vec(vec![]),
		// 	|sum, term| {
		// 		let curr = cb(term);
		// 		sum + curr
		// 	}
		// );
        // f_xy_unipoly.mul(&f_yz_unipoly).mul(&f_xz_unipoly)
    }

    fn num_vars(&self) -> usize {
        (self.f_xy.num_vars() + self.f_yz.num_vars() + self.f_xz.num_vars()) / 2
    }

    fn evaluate(&self, point: &Vec<F>) -> F {
        let y_start_index = self.matrix.one_dimension_size();
        let z_start_index = self.matrix.one_dimension_size() * 2;
        let point_xy = &point[0..(z_start_index)];
        let point_yz = &point[y_start_index..(z_start_index + self.matrix.one_dimension_size())];
        let point_xz = [&point[0..(y_start_index)], &point[z_start_index..]].concat();

        // let xy_eval = eval_slow_mle(&self.matrix.flatten(), &point_xy.to_vec());
        // let yz_eval = eval_slow_mle(&self.matrix.flatten(), &point_yz.to_vec());
        // let xz_eval = eval_slow_mle(&self.matrix.flatten(), &point_xz.to_vec());
        // let xyp = point_xy.iter().map(|e| poly_constant(*e)).collect();
        // let yzp = point_yz.iter().map(|e| poly_constant(*e)).collect();
        // let xzp = point_xz.iter().map(|e| poly_constant(*e)).collect();
        let mle = E::new(self.matrix.flatten());
        let xy_eval = mle.evaluate(&point_xy.to_vec());
        let yz_eval = mle.evaluate(&point_yz.to_vec());
        let xz_eval = mle.evaluate(&point_xz);

        xy_eval * yz_eval * xz_eval

        // match self.eval_type {
        //     MLEAlgorithm::SlowMLE => {
        //         let xy_eval = mle.slow_eval(&point_xy.to_vec());
        //         let yz_eval = mle.slow_eval(&point_yz.to_vec());
        //         let xz_eval = mle.slow_eval(&point_xz);
        
        //         xy_eval * yz_eval * xz_eval
        //     }
        //     MLEAlgorithm::DynamicMLE => {
        //         let xy_eval = mle.dynamic_eval(&point_xy.to_vec());
        //         let yz_eval = mle.dynamic_eval(&point_yz.to_vec());
        //         let xz_eval = mle.dynamic_eval(&point_xz);
        
        //         xy_eval * yz_eval * xz_eval
        //     }
        //     MLEAlgorithm::StreamMLE => {
        //         let xy_eval = mle.stream_eval(&point_xy.to_vec());
        //         let yz_eval = mle.stream_eval(&point_yz.to_vec());
        //         let xz_eval = mle.stream_eval(&point_xz);
        
        //         xy_eval * yz_eval * xz_eval
        //     },
        // }
    }
}
