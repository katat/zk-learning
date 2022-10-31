use std::{iter};
use ark_ff::{Field};
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

        let xy_mle = E::new(a.clone(), Option::Some(xy_indexes.clone()));
        let yz_mle = E::new(a.clone(), Option::Some(yz_indexes.clone()));
        let xz_mle = E::new(a.clone(), Option::Some(xz_indexes));

        TriangleMLE {
            matrix,
            f_xy: xy_mle.clone(),
            f_yz: yz_mle.clone(),
            f_xz: xz_mle.clone(),
        }
    }
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

        let mut xy_var= vec![];
        let mut yz_var= vec![];
        let mut xz_var= vec![];
        if var < z_start_index {
            xy_var.push(var);
            if var < y_start_index {
                xz_var.push(var);
            }
            else {
                yz_var.push(var);
            }
        }
        else {
            yz_var.push(var);
            xz_var.push(var);
        }

        // println!("xy var {:?}", xy_var);
        // println!("yz var {:?}", yz_var);
        // println!("xz var {:?}", xz_var);

        let mut f_xy = self.f_xy.clone();
        let mut f_yz = self.f_yz.clone();
        let mut f_xz = self.f_xz.clone();

        println!("f xy var");
        f_xy.fix_vars(&xy_var, point_xy.to_vec());
        println!("f yz var");
        f_yz.fix_vars(&yz_var, point_yz.to_vec());
        println!("f xz var");
        f_xz.fix_vars(&xz_var, point_xz.to_vec());

        // println!("xy mle fixed vars, evals {:?}", xy_mle.to_evals());
        // println!("yz mle fixed vars, evals {:?}", yz_mle.to_evals());
        // println!("xz mle fixed vars, evals {:?}", xz_mle.to_evals());

        let xy_poly = f_xy.interpolate();
        let yz_poly = f_yz.interpolate();
        let xz_poly = f_xz.interpolate();

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

        let xy_eval = self.f_xy.evaluate(&point_xy.to_vec());
        let yz_eval = self.f_yz.evaluate(&point_yz.to_vec());
        let xz_eval = self.f_xz.evaluate(&point_xz);

        xy_eval * yz_eval * xz_eval
    }
}
