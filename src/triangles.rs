use std::{iter};

use ark_ff::{Field};
use ark_poly::{multivariate::{SparsePolynomial, SparseTerm, Term}};

use crate::{lagrange::{poly_slow_mle, naive_mul, dynamic_mle}, sumcheck::{SumCheckPolynomial, UniPoly}};

#[derive(Debug, Clone)]
pub struct Matrix <F: Field>  {
    vec: Vec<Vec<F>>
}

#[derive(Debug, Clone)]
pub struct Triangles <F: Field> {
	pub matrix: Matrix<F>,
    f_xy: SparsePolynomial<F, SparseTerm>,
    f_yz: SparsePolynomial<F, SparseTerm>,
    f_xz: SparsePolynomial<F, SparseTerm>,
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

impl <F: Field> Matrix<F> {
    pub fn new(m: Vec<Vec<F>>) -> Self {
        Matrix {
            vec: m
        }
    }
    pub fn flatten(&self) -> Vec<F>{
        self.vec.iter().flatten().cloned().collect()
    }

    pub fn size(&self) -> usize{
        (self.flatten().len() as f32).sqrt() as usize
    }

    pub fn var_num(&self) -> usize {
        (self.size() as f32).log2() as usize
    }

    pub fn get(&self, x: usize, y: usize) -> F {
        self.vec[x][y]
    }
}

impl <F: Field> Triangles <F> {
    pub fn new(matrix: Vec<Vec<F>>) -> Self {
        let _matrix = Matrix { vec: matrix };
        let a: Vec<F> = _matrix.flatten();
        let var_num = _matrix.var_num();
    
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

        //optimize these polynomial representations
        let poly_exist_xy = poly_slow_mle(&a, &xy_indexes);
        let poly_exist_yz = poly_slow_mle(&a, &yz_indexes);
        let poly_exist_xz = poly_slow_mle(&a, &xz_indexes);

        Triangles {
            matrix: _matrix,
            f_xy: poly_exist_xy,
            f_yz: poly_exist_yz,
            f_xz: poly_exist_xz,
        }
    }

    pub fn multiply(&self, matrix: Matrix<F>) -> Matrix<F> {
        let size = self.matrix.size();
        let row: Vec<F> = iter::repeat(F::zero()).take(size).collect();
        let mut result_matrix: Vec<Vec<F>> = iter::repeat(row).take(size).collect();
        for i in 0..size {
            for j in 0..size {
                let mut elm = F::zero();
                for k in 0..size {
                    elm += self.matrix.get(i, k) * matrix.get(k, j);
                }
                result_matrix[i][j] = elm;
            }
        }

        Matrix {
            vec: result_matrix
        }
    }

    pub fn count(&self) -> F {
        let a2 = self.multiply(self.matrix.clone());
        let a3 = self.multiply(a2);
        let mut number = F::zero();

        for i in 0..a3.size() {
            number += a3.get(i, i);
        }

        number.div(F::from(6u32))
    }

    pub fn poly_count_triangles(&self) -> SparsePolynomial<F, SparseTerm> {
        let a: Vec<F> = self.matrix.flatten();
    
        let var_num = self.matrix.var_num();
    
        // encapsulate them as first/second/third...
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

        //clean up
        
        let poly_exist_xy = poly_slow_mle(&a, &xy_indexes);
        let poly_exist_yz = poly_slow_mle(&a, &yz_indexes);
        let poly_exist_xz = poly_slow_mle(&a, &xz_indexes);

        println!("poly xz {}", poly_exist_xz.terms.len());
        naive_mul(&naive_mul(&poly_exist_xy, &poly_exist_yz), &poly_exist_xz)
    }

    pub fn count_by_mle(&self) -> F {
        let poly_exist = self.poly_count_triangles();
        let len = self.matrix.size();
        let var_num = self.matrix.var_num();
        let mut total_triangles = F::zero();
    
        for x in 0..len {
            for y in 0..len {
                for z in 0..len {
                    let xyz_bin = convert_bin_z(x, y, z, var_num);
                    let r: Vec<F> = xyz_bin.iter().map(|i| F::from(*i)).collect();
    
                    let result = ark_poly::Polynomial::evaluate(&poly_exist, &r);
                    total_triangles += result;
                }
            }
        }
    
        total_triangles.div(F::from(6u32))
    
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

impl <F: Field> SumCheckPolynomial<F> for Triangles<F> {
    fn terms(&self) -> Vec<(F, SparseTerm)> {
        mul(mul(self.f_xy.terms.clone(), self.f_yz.terms.clone()), self.f_xz.terms.clone())
    }

    fn var_fixed_evaluate<C>(&self, mut cb: C) -> UniPoly where C: FnMut((F, SparseTerm)) -> UniPoly {

        let f_xy_unipoly: UniPoly = self.f_xy.terms.clone().into_iter().fold(
			UniPoly::from_coefficients_vec(vec![]),
			|sum, term| {
				let curr = cb(term);
				sum + curr
			}
		);
        let f_yz_unipoly: UniPoly = self.f_yz.terms.clone().into_iter().fold(
			UniPoly::from_coefficients_vec(vec![]),
			|sum, term| {
				let curr = cb(term);
				sum + curr
			}
		);
        let f_xz_unipoly: UniPoly = self.f_xz.terms.clone().into_iter().fold(
			UniPoly::from_coefficients_vec(vec![]),
			|sum, term| {
				let curr = cb(term);
				sum + curr
			}
		);
        f_xy_unipoly.mul(&f_yz_unipoly).mul(&f_xz_unipoly)
    }

    fn num_vars(&self) -> usize {
        *[
            self.f_xy.num_vars,
            self.f_yz.num_vars,
            self.f_xz.num_vars,
        ].iter().max().unwrap()
    }

    fn evaluate(&self, point: &Vec<F>) -> F {
        let y_start_index = self.matrix.var_num();
        let z_start_index = self.matrix.var_num() * 2;
        let point_xy = &point[0..(z_start_index)];
        let point_yz = &point[y_start_index..(z_start_index + self.matrix.var_num())];
        let point_xz = [&point[0..(y_start_index)], &point[z_start_index..]].concat();

        // println!("point {:?}", point);
        // println!("xy point {:?} {}", point_xy, self.f_xy.num_vars);
        // println!("f xy terms {}", self.f_xy.terms.len());
        // println!("point len {}", point.len());
        // let xy_evaluation = ark_poly::Polynomial::evaluate(&self.f_xy, point);
        // let yz_evaluation = ark_poly::Polynomial::evaluate(&self.f_yz, point);
        // let xz_evaluation = ark_poly::Polynomial::evaluate(&self.f_xz, point);

        // let xy_eval = eval_slow_mle(&self.matrix.flatten(), &point_xy.to_vec());
        // let yz_eval = eval_slow_mle(&self.matrix.flatten(), &point_yz.to_vec());
        // let xz_eval = eval_slow_mle(&self.matrix.flatten(), &point_xz.to_vec());
        // let xyp = point_xy.iter().map(|e| poly_constant(*e)).collect();
        // let yzp = point_yz.iter().map(|e| poly_constant(*e)).collect();
        // let xzp = point_xz.iter().map(|e| poly_constant(*e)).collect();
        let xy_eval = dynamic_mle(&self.matrix.flatten(), &point_xy.to_vec());
        let yz_eval = dynamic_mle(&self.matrix.flatten(), &point_yz.to_vec());
        let xz_eval = dynamic_mle(&self.matrix.flatten(), &point_xz);

        xy_eval * yz_eval * xz_eval
        // ark_poly::Polynomial::evaluate(&xy_eval, &vec![]) * 
        // ark_poly::Polynomial::evaluate(&yz_eval, &vec![]) * 
        // ark_poly::Polynomial::evaluate(&xz_eval, &vec![]) 

        // xy_evaluation * yz_evaluation * xz_evaluation
        // xy_evaluation
    }
}
