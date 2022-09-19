use std::{iter, ops::Div};

use ark_ff::{Zero, PrimeField};
use ark_poly::multivariate::{SparsePolynomial, SparseTerm};

use crate::{small_fields::F251, lagrange::{poly_slow_mle, naive_mul}};

#[derive(Debug, Clone)]
pub struct Triangles {
	pub matrix: Vec<Vec<F251>>,
}

fn println_matrix(matrix: Vec<Vec<F251>>) {
    for i in 0..matrix.len() {
        for j in 0..matrix[i].len() {
            print!("{} ", matrix[i][j].into_repr());
        }
        println!("");
    }
}

impl Triangles {
    pub fn new(matrix: Vec<Vec<F251>>) -> Self {
        Triangles {
            matrix: matrix.clone()
        }
    }

    pub fn multiply(&self, matrix: Vec<Vec<F251>>) -> Vec<Vec<F251>> {
        let size = self.matrix.len();
        let row: Vec<F251> = iter::repeat(F251::zero()).take(size).collect();
        let mut result_matrix: Vec<Vec<F251>> = iter::repeat(row).take(size).collect();
        for i in 0..size {
            for j in 0..size {
                let mut elm = F251::zero();
                for k in 0..size {
                    elm += self.matrix[i][k] * matrix[k][j];
                    println!("{}", matrix[k][j]);
                }
                result_matrix[i][j] = elm;
            }
        }
        println!("matrix a:");
        println_matrix(self.matrix.clone());

        println!("matrix b:");
        println_matrix(matrix);

        println!("matrix c:");
        println_matrix(result_matrix.clone());

        result_matrix
    }

    pub fn count(&self) -> F251 {
        let a2 = self.multiply(self.matrix.clone());
        let a3 = self.multiply(a2);
        let mut number = F251::zero();

        for i in 0..a3.len() {
            number += a3[i][i];
        }

        number.div(F251::from(6))
    }

    pub fn gen_var_indexes (start_index: usize, var_num: usize) -> Vec<usize> {
        let arr: Vec<usize> = (0..var_num).map(|x| x + start_index).collect();
        arr
    }
    
    pub fn poly_count_triangles(&self, matrix: &Vec<i128>) -> SparsePolynomial<F251, SparseTerm> {
        let a = matrix.clone();
    
        let len: usize = (matrix.len() as f32).sqrt() as usize;
    
    
        let var_num = (len as f32).log2() as usize;
    
        let x_indexes = Triangles::gen_var_indexes(0, var_num);
        let y_indexes = Triangles::gen_var_indexes(x_indexes.last().unwrap() + 1, var_num);
        println!("x indexes {:?}", x_indexes);
        println!("y indexes {:?}", y_indexes);
        let mut xy_indexes: Vec<usize> = x_indexes.clone();
        xy_indexes.append(&mut y_indexes.clone());
        println!("xy indexes {:?}", xy_indexes);
    
        let z_indexes = Triangles::gen_var_indexes(y_indexes.last().unwrap() + 1, var_num);
        println!("z indexes {:?}", z_indexes);
    
        let mut yz_indexes: Vec<usize> = y_indexes.clone();
        yz_indexes.append(&mut z_indexes.clone());
        println!("yz indexes {:?}", yz_indexes);
        
        let mut xz_indexes: Vec<usize> = x_indexes.clone();
        xz_indexes.append(&mut z_indexes.clone());
        println!("xz indexes {:?}", xz_indexes);
        let poly_exist_xy = poly_slow_mle(&a, &xy_indexes);
        let poly_exist_yz = poly_slow_mle(&a, &yz_indexes);
        let poly_exist_xz = poly_slow_mle(&a, &xz_indexes);
        let poly_exist = naive_mul(&naive_mul(&poly_exist_xy, &poly_exist_yz), &poly_exist_xz);
        poly_exist
    }
}

// bench naive calculation
// bench MLE 

// bench sum check
// prover
// verifier