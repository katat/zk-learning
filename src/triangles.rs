use std::{iter, ops::Div};

use ark_ff::{Zero, PrimeField};
use ark_poly::{multivariate::{SparsePolynomial, SparseTerm}, Polynomial};

use crate::{small_fields::{F251}, lagrange::{poly_slow_mle, naive_mul}};

#[derive(Debug, Clone)]
pub struct Matrix {
    vec: Vec<Vec<F251>>
}

#[derive(Debug, Clone)]
pub struct Triangles {
	pub matrix: Matrix,
}

fn println_matrix(matrix: Vec<Vec<F251>>) {
    (0..matrix.len()).for_each(|i| {
        for _j in 0..matrix[i].len() {
            // print!("{} ", matrix[i][j].into_bigint().as_ref()[0]);
        }
        //println!("");
    });
}

impl Matrix {
    pub fn new(m: Vec<Vec<F251>>) -> Self {
        Matrix {
            vec: m
        }
    }
    pub fn flatten(&self) -> Vec<F251>{
        self.vec.iter().flatten().cloned().collect()
    }

    pub fn size(&self) -> usize{
        (self.flatten().len() as f32).sqrt() as usize
    }

    pub fn var_num(&self) -> usize {
        (self.size() as f32).log2() as usize
    }

    pub fn get(&self, x: usize, y: usize) -> F251 {
        self.vec[x][y]
    }
}

// todo split into triangles_mle and triangles_matrix traits?
impl Triangles {
    //todo make Vec<Vec<F251>> type generic
    pub fn new(matrix: Vec<Vec<F251>>) -> Self {
        Triangles {
            matrix: Matrix { vec: matrix }
        }
    }

    pub fn multiply(&self, matrix: Matrix) -> Matrix {
        let size = self.matrix.size();
        let row: Vec<F251> = iter::repeat(F251::zero()).take(size).collect();
        let mut result_matrix: Vec<Vec<F251>> = iter::repeat(row).take(size).collect();
        for i in 0..size {
            for j in 0..size {
                let mut elm = F251::zero();
                for k in 0..size {
                    elm += self.matrix.get(i, k) * matrix.get(k, j);
                    // println!("{}", matrix[k][j]);
                }
                result_matrix[i][j] = elm;
            }
        }

        Matrix {
            vec: result_matrix
        }
    }

    pub fn count(&self) -> F251 {
        let a2 = self.multiply(self.matrix.clone());
        let a3 = self.multiply(a2);
        let mut number = F251::zero();

        for i in 0..a3.size() {
            number += a3.get(i, i);
        }

        number.div(F251::from(6))
    }

    pub fn gen_var_indexes (start_index: usize, var_num: usize) -> Vec<usize> {
        let arr: Vec<usize> = (0..var_num).map(|x| x + start_index).collect();
        arr
    }

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
            .map(|i| i.to_string().parse::<i128>().unwrap())
            .collect()
    }
    
    pub fn poly_count_triangles(&self) -> SparsePolynomial<F251, SparseTerm> {
        let a: Vec<F251> = self.matrix.flatten();
    
        let var_num = self.matrix.var_num();
    
        // encapsulate them as first/second/third...
        let x_start_index = 0;
        let y_start_index = var_num;
        let z_start_index = var_num * 2;
        let x_indexes = Triangles::gen_var_indexes(x_start_index, var_num);
        let y_indexes = Triangles::gen_var_indexes(y_start_index, var_num);
        let mut xy_indexes: Vec<usize> = x_indexes.clone();
        xy_indexes.append(&mut y_indexes.clone());
    
        let mut z_indexes = Triangles::gen_var_indexes(z_start_index, var_num);
    
        let mut yz_indexes: Vec<usize> = y_indexes;
        yz_indexes.append(&mut z_indexes.clone());
        
        let mut xz_indexes: Vec<usize> = x_indexes;
        xz_indexes.append(&mut z_indexes);

        //clean up
        
        let converted_a = a.into_iter().map(|e| e.into_bigint().as_ref()[0] as i128).collect::<Vec<i128>>();
        let poly_exist_xy = poly_slow_mle(&converted_a, &xy_indexes);
        let poly_exist_yz = poly_slow_mle(&converted_a, &yz_indexes);
        let poly_exist_xz = poly_slow_mle(&converted_a, &xz_indexes);

        println!("poly xz {}", poly_exist_xz.terms.len());
        naive_mul(&naive_mul(&poly_exist_xy, &poly_exist_yz), &poly_exist_xz)
    }

    pub fn count_by_mle(&self) -> i128 {
        let poly_exist = self.poly_count_triangles();
        let len = self.matrix.size();
        let var_num = self.matrix.var_num();
        let mut total_triangles = 0;
    
        for x in 0..len {
            for y in 0..len {
                for z in 0..len {
                    let xyz_bin = Triangles::convert_bin_vec(Triangles::convert_bin_z(x, y, z, var_num));
                    let r: Vec<F251> = xyz_bin.iter().map(|i| F251::from(*i)).collect();
    
                    let result = poly_exist.evaluate(&r);
                    let exist = result.into_bigint().as_ref()[0];
    
                    if exist != 0 {
                        //println!("exist {} at x: {}, y: {}, z: {}", exist, x, y, z);
                    }
                    total_triangles += exist;
                }
            }
        }
    
        total_triangles as i128 / 6
    
    }
}

// bench naive calculation
// bench MLE 

// bench sum check
// prover
// verifier