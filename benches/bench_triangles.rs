#![feature(test)]

#[macro_use]
extern crate lazy_static;

extern crate test;
use ark_ff::{Zero, One, PrimeField};
use ark_poly::{MVPolynomial, Polynomial};
use test::Bencher;
use thaler::small_fields::{F251, to_u64};
use thaler::{sumcheck, lagrange};
use thaler::triangles::Triangles;

fn convert_vec (m: &[i32]) -> Vec<Vec<F251>> {
    let len = (m.len() as f64).sqrt() as usize;

    let mut matrix = Vec::new();
    for i in 0..m.len() {
        if i % len == 0 {
            let row = Vec::new();
            matrix.push(row);
        }
        let v = F251::from(m[i]);
        let row = matrix.last_mut().unwrap();
        
        row.push(v);
    }

    // //println!("a3 {:#?}", matrix);

    matrix
}

fn gen_matrix (len: usize) -> Vec<Vec<F251>> {
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

lazy_static! {
	// static ref M_1: Vec<Vec<F251>> = convert_vec(&[
    //     0, 1, 1, 0,
    //     1, 0, 1, 1,
    //     1, 1, 0, 1,
    //     0, 1, 1, 0,
    // ]);
	static ref M_1: Vec<Vec<F251>> = gen_matrix(4);
}

#[bench]
fn naive_count(b: &mut Bencher) {
    let triangles = Triangles::new(M_1.to_vec());
	
	b.iter(|| {
		triangles.count();
	});
}

#[bench]
fn count_in_mle(b: &mut Bencher) {
    let triangles = Triangles::new(M_1.to_vec());
	let mle = triangles.poly_count_triangles();
	//println!("mle var {}", mle.degree());
	let points = triangles.gen_points();
	
	b.iter(|| {
		triangles.count_by_mle(&mle, &points);
	});
}

fn //println_matrix(matrix: Vec<Vec<F251>>) {
    for i in 0..matrix.len() {
        for j in 0..matrix[i].len() {
            //print!("{} ", matrix[i][j].into_repr());
        }
        //println!("");
    }
}

#[bench]
fn sumcheck(b: &mut Bencher) {
    let triangles = Triangles::new(M_1.to_vec());
	let mle = triangles.poly_count_triangles();
	let points = triangles.gen_points();
	let expected = triangles.count_by_mle(&mle, &points);
	let expected_1 = triangles.count();
	
	b.iter(|| {
		sumcheck::verify(&mle, expected * F251::from(6));
	});
}

#[bench]
fn poly_count_triangles(b: &mut Bencher) {
    let triangles = Triangles::new(M_1.to_vec());
	
	b.iter(|| {
		triangles.poly_count_triangles();
	});
}
