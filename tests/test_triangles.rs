#[macro_use]
extern crate lazy_static;

use rstest::rstest;
use thaler::{small_fields::F251, triangles::{TriangleGraph, MLEAlgorithm, TriangleMLE}, utils::convert_field};
use thaler::sumcheck::SumCheckPolynomial;

fn convert_vec (m: &[i32]) -> Vec<Vec<F251>> {
    let len = (m.len() as f64).sqrt() as usize;

    let mut matrix = Vec::new();
    for (i, e) in m.iter().enumerate() {
        if i % len == 0 {
            let row = Vec::new();
            matrix.push(row);
        }
        let v = F251::from(*e);
        let row = matrix.last_mut().unwrap();
        
        row.push(v);
    }

    matrix
}

lazy_static! {
	static ref M_1: Vec<Vec<F251>> = convert_vec(&[
        0, 1, 1, 0,
        1, 0, 1, 1,
        1, 1, 0, 1,
        0, 1, 1, 0,
    ]);
}

#[rstest]
#[case(&M_1, F251::from(2))]
fn naive_count_test(
	#[case] m: &[Vec<F251>],
	#[case] expected: F251,
) {
    let matrix = TriangleGraph::new(m.to_vec());
	assert_eq!(matrix.count(), expected);
}

#[rstest]
fn var_fixed_evaluate_test() {
    let m = convert_vec(&[
        0, 1, 
        1, 0, 
    ]);
    let matrix = TriangleGraph::new(m.to_vec());
    let triangle: TriangleMLE<F251> = matrix.derive_mle(MLEAlgorithm::SlowMLE);
    let point = convert_field(&[/*x*/0,/*y*/0,/*z*/1]);
    let p = triangle.var_fixed_evaluate(0, point);
    println!("poly {:?}", p);
}
