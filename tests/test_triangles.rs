#[macro_use]
extern crate lazy_static;

use rstest::rstest;
use thaler::{small_fields::F251, triangles::Triangles};

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
	#[case] m: &Vec<Vec<F251>>,
	#[case] expected: F251,
) {
    let triangles = Triangles::new(m.to_vec());
	assert_eq!(triangles.count(), F251::from(expected));
}

#[rstest]
#[case(&M_1, F251::from(2))]
fn slow_mle_count_test(
	#[case] m: &Vec<Vec<F251>>,
	#[case] expected: F251,
) {
    let triangles = Triangles::new(m.to_vec());
	assert_eq!(triangles.count_by_mle(), expected);
}

