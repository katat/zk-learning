#[macro_use]
extern crate lazy_static;

use rstest::rstest;
use thaler::{
    small_fields::F251, 
    triangles::{TriangleGraph, TriangleMLE}, 
    utils::convert_field, 
    mles::{
        ValueBasedMultilinearExtension, 
        value_mle::{methods::DynamicEvaluationMethod}
    }, lagrange::UniPoly, 
};
use thaler::sumcheck::SumCheckPolynomial;
use ark_poly::polynomial::Polynomial;

type TestField = F251;

fn convert_vec (m: &[i32]) -> Vec<Vec<TestField>> {
    let len = (m.len() as f64).sqrt() as usize;

    let mut matrix = Vec::new();
    for (i, e) in m.iter().enumerate() {
        if i % len == 0 {
            let row = Vec::new();
            matrix.push(row);
        }
        let v = TestField::from(*e);
        let row = matrix.last_mut().unwrap();
        
        row.push(v);
    }

    matrix
}

lazy_static! {
	static ref M_1: Vec<Vec<TestField>> = convert_vec(&[
        0, 1, 1, 0,
        1, 0, 1, 1,
        1, 1, 0, 1,
        0, 1, 1, 0,
    ]);
}

#[rstest]
#[case(&M_1, TestField::from(2))]
fn naive_count_test(
	#[case] m: &[Vec<TestField>],
	#[case] expected: TestField,
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
    let triangle: TriangleMLE<TestField, ValueBasedMultilinearExtension<TestField, DynamicEvaluationMethod>> = matrix.derive_mle();
    let point = convert_field(&[/*x*/0,/*y*/0,/*z*/1]);
    let p: UniPoly<TestField> = triangle.var_fixed_evaluate(0, point);
    assert_eq!(p.evaluate(&TestField::from(1)), TestField::from(0));
}
