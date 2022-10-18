#![feature(test)]

#[macro_use]
extern crate lazy_static;

extern crate test;
use ark_ff::{Zero, One};
use ark_poly::polynomial::{Polynomial};
use test::Bencher;
use thaler::small_fields::{F251};
use thaler::sumcheck::{self, SumCheckPolynomial, Prover, UniPoly, Verifier};
use thaler::triangles::Triangles;

lazy_static! {
	static ref M_1: Vec<Vec<F251>> = thaler::utils::gen_matrix(16);
	static ref G_1: Triangles = thaler::triangles::Triangles::new(M_1.clone());
	static ref P_1: Prover<Triangles> = sumcheck::Prover::new(&G_1);
	static ref G_1_SUM: F251 = P_1.slow_sum_g();
}

// a gi lookup table
fn build_gi_lookup() -> Vec<sumcheck::UniPoly> {
	let r: Option<F251> = Some(2u32.into());
	let mut lookup = vec![];
	let mut p: Prover<Triangles> = sumcheck::Prover::<Triangles>::new(&G_1);
	// OVERHEAD
	let mut gi = p.gen_uni_polynomial(None);
	lookup.push(gi.clone());
	for _ in 1..p.g.num_vars() {
		gi = p.gen_uni_polynomial(r);
		lookup.push(gi.clone());
	}
	lookup
}

// Steps being benchmarked
fn verifier_steps_only(p: &sumcheck::Prover<Triangles>, gi_lookup: &Vec<sumcheck::UniPoly>, r: Option<F251>) {
	let mut v: Verifier<UniPoly, Triangles> = Verifier::new(*G_1_SUM, G_1.clone());
	while v.current_round != Some(sumcheck::Round::Final()) {
		let gi = gi_lookup[v.r_vec.len()].clone();
		v.verify(Some(gi));
	}

	// final round
	v.verify(None);
}

// Verifier benchmark
#[bench]
fn sumcheck_test(b: &mut Bencher) {
	println!("g1 num vars {}", G_1.num_vars());
	let gi_lookup = build_gi_lookup();
	println!("built lookup");
	let r: Option<F251> = Some(2u32.into());
	let p = sumcheck::Prover::<Triangles>::new(&G_1);
	
	b.iter(|| verifier_steps_only(&p, &gi_lookup, r));

	// Prover::init

	// random = SumCheck::Verifier.generate_random()
	// gj = SumCheck::Prover.generate_unipoly(random)
	// enum round = SumCheck::Verifier.verify(gj)
	
	// make it easy to plug in implementations for the following as they are source of overheads
		// SumCheck::Prover.generate_unipoly
		// last round of verifier evaluation
			// make it easy to bench different mle algos
}

// #[bench]
// fn prover_lookup_build(b: &mut Bencher) {
// 	b.iter(|| {
// 		build_gi_lookup();
// 	});
// }

#[bench]
fn sumcheck_last_round_test(b: &mut Bencher) {
	let r: Option<F251> = Some(2u32.into());
	
	b.iter(|| {
		G_1.evaluate(&vec![r.unwrap(); G_1.num_vars()]);
	});
}

// #[bench]
// fn slow_sumcheck_test(b: &mut Bencher) {
// 	let p = sumcheck::Prover::<Triangles>::new(&G_1);
// 	b.iter(|| p.slow_sum_g());
// }

#[bench]
fn naive_count(b: &mut Bencher) {
	b.iter(|| {
		G_1.count();
	});
}