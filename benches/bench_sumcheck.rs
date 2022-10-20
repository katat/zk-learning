use std::rc::Rc;
use thaler::small_fields::{F251};
use thaler::sumcheck::{self, SumCheckPolynomial, Prover, UniPoly, Verifier};
use thaler::triangles::Triangles;
use criterion::{criterion_group, criterion_main, Criterion, BenchmarkId};

// a gi lookup table
fn build_gi_lookup(g: &Triangles<F251>) -> Vec<sumcheck::UniPoly> {
	let r: Option<F251> = Some(1u32.into());
	let mut lookup= vec![];
	let mut p: Prover<Triangles<F251>> = sumcheck::Prover::<Triangles<F251>>::new(g);
	// OVERHEAD
	let mut gi = p.gen_uni_polynomial(None);
	lookup.push(gi.clone());
	for _ in 1..p.g.num_vars() {
		gi = p.gen_uni_polynomial(r);
		lookup.push(gi.clone());
	}
	lookup
}

fn bench_verifier(c: &mut Criterion) {
	let mut group = c.benchmark_group("verifier");

	let matrix_sizes = [4, 8, 16, 32];
	for size in matrix_sizes {
		let matrix: Vec<Vec<F251>> = thaler::utils::gen_matrix(size);
		let g: Triangles<F251> = thaler::triangles::Triangles::new(matrix.clone());
		let p: Prover<Triangles<F251>> = sumcheck::Prover::new(&g.clone());
		let g_sum = p.slow_sum_g();
		let lookup = build_gi_lookup(&g);
	
		let mut v: Verifier<UniPoly, Triangles<F251>> = Verifier::new(g_sum, Rc::new(g.clone()));
		v.random_func(|| F251::from(1));
		
		//todo figure an easy way to add bench for different algo implementations
		group.bench_function(
			BenchmarkId::new::<&str, usize>("verifier with vars", lookup.len()), 
			|b| {
				b.iter(|| {
					let mut v = v.clone();
					(0..lookup.len()).for_each(|i| {
						let gi = lookup[i].clone();
						v.verify(Some(gi))
					});
					v.verify(None);
				});
			}
		);
		group.bench_function(
			BenchmarkId::new::<&str, usize>("count triangles with size", lookup.len()), 
			|b| {
				b.iter(|| {
					g.count()
				});
			}
		);
	}


	group.finish();
}

// #[bench]
// fn prover_lookup_build(b: &mut Bencher) {
// 	b.iter(|| {
// 		build_gi_lookup();
// 	});
// }

// #[bench]
// fn sumcheck_last_round_test(b: &mut Bencher) {
// 	let r: Option<F251> = Some(2u32.into());
	
// 	b.iter(|| {
// 		G_1.evaluate(&vec![r.unwrap(); G_1.num_vars()]);
// 	});
// }

// #[bench]
// fn slow_sumcheck_test(b: &mut Bencher) {
// 	let p = sumcheck::Prover::<Triangles>::new(&G_1);
// 	b.iter(|| p.slow_sum_g());
// }

// #[bench]
// fn naive_count(b: &mut Bencher) {
// 	b.iter(|| {
// 		G_1.count();
// 	});
// }

criterion_group!(benches, bench_verifier);
criterion_main!(benches);