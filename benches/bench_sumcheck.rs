use std::rc::Rc;
use thaler::small_fields::{F251};
use thaler::sumcheck::{self, SumCheckPolynomial, Prover, UniPoly, Verifier};
use thaler::triangles::{Triangles, MLEAlgorithm, Matrix};
use criterion::{criterion_group, criterion_main, Criterion, BenchmarkId, BenchmarkGroup};

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

fn bench_verifier_steps(c: &mut Criterion) {
	let mut group = c.benchmark_group("verifier");

	let matrix_sizes = [4, 8, 16, 32];
	let eval_types = &[
		MLEAlgorithm::SlowMLE,
		MLEAlgorithm::DynamicMLE,
		MLEAlgorithm::StreamMLE,
	];
	for size in matrix_sizes {
		let matrix = Matrix::new(thaler::utils::gen_matrix(size));

		let g: Triangles<F251> = matrix.derive_mle(MLEAlgorithm::DynamicMLE);
		group.bench_function(
			BenchmarkId::new::<&str, usize>("count triangles with size", g.num_vars()), 
			|b| {
				b.iter(|| {
					g.count()
				});
			}
		);

		for eval_type in eval_types {
			let g: Triangles<F251> = matrix.derive_mle(eval_type.clone());
			let p: Prover<Triangles<F251>> = sumcheck::Prover::new(&g.clone());
			let g_sum = p.slow_sum_g();
			let lookup = build_gi_lookup(&g);
		
			let mut v: Verifier<UniPoly, Triangles<F251>> = Verifier::new(g_sum, Rc::new(g.clone()));
			v.random_func(|| F251::from(1));
			
			group.bench_function(
				BenchmarkId::new::<&str, usize>(&format!("verifier with eval {:?} in vars", eval_type), lookup.len()), 
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
		}
	}

	group.finish()
}

fn benchmarks(c: &mut Criterion) {
	bench_verifier_steps(c);

	// bench_prover_lookup_build(c);
}

fn bench_prover_lookup_build(c: &mut Criterion) {
	let mut group = c.benchmark_group("prover");

	let matrix_sizes = [4];
	for size in matrix_sizes {
		let matrix = Matrix::new(thaler::utils::gen_matrix(size));
		let g: Triangles<F251> = matrix.derive_mle(MLEAlgorithm::DynamicMLE);
		group.bench_function(
			BenchmarkId::new::<&str, usize>("prover for size", size), 
			|b| {
				b.iter(|| {
					build_gi_lookup(&g);
				});
			}
		);
	}
	group.finish();
}

criterion_group!(benches, benchmarks);
criterion_main!(benches);