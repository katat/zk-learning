use std::rc::Rc;
use ark_ff::Field;
use criterion::measurement::WallTime;
use thaler::lagrange::{MultilinearExtension};
use thaler::mles::{
	value_mle::{
		ValueBasedMultilinearExtension,
		methods::{
			DynamicEvaluationMethod,
			SlowEvaluationMethod,
			StreamEvaluationMethod,
		}
	},
	PolyMultilinearExtension
};
use thaler::small_fields::{F251};
use thaler::sumcheck::{self, SumCheckPolynomial, Prover, UniPoly, Verifier};
use thaler::triangles::{TriangleMLE, TriangleGraph};
use criterion::{criterion_group, criterion_main, Criterion, BenchmarkId, BenchmarkGroup};

type TestField = F251;

// a gi lookup table
fn build_gi_lookup<F, E>(g: &TriangleMLE<F, E>) -> Vec<sumcheck::UniPoly<F>> where F: Field, E: MultilinearExtension<F> {
	let r: Option<F> = Some(1u32.into());
	let mut lookup= vec![];
	let mut p: Prover<F, TriangleMLE<F, E>> = sumcheck::Prover::<F, TriangleMLE<F, E>>::new(g);

	let mut gi = p.gen_uni_polynomial(None);
	lookup.push(gi.clone());
	for _ in 1..p.g.num_vars() {
		gi = p.gen_uni_polynomial(r);
		lookup.push(gi.clone());
	}
	lookup
}

fn bench_verifier_steps<'a, F, E>(mut group: BenchmarkGroup<'a, WallTime>, eval_type: &'a str) -> BenchmarkGroup<'a, WallTime> where F: Field, E: MultilinearExtension<F> {

	let matrix_sizes = [4, 8, 16, 32];

	for size in matrix_sizes {
		let matrix = TriangleGraph::new(thaler::utils::gen_matrix(size));
		let num_vars = matrix.one_dimension_size() * 3;

		let g: TriangleMLE<F, E> = matrix.derive_mle();
		let p: Prover<F, TriangleMLE<F, E>> = sumcheck::Prover::new(&g.clone());
		let g_sum = p.slow_sum_g();
		let lookup = build_gi_lookup(&g);
	
		let mut v: Verifier<F, UniPoly<F>, TriangleMLE<F, E>> = Verifier::new(g_sum, Rc::new(g.clone()));
		v.random_func(|| F::one());

		assert_eq!(num_vars, lookup.len());
		
		group.bench_function(
			BenchmarkId::new::<&str, usize>(&format!("verifier with mle algo {:?} in vars", eval_type), lookup.len()), 
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

	group

}

fn benchmarks(c: &mut Criterion) {
	let mut group = c.benchmark_group("verifier");
	let matrix_sizes = [4, 8, 16, 32];

	for size in matrix_sizes {
		let matrix = TriangleGraph::<TestField>::new(thaler::utils::gen_matrix(size));

		let num_vars = matrix.one_dimension_size() * 3;
		group.bench_function(
			BenchmarkId::new::<&str, usize>("count triangles in vars", num_vars), 
			|b| {
				b.iter(|| {
					matrix.count()
				});
			}
		);
	}
	group = bench_verifier_steps::<TestField, ValueBasedMultilinearExtension<TestField, SlowEvaluationMethod>>(group, "slow");
	group = bench_verifier_steps::<TestField, ValueBasedMultilinearExtension<TestField, StreamEvaluationMethod>>(group, "stream");
	group = bench_verifier_steps::<TestField, ValueBasedMultilinearExtension<TestField, DynamicEvaluationMethod>>(group, "dynamic");
	group = bench_verifier_steps::<TestField, PolyMultilinearExtension<TestField>>(group, "polynomial");
	group.finish();

	let mut group = c.benchmark_group("prover");
	group = bench_prover_lookup_build::<TestField, ValueBasedMultilinearExtension<TestField, SlowEvaluationMethod>>(group, "slow");
	group = bench_prover_lookup_build::<TestField, ValueBasedMultilinearExtension<TestField, StreamEvaluationMethod>>(group, "stream");
	group = bench_prover_lookup_build::<TestField, ValueBasedMultilinearExtension<TestField, DynamicEvaluationMethod>>(group, "dynamic");
	group = bench_prover_lookup_build::<TestField, PolyMultilinearExtension<TestField>>(group, "polynomial");
	group.finish();
}

fn bench_prover_lookup_build<'a, F, E>(mut group: BenchmarkGroup<'a, WallTime>, eval_type: &'a str) -> BenchmarkGroup<'a, WallTime> where F: Field, E: MultilinearExtension<F> {
	let matrix_sizes = [4, 8, 16];
	for size in matrix_sizes {
		let matrix = TriangleGraph::new(thaler::utils::gen_matrix(size));
		let g: TriangleMLE<F, E> = matrix.derive_mle();
		group.bench_function(
			BenchmarkId::new::<&str, usize>(&format!("prover with mle algo {:?} in matrix size", eval_type), size), 
			|b| {
				b.iter(|| {
					build_gi_lookup(&g);
				});
			}
		);
	}
	group
}

criterion_group!(benches, benchmarks);
criterion_main!(benches);