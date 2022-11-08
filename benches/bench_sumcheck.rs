use ark_ff::{Field, One};
use criterion::measurement::WallTime;
use thaler::lagrange::{MultilinearExtension, UniPoly};
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
use thaler::sumcheck::{self, SumCheckPolynomial, Prover, Verifier};
use thaler::triangles::{TriangleMLE, TriangleGraph};
use criterion::{criterion_group, criterion_main, Criterion, BenchmarkId, BenchmarkGroup};

type TestField = F251;

type SlowMLE = ValueBasedMultilinearExtension<TestField, SlowEvaluationMethod>;
type StreamMLE = ValueBasedMultilinearExtension<TestField, StreamEvaluationMethod>;
type DynamicMLE = ValueBasedMultilinearExtension<TestField, DynamicEvaluationMethod>;

type PolyMLE = PolyMultilinearExtension<TestField>;

// a gi lookup table
fn build_gi_lookup<F, E>(g: &TriangleMLE<F, E>) -> Vec<UniPoly<F>> where F: Field, E: MultilinearExtension<F> {
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

fn bench_verifier_steps<'a, E>(mut group: BenchmarkGroup<'a, WallTime>, matrix_sizes: &Vec<usize>) 
-> BenchmarkGroup<'a, WallTime> where E: MultilinearExtension<TestField> {
	for size in matrix_sizes {
		let matrix = TriangleGraph::new(thaler::utils::gen_matrix(*size));
		let num_vars = matrix.one_dimension_size() * 3;

		let g: TriangleMLE<TestField, E> = matrix.derive_mle();
		let g_sum = g.hypercube_sum();
		let lookup = build_gi_lookup(&g);
	
		let mut v: Verifier<TestField, UniPoly<TestField>, TriangleMLE<TestField, E>> = Verifier::new(g_sum, g.clone());
		v.random_func(|| TestField::one());

		assert_eq!(num_vars, lookup.len());
		
		group.bench_function(
			BenchmarkId::new::<&str, usize>(&format!("verifier with mle algo {:?} in vars", std::any::type_name::<E>()), lookup.len()), 
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

fn bench_prover_lookup_build<'a, E>(mut group: BenchmarkGroup<'a, WallTime>, matrix_sizes: &Vec<usize>) 
-> BenchmarkGroup<'a, WallTime> where E: MultilinearExtension<TestField> {
	for size in matrix_sizes {
		let matrix = TriangleGraph::new(thaler::utils::gen_matrix(*size));
		let g: TriangleMLE<TestField, E> = matrix.derive_mle();
		group.bench_function(
			BenchmarkId::new::<&str, usize>(&format!("prover with mle algo {:?} in matrix size", std::any::type_name::<E>()), *size), 
			|b| {
				b.iter(|| {
					build_gi_lookup(&g);
				});
			}
		);
	}
	group
}

fn bench_triangle_graph<'a>(mut group: BenchmarkGroup<'a, WallTime>, matrix_sizes: &Vec<usize>) 
-> BenchmarkGroup<'a, WallTime> {
	for size in matrix_sizes {
		let matrix = TriangleGraph::<TestField>::new(thaler::utils::gen_matrix(*size));

		let num_vars = matrix.one_dimension_size() * 3;
		group.bench_function(
			BenchmarkId::new::<&str, usize>("triangle graph with vars", num_vars), 
			|b| {
				b.iter(|| {
					matrix.count()
				});
			}
		);
	}
	group
}

fn bench_verifier(c: &mut Criterion) {
	let mut group = c.benchmark_group("verifier");
	let matrix_sizes = vec![4, 8, 16, 32];

	group = bench_triangle_graph(group, &matrix_sizes);

	group = bench_verifier_steps::<SlowMLE>(group, &matrix_sizes);
	group = bench_verifier_steps::<StreamMLE>(group, &matrix_sizes);
	group = bench_verifier_steps::<DynamicMLE>(group, &matrix_sizes);
	group = bench_verifier_steps::<PolyMLE>(group, &matrix_sizes);
	group.finish();
}

fn bench_prover(c: &mut Criterion) {
	let mut group = c.benchmark_group("prover");
	let matrix_sizes = vec![4, 8, 16];

	group = bench_prover_lookup_build::<SlowMLE>(group, &matrix_sizes);
	group = bench_prover_lookup_build::<StreamMLE>(group, &matrix_sizes);
	group = bench_prover_lookup_build::<DynamicMLE>(group, &matrix_sizes);
	group = bench_prover_lookup_build::<PolyMLE>(group, &matrix_sizes);
	group.finish();
}

fn benchmarks(c: &mut Criterion) {
	bench_verifier(c);
	bench_prover(c);
}

criterion_group!(benches, benchmarks);
criterion_main!(benches);