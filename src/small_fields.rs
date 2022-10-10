use ark_ff::{Fp64, MontBackend, MontConfig};

#[derive(MontConfig)]
#[modulus = "251"]
#[generator = "2"]
pub struct FrConfig;

pub type F251 = Fp64<MontBackend<FrConfig, 1>>;
