use ark_ff::{Fp64, MontBackend, MontConfig};

#[derive(MontConfig)]
#[modulus = "251"]
#[generator = "2"]
pub struct FrConfigF251;

#[derive(MontConfig)]
#[modulus = "5"]
#[generator = "2"]
pub struct FrConfigF5;

pub type F5 = Fp64<MontBackend<FrConfigF5, 1>>;
pub type F251 = Fp64<MontBackend<FrConfigF251, 1>>;
