//! CLI binary for SHA256 STWO prover compatible with csp-benchmarks framework.
//!
//! Subcommands:
//! - `prepare`: Initialize prover state for a given input size
//! - `prove`: Generate a proof
//! - `verify`: Verify a proof
//! - `measure`: Report proof and preprocessing sizes

#![feature(portable_simd)]

use std::{fs, path::PathBuf};

use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use serde::{Deserialize, Serialize};
use sha256::{components::ClaimedSum, prove_sha256, verify_sha256};
use stwo::{core::pcs::PcsConfig, prover::backend::Column};
use tracing::info;
use tracing_subscriber::{fmt, prelude::*, EnvFilter};

/// SHA256 STWO Prover CLI for csp-benchmarks
#[derive(Parser)]
#[command(name = "sha256_prover")]
#[command(about = "SHA256 STWO prover for csp-benchmarks framework")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Prepare prover state for a given input size
    Prepare {
        /// Input size (number of SHA256 hash inputs)
        #[arg(long, env = "INPUT_SIZE")]
        input_size: u64,

        /// Path to write state JSON
        #[arg(long, env = "STATE_JSON")]
        state_json: PathBuf,
    },

    /// Generate a proof
    Prove {
        /// Path to state JSON
        #[arg(long, env = "STATE_JSON")]
        state_json: PathBuf,
    },

    /// Verify a proof
    Verify {
        /// Path to state JSON
        #[arg(long, env = "STATE_JSON")]
        state_json: PathBuf,
    },

    /// Report proof and preprocessing sizes
    Measure {
        /// Path to state JSON
        #[arg(long, env = "STATE_JSON")]
        state_json: PathBuf,

        /// Path to write sizes JSON
        #[arg(long, env = "SIZES_JSON")]
        sizes_json: PathBuf,
    },
}

/// State persisted between benchmark phases
#[derive(Serialize, Deserialize)]
struct BenchmarkState {
    /// Log2 of the number of SHA256 instances
    log_size: u32,
    /// Original input size
    input_size: u64,
    /// Path to the serialized proof
    proof_path: Option<PathBuf>,
    /// Claimed sums from the prover (needed for verification)
    claimed_sum: Option<ClaimedSum>,
}

/// Sizes reported by the measure command
#[derive(Serialize, Deserialize)]
struct BenchmarkSizes {
    /// Size of the proof in bytes
    proof_size: u64,
    /// Size of preprocessing data in bytes
    preprocessing_size: u64,
}

fn main() -> Result<()> {
    // Initialize tracing
    tracing_subscriber::registry()
        .with(fmt::layer())
        .with(EnvFilter::from_default_env())
        .init();

    let cli = Cli::parse();

    match cli.command {
        Commands::Prepare {
            input_size,
            state_json,
        } => cmd_prepare(input_size, state_json),

        Commands::Prove { state_json } => cmd_prove(state_json),

        Commands::Verify { state_json } => cmd_verify(state_json),

        Commands::Measure {
            state_json,
            sizes_json,
        } => cmd_measure(state_json, sizes_json),
    }
}

/// Prepare prover state for a given input size.
///
/// Converts INPUT_SIZE to log_size: log_size = ceil(log2(INPUT_SIZE))
fn cmd_prepare(input_size: u64, state_json: PathBuf) -> Result<()> {
    info!("Preparing state for input_size={}", input_size);

    // Convert input_size to log_size
    // log_size = ceil(log2(input_size))
    let log_size = if input_size == 0 {
        0
    } else {
        (input_size as f64).log2().ceil() as u32
    };

    // Minimum log_size is 4 (for SIMD lanes)
    let log_size = log_size.max(4);

    info!("Computed log_size={}", log_size);

    let state = BenchmarkState {
        log_size,
        input_size,
        proof_path: None,
        claimed_sum: None,
    };

    let state_json_content =
        serde_json::to_string_pretty(&state).context("Failed to serialize state")?;
    fs::write(&state_json, state_json_content).context("Failed to write state JSON")?;

    info!("State written to {:?}", state_json);
    Ok(())
}

/// Generate a proof.
fn cmd_prove(state_json: PathBuf) -> Result<()> {
    info!("Reading state from {:?}", state_json);

    let state_content = fs::read_to_string(&state_json).context("Failed to read state JSON")?;
    let mut state: BenchmarkState =
        serde_json::from_str(&state_content).context("Failed to parse state JSON")?;

    info!("Proving with log_size={}", state.log_size);

    // Generate proof
    sha256::print_enabled_features();
    let (proof, claimed_sum) = prove_sha256(state.log_size, PcsConfig::default());

    // Serialize proof
    let proof_bytes = bincode::serialize(&proof).context("Failed to serialize proof")?;

    // Write proof to file
    let proof_path = state_json.with_extension("proof.bin");
    fs::write(&proof_path, &proof_bytes).context("Failed to write proof")?;

    info!(
        "Proof written to {:?} ({} bytes)",
        proof_path,
        proof_bytes.len()
    );

    // Update state
    state.proof_path = Some(proof_path);
    state.claimed_sum = Some(claimed_sum);

    let state_json_content =
        serde_json::to_string_pretty(&state).context("Failed to serialize state")?;
    fs::write(&state_json, state_json_content).context("Failed to write state JSON")?;

    info!("State updated");
    Ok(())
}

/// Verify a proof.
fn cmd_verify(state_json: PathBuf) -> Result<()> {
    info!("Reading state from {:?}", state_json);

    let state_content = fs::read_to_string(&state_json).context("Failed to read state JSON")?;
    let state: BenchmarkState =
        serde_json::from_str(&state_content).context("Failed to parse state JSON")?;

    let proof_path = state
        .proof_path
        .context("No proof path in state - run prove first")?;
    let claimed_sum = state
        .claimed_sum
        .context("No claimed_sum in state - run prove first")?;

    info!("Verifying proof from {:?}", proof_path);

    // Read and deserialize proof
    let proof_bytes = fs::read(&proof_path).context("Failed to read proof")?;
    let proof = bincode::deserialize(&proof_bytes).context("Failed to deserialize proof")?;

    // Verify
    verify_sha256(proof, state.log_size, &claimed_sum).context("Proof verification failed")?;

    info!("Proof verified successfully!");
    Ok(())
}

/// Report proof and preprocessing sizes.
fn cmd_measure(state_json: PathBuf, sizes_json: PathBuf) -> Result<()> {
    info!("Reading state from {:?}", state_json);

    let state_content = fs::read_to_string(&state_json).context("Failed to read state JSON")?;
    let state: BenchmarkState =
        serde_json::from_str(&state_content).context("Failed to parse state JSON")?;

    let proof_path = state
        .proof_path
        .context("No proof path in state - run prove first")?;

    // Get proof size
    let proof_metadata = fs::metadata(&proof_path).context("Failed to get proof metadata")?;
    let proof_size = proof_metadata.len();

    // Compute preprocessing size
    // The preprocessing consists of the preprocessed trace columns
    // For now, estimate based on log_size
    let preprocessing_size = estimate_preprocessing_size(state.log_size);

    info!("Proof size: {} bytes", proof_size);
    info!("Preprocessing size: {} bytes", preprocessing_size);

    let sizes = BenchmarkSizes {
        proof_size,
        preprocessing_size,
    };

    let sizes_json_content =
        serde_json::to_string_pretty(&sizes).context("Failed to serialize sizes")?;
    fs::write(&sizes_json, sizes_json_content).context("Failed to write sizes JSON")?;

    info!("Sizes written to {:?}", sizes_json);
    Ok(())
}

/// Estimate the preprocessing size based on log_size.
///
/// The preprocessing consists of static lookup tables for SHA256 operations.
/// This is an approximation based on the number and size of preprocessed columns.
fn estimate_preprocessing_size(log_size: u32) -> u64 {
    use sha256::preprocessed::PreProcessedTrace;

    // Generate the preprocessed trace to count columns and sizes
    let preprocessed = PreProcessedTrace::new(log_size);

    // Each column is a CircleEvaluation with domain size 2^log_size
    // Each element is a BaseField (4 bytes = 32 bits for M31)
    let total_elements: u64 = preprocessed
        .trace
        .iter()
        .map(|eval| eval.values.len() as u64)
        .sum();

    // Each BaseField is 4 bytes
    total_elements * 4
}
