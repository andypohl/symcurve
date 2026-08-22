//! Integration tests driving the built binary.
use std::process::Command;

/// Cargo sets this to the path of the binary under test, which is more robust than
/// assuming a debug profile and a fixed target directory.
const EXE: &str = env!("CARGO_BIN_EXE_symcurve");

fn tmp_path(name: &str) -> std::path::PathBuf {
    let mut p = std::env::temp_dir();
    p.push(format!("symcurve-it-{}-{}", std::process::id(), name));
    p
}

/// A record with a soft-masked stretch, an N gap and an ambiguity code, so the run
/// exercises gap splitting rather than a single clean piece.
///
/// `tag` must be unique per test: these run in parallel threads of one process, so a
/// shared path would let one test delete the input another is still reading.
fn write_fasta(tag: &str) -> std::path::PathBuf {
    let unit = "CCAACATTTTGACTTTTTGGGAGGGCACTAGCACCTATCTACCCTGAATC";
    let piece_a = unit.repeat(3);
    let piece_b = unit.repeat(2).to_lowercase();
    let piece_c = unit.repeat(2);
    let path = tmp_path(&format!("{tag}-in.fa"));
    std::fs::write(&path, format!(">chrIT\n{piece_a}NNNN{piece_b}R{piece_c}\n")).unwrap();
    path
}

#[test]
fn test_app_runs() {
    let output = Command::new(EXE)
        .arg("-V")
        .output()
        .expect("Failed to execute command");
    assert!(String::from_utf8_lossy(&output.stdout).starts_with("symcurve"));
}

#[test]
fn test_end_to_end_bedgraph() {
    let input = write_fasta("bedgraph");
    let out = tmp_path("bedgraph-out.bedGraph");
    let output = Command::new(EXE)
        .args([input.to_str().unwrap(), out.to_str().unwrap(), "--verbose"])
        .output()
        .expect("failed to run");
    assert!(output.status.success(), "{:?}", output);

    let text = std::fs::read_to_string(&out).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    // Three pieces of 150, 100 and 100 bases, each losing 21 at both ends.
    assert_eq!(lines.len(), (150 - 42) + (100 - 42) + (100 - 42));

    let fields: Vec<&str> = lines[0].split('\t').collect();
    assert_eq!(fields.len(), 4);
    assert_eq!(fields[0], "chrIT");
    assert_eq!(fields[1], "21"); // zero-based start of the first scored base
    assert_eq!(fields[2], "22");
    assert!(fields[3].parse::<f64>().unwrap() > 0.0);

    // Positions must be strictly increasing across the whole file.
    let starts: Vec<u64> = lines
        .iter()
        .map(|l| l.split('\t').nth(1).unwrap().parse().unwrap())
        .collect();
    assert!(
        starts.windows(2).all(|w| w[0] < w[1]),
        "positions not sorted"
    );

    std::fs::remove_file(&input).ok();
    std::fs::remove_file(&out).ok();
}

#[test]
fn test_end_to_end_bigwig() {
    let input = write_fasta("bigwig");
    let out = tmp_path("bigwig-out.bw");
    let status = Command::new(EXE)
        .args([input.to_str().unwrap(), out.to_str().unwrap()])
        .status()
        .expect("failed to run");
    assert!(status.success());

    // bigWig files begin with the magic 0x888FFC26, little-endian on this platform.
    let bytes = std::fs::read(&out).unwrap();
    assert!(bytes.len() > 64);
    assert_eq!(&bytes[..4], &[0x26, 0xFC, 0x8F, 0x88]);

    std::fs::remove_file(&input).ok();
    std::fs::remove_file(&out).ok();
}

#[test]
fn test_unknown_output_extension_is_rejected() {
    let input = write_fasta("reject");
    let out = tmp_path("reject-out.txt");
    let output = Command::new(EXE)
        .args([input.to_str().unwrap(), out.to_str().unwrap()])
        .output()
        .expect("failed to run");
    assert!(!output.status.success());
    let err = String::from_utf8_lossy(&output.stderr);
    assert!(err.contains("cannot tell the output format"), "{err}");
    assert!(
        !out.exists(),
        "a rejected format must not leave a file behind"
    );
    std::fs::remove_file(&input).ok();
}
