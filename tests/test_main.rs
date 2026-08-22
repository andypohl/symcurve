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

#[test]
fn test_max_memory_changes_chunking_without_changing_the_answer() {
    // A small budget makes chunks smaller, which changes how far twist accumulates
    // within a chunk and so perturbs the last digit of some scores. Positions must be
    // unaffected and values must agree to well within the precision bigWig stores.
    let input = write_fasta("budget");
    let big = tmp_path("budget-big.bedGraph");
    let small = tmp_path("budget-small.bedGraph");

    for (out, budget) in [(&big, "8G"), (&small, "1M")] {
        let status = Command::new(EXE)
            .args([
                input.to_str().unwrap(),
                out.to_str().unwrap(),
                "--max-memory",
                budget,
            ])
            .status()
            .expect("failed to run");
        assert!(status.success(), "run with --max-memory {budget} failed");
    }

    let a = std::fs::read_to_string(&big).unwrap();
    let b = std::fs::read_to_string(&small).unwrap();
    let a: Vec<&str> = a.lines().collect();
    let b: Vec<&str> = b.lines().collect();
    assert_eq!(a.len(), b.len(), "different number of scores");
    assert!(!a.is_empty());

    for (la, lb) in a.iter().zip(&b) {
        let fa: Vec<&str> = la.split('\t').collect();
        let fb: Vec<&str> = lb.split('\t').collect();
        assert_eq!(fa[..3], fb[..3], "positions differ");
        let (va, vb) = (fa[3].parse::<f64>().unwrap(), fb[3].parse::<f64>().unwrap());
        let rel = (va - vb).abs() / va.abs().max(1e-12);
        assert!(rel < 1e-5, "{va} vs {vb} differ by {rel}");
    }

    std::fs::remove_file(&input).ok();
    std::fs::remove_file(&big).ok();
    std::fs::remove_file(&small).ok();
}

#[test]
fn test_bad_max_memory_is_rejected() {
    let out = tmp_path("badmem-out.bedGraph");
    let output = Command::new(EXE)
        .args(["in.fa", out.to_str().unwrap(), "--max-memory", "lots"])
        .output()
        .expect("failed to run");
    assert!(!output.status.success());
    let err = String::from_utf8_lossy(&output.stderr);
    assert!(
        err.contains("8G"),
        "error should show the expected form: {err}"
    );
}

#[test]
fn test_symmetry_stage_end_to_end() {
    let input = write_fasta("sym");
    let curv = tmp_path("sym-curv.bedGraph");
    let symm = tmp_path("sym-symm.bedGraph");

    for (out, extra) in [
        (&curv, vec![]),
        (&symm, vec!["--stage", "symmetry", "--symcurve-win", "20"]),
    ] {
        let mut args = vec![input.to_str().unwrap(), out.to_str().unwrap()];
        args.extend(extra);
        let output = Command::new(EXE)
            .args(&args)
            .output()
            .expect("failed to run");
        assert!(output.status.success(), "{:?}", output);
    }

    let curv_lines: Vec<String> = std::fs::read_to_string(&curv)
        .unwrap()
        .lines()
        .map(str::to_string)
        .collect();
    let symm_lines: Vec<String> = std::fs::read_to_string(&symm)
        .unwrap()
        .lines()
        .map(str::to_string)
        .collect();

    // Symmetry costs 20 more on each side of every piece than curvature does. Three
    // pieces of 150, 100 and 100 bases: the 100s lose everything beyond the margin.
    assert!(!symm_lines.is_empty(), "symmetry produced nothing");
    assert!(
        symm_lines.len() < curv_lines.len(),
        "symmetry should yield fewer scores than curvature"
    );

    // The first symmetry score sits 20 further in than the first curvature score.
    let first_pos =
        |lines: &[String]| -> u64 { lines[0].split('\t').nth(1).unwrap().parse().unwrap() };
    assert_eq!(first_pos(&symm_lines), first_pos(&curv_lines) + 20);

    // Scores are non-negative, and at least one dyad scored: symmetry is zero except at
    // strict local minima, so an all-zero file would mean the minimum test never fired.
    let values: Vec<f64> = symm_lines
        .iter()
        .map(|l| l.split('\t').nth(3).unwrap().parse().unwrap())
        .collect();
    assert!(values.iter().all(|v| *v >= 0.0), "negative symmetry score");
    assert!(values.iter().any(|v| *v > 0.0), "no dyad scored at all");

    std::fs::remove_file(&input).ok();
    std::fs::remove_file(&curv).ok();
    std::fs::remove_file(&symm).ok();
}

#[test]
fn test_roll_matrix_is_selectable() {
    let input = write_fasta("roll");
    let simple = tmp_path("roll-simple.bedGraph");
    let active = tmp_path("roll-active.bedGraph");
    for (out, roll) in [(&simple, "simple"), (&active, "active")] {
        let status = Command::new(EXE)
            .args([
                input.to_str().unwrap(),
                out.to_str().unwrap(),
                "--roll",
                roll,
            ])
            .status()
            .expect("failed to run");
        assert!(status.success());
    }
    let a = std::fs::read_to_string(&simple).unwrap();
    let b = std::fs::read_to_string(&active).unwrap();
    assert_ne!(a, b, "the two roll matrices produced identical output");
    std::fs::remove_file(&input).ok();
    std::fs::remove_file(&simple).ok();
    std::fs::remove_file(&active).ok();
}
