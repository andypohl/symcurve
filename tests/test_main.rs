//! Integration test on main() function.
use std::process::Command;

#[test]
fn test_app_runs() {
    let output = Command::new("target/debug/symcurve")
        .arg("-V")
        .output()
        .expect("Failed to execute command");
    assert!(String::from_utf8_lossy(&output.stdout).starts_with("symcurve"));
}
