use std::process::Command;

fn main() {
    // Get the short git hash
    let git_hash = Command::new("git")
        .args(["rev-parse", "--short", "HEAD"])
        .output()
        .ok()
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .map(|s| s.trim().to_string())
        .unwrap_or_else(|| "unknown".to_string());

    // Check for uncommitted changes
    let is_dirty = Command::new("git")
        .args(["status", "--porcelain"])
        .output()
        .map(|o| !o.stdout.is_empty())
        .unwrap_or(false);

    let version_suffix = if is_dirty {
        format!("{}-dev", git_hash)
    } else {
        git_hash
    };

    // Make it available as GIT_VERSION env var at compile time
    println!("cargo:rustc-env=GIT_VERSION={}", version_suffix);

    // Re-run if git HEAD changes or if files change
    println!("cargo:rerun-if-changed=.git/HEAD");
    println!("cargo:rerun-if-changed=.git/index");
}
