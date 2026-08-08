//! The version this crate tells people to depend on, checked against the version it is.
//!
//! `AGENTS.md` and `README.md` both print an install line. Those lines are prose, and prose
//! stating a number that nothing checks is how this repository has shipped stale counts more
//! than once — the FRICTION totals, the example count, the test count, and `AGENTS.md` itself
//! saying the Python bindings were "not yet on PyPI" for a day after they were.
//!
//! A release bumps `Cargo.toml` and has no reason to touch a markdown file, so this is the
//! failure mode with the least friction of all of them: nothing goes wrong, the number is simply
//! a release behind, and the first person to notice is someone who copied it.

/// Where the repository root is, relative to this crate.
///
/// `None` when the files are absent, which is the packaged case: `cargo package` builds the
/// crate from a tarball that contains no `AGENTS.md`, and a test that failed there would make
/// publishing impossible for a reason that has nothing to do with the crate. Absent is
/// therefore skipped and *present but wrong* is a failure — the distinction that matters.
fn repo_file(name: &str) -> Option<String> {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(name);
    std::fs::read_to_string(path).ok()
}

/// The `major.minor` a caller should write, which is what the install lines quote.
fn series() -> String {
    let v = env!("CARGO_PKG_VERSION");
    let mut parts = v.split('.');
    let major = parts.next().expect("a version has a major");
    let minor = parts.next().expect("a version has a minor");
    format!("{major}.{minor}")
}

/// **Every `dualis = "x.y"` in the documentation is the version this crate actually is.**
///
/// Matches the dependency line wherever it appears rather than at a fixed location, so moving
/// the snippet does not silently switch the check off — which is the other way a check like this
/// dies.
#[test]
fn the_install_lines_quote_the_current_series() {
    let want = series();
    let needle = "dualis = \"";
    let mut checked = 0;

    for file in ["AGENTS.md", "README.md", "CONTRIBUTING.md", "CLAUDE.md"] {
        let Some(text) = repo_file(file) else {
            continue;
        };
        for (line_no, line) in text.lines().enumerate() {
            let Some(at) = line.find(needle) else {
                continue;
            };
            let rest = &line[at + needle.len()..];
            let quoted = rest.split('"').next().unwrap_or("");
            // Only the bare series form is on trial. A line pinning an exact patch, or a path
            // dependency with a version beside it, is saying something else deliberately.
            if quoted.chars().filter(|c| *c == '.').count() != 1 {
                continue;
            }
            checked += 1;
            assert_eq!(
                quoted,
                want,
                "{file}:{} says dualis = {quoted:?} but this crate is {}. Bump the prose \
                 with the release.",
                line_no + 1,
                env!("CARGO_PKG_VERSION")
            );
        }
    }

    // A check that stopped finding anything would pass forever. If the install lines are
    // reworded out of this shape, this should fail and be rewritten rather than quietly retire.
    if repo_file("AGENTS.md").is_some() {
        assert!(
            checked > 0,
            "no `dualis = \"x.y\"` line found in the documentation — either it was reworded, \
             in which case update this test, or it was deleted, in which case a caller no \
             longer has one to copy"
        );
    }
}
