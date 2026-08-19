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

/// **Every `pantometry = "x.y"` in the documentation is the version this crate actually is.**
///
/// Matches the dependency line wherever it appears rather than at a fixed location, so moving
/// the snippet does not silently switch the check off — which is the other way a check like this
/// dies.
#[test]
fn the_install_lines_quote_the_current_series() {
    let want = series();
    let needle = "pantometry = \"";
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
                "{file}:{} says pantometry = {quoted:?} but this crate is {}. Bump the prose \
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
            "no `pantometry = \"x.y\"` line found in the documentation — either it was reworded, \
             in which case update this test, or it was deleted, in which case a caller no \
             longer has one to copy"
        );
    }
}

/// **The subagents' own statements of the current version are current.**
///
/// `invariant-guard` has a section on the one invariant that cannot be fixed after the fact — a
/// published version is permanent — and it opens by stating what is on crates.io and what the
/// tree is. That line was a release behind at 0.3.0 and a release behind again at 0.4.0: the
/// document whose subject is checking things was the last place in the repository stating a
/// version with nothing standing behind it.
///
/// Only the *tree's* version is checked. What is on crates.io is a fact about the outside world
/// that no test here can know, and the agent is told to run `curl` for it rather than trust the
/// prose — which is the right division, and the reason this checks one number and not two.
#[test]
fn the_agents_know_what_version_the_tree_is() {
    let full = env!("CARGO_PKG_VERSION");
    let needle = "the tree is ";
    let mut checked = 0;

    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../.claude/agents");
    let Ok(entries) = std::fs::read_dir(&dir) else {
        return; // packaged build: the agents are not part of any crate
    };
    let mut files: Vec<_> = entries.filter_map(|e| e.ok().map(|e| e.path())).collect();
    files.sort(); // deterministic order, as everything here must be

    for path in files {
        let Ok(text) = std::fs::read_to_string(&path) else {
            continue;
        };
        for (line_no, line) in text.lines().enumerate() {
            let Some(at) = line.find(needle) else {
                continue;
            };
            let stated: String = line[at + needle.len()..]
                .chars()
                .take_while(|c| c.is_ascii_digit() || *c == '.')
                .collect();
            let stated = stated.trim_end_matches('.');
            if stated.is_empty() {
                continue;
            }
            checked += 1;
            assert_eq!(
                stated,
                full,
                "{}:{} says the tree is {stated:?}, and it is {full:?}",
                path.file_name().unwrap_or_default().to_string_lossy(),
                line_no + 1
            );
        }
    }

    assert!(
        checked > 0,
        "no agent states what version the tree is — if that sentence was reworded, update this \
         test rather than letting it pass on finding nothing"
    );
}
