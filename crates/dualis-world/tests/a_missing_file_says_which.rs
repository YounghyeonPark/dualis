//! A file this binary cannot read is named, and the naming is checked.
//!
//! The first thing a stranger typed at this binary was `dualis-world scene.json` in a directory
//! with no `scene.json`, and it answered
//!
//! ```text
//! Error: Os { code: 2, kind: NotFound, message: "The system cannot find the file specified." }
//! ```
//!
//! which is a complete diagnosis of the wrong question: it says what happened and not what it
//! happened to. Every other error this binary produces had carried its path for months — the
//! parse path even carries `file:line:column` — so the two halves of reading a file disagreed
//! about whether the reader deserved to know the file's name.
//!
//! These tests run the built binary rather than calling a function, because the defect was in
//! the seam between a `?` and `main`'s return type, and nothing below that seam could see it.

use std::process::Command;

/// The binary under test, as cargo built it beside this test.
fn binary() -> std::path::PathBuf {
    let mut p = std::env::current_exe().expect("the test binary knows where it is");
    p.pop(); // the test's own file name
    if p.ends_with("deps") {
        p.pop();
    }
    p.join(format!("dualis-world{}", std::env::consts::EXE_SUFFIX))
}

/// Run the binary with these arguments and give back `(stdout + stderr, success)`.
fn run(args: &[&str]) -> (String, bool) {
    let out = Command::new(binary())
        .args(args)
        .output()
        .expect("the binary runs");
    let mut text = String::from_utf8_lossy(&out.stdout).into_owned();
    text.push_str(&String::from_utf8_lossy(&out.stderr));
    (text, out.status.success())
}

/// **A missing scene is named, and the message says how to get one.**
///
/// "Not found" is a complete diagnosis and an incomplete answer: nothing in this repository is
/// called `scene.json` until somebody writes one, and the two ways to get one are worth the
/// line they cost.
#[test]
fn a_missing_scene_is_named_and_answered() {
    let (text, ok) = run(&["definitely-not-here.json"]);
    assert!(!ok, "a missing file must not exit zero");
    assert!(
        text.contains("definitely-not-here.json"),
        "the error does not name the file: {text}"
    );
    assert!(
        text.contains("--emit-default") && text.contains("scenes/"),
        "the error does not say how to get a scene: {text}"
    );
}

/// **The message is printed, not `Debug`-printed.** Returning `Result` from `main` prints the
/// Debug form, so a `String` error arrives wrapped in quotes with its newlines escaped — the
/// hint above came out as a literal backslash-n on one line. Checked by its absence.
#[test]
fn the_error_is_a_message_rather_than_a_debug_dump() {
    let (text, _) = run(&["definitely-not-here.json"]);
    assert!(
        !text.contains("Os {") && !text.contains("kind: NotFound"),
        "the raw io::Error is showing through: {text}"
    );
    assert!(
        !text.contains(r"\n"),
        "an escaped newline means this went through Debug: {text}"
    );
    // And the hint really is on its own line, which is what makes it readable.
    assert!(
        text.lines().count() >= 2,
        "the hint should be a second line: {text}"
    );
}

/// **Every path the binary reads or writes is named on failure, not just the run path.**
///
/// The fix routed nine call sites through one helper; a tenth added later that goes around it
/// would pass the test above and fail this one, because this one asks the other entry points.
#[test]
fn the_other_entry_points_name_their_files_too() {
    for args in [
        vec!["--check", "definitely-not-here.json"],
        vec!["verify", "definitely-not-here.json"],
    ] {
        let (text, ok) = run(&args);
        assert!(!ok, "{args:?} must not exit zero");
        assert!(
            text.contains("definitely-not-here.json"),
            "{args:?} does not name the file: {text}"
        );
    }
}
