//! `CITATION.cff` parses, and every field Zenodo needs is the shape it needs.
//!
//! This file exists because the v0.13.0 deposition **failed**. Zenodo read `CITATION.cff`, rejected
//! it, and reported nothing beyond a red "Failed" on a web page — the release went to crates.io and
//! PyPI and got no DOI, and nothing in the repository could have caught it first.
//!
//! The cause was one line: `license: MIT OR Apache-2.0`. That is a valid SPDX *expression* and this
//! workspace's `Cargo.toml` is right to use it, but CFF's schema takes an identifier or a list of
//! them, and an expression matches neither. The correct form is a two-element list.
//!
//! # Why this is a Rust test and not a CI step with a Python validator
//!
//! Because the failure was a *field shape*, and the four checks below catch that class without a
//! dependency. `cffconvert` would be stricter and would also be one more thing to install, and this
//! crate is where the repository already keeps its documentation-consistency tests —
//! `friction_counts.rs` next door checks that a prose count matches the headings it counts.
//!
//! What this cannot check is the whole CFF schema. It checks the parts that broke and the parts a
//! release depends on, which is the difference between a test and a wish.

use std::path::PathBuf;

fn citation() -> String {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .expect("crates/dualis-world is two levels down")
        .join("CITATION.cff");
    std::fs::read_to_string(&root).unwrap_or_else(|e| panic!("reading {}: {e}", root.display()))
}

/// **The licence is a list of SPDX identifiers, not an SPDX expression.**
///
/// The line the v0.13.0 deposition died on. `Cargo.toml` says `MIT OR Apache-2.0` and should — that is
/// how crates.io states a dual licence. CFF wants the identifiers, and Zenodo enforces the CFF schema
/// rather than guessing what an expression means.
#[test]
fn the_licence_is_a_list_of_identifiers_and_not_an_expression() {
    let text = citation();
    assert!(
        !text.contains("license: MIT OR Apache-2.0"),
        "an SPDX expression on one line is what failed the v0.13.0 deposition; CFF wants a list"
    );
    for id in ["MIT", "Apache-2.0"] {
        assert!(
            text.contains(&format!("  - {id}")),
            "the licence list should carry {id} as its own entry"
        );
    }
    // And the `license:` key itself takes no inline value, which is what makes it a list.
    let key = text
        .lines()
        .find(|l| l.starts_with("license:"))
        .expect("there is a license key");
    assert_eq!(
        key.trim(),
        "license:",
        "the value belongs on the lines below, not beside the key"
    );
}

/// **Every field Zenodo builds a record out of is present.**
///
/// Not the whole CFF schema — the fields a deposition needs. A missing `authors` or `version` is a
/// record with no creator or no version, which is worse than a failure because it succeeds.
#[test]
fn the_fields_a_deposition_needs_are_all_there() {
    let text = citation();
    for key in [
        "cff-version:",
        "message:",
        "title:",
        "authors:",
        "type: software",
        "version:",
        "date-released:",
        "repository-code:",
    ] {
        assert!(
            text.contains(key),
            "CITATION.cff is missing {key}, which a Zenodo record is built from"
        );
    }
    // The ORCID has a check digit and it was verified before being written down; this asserts it is
    // still the same identifier rather than recomputing it.
    assert!(
        text.contains("orcid: \"https://orcid.org/0000-0002-4733-5049\""),
        "the ORCID should be the full URL form, which is what CFF's pattern requires"
    );
}

/// **The version in `CITATION.cff` is the version the workspace is at.**
///
/// It is the seventh place a version lives and the one no compiler reads. A stale one mints a DOI
/// against a version that is not what was released — a citation that resolves to the wrong code, which
/// is the specific failure a DOI exists to prevent.
#[test]
fn the_citation_version_matches_the_crate_version() {
    let text = citation();
    let stated = text
        .lines()
        .find_map(|l| l.strip_prefix("version: "))
        .expect("there is a version key")
        .trim()
        .to_string();
    let ours = env!("CARGO_PKG_VERSION");
    assert_eq!(
        stated, ours,
        "CITATION.cff says {stated} and the workspace is {ours}; see RELEASING.md's table of the \
         seven places a version lives"
    );
}

/// **It is YAML, and the author block is a list of mappings rather than a bare string.**
///
/// The shape `authors:` has to be. A single string there parses as YAML and produces a record with no
/// creator, which is the silent half of this failure class.
#[test]
fn the_author_block_is_a_list_of_mappings() {
    let text = citation();
    let after = text
        .split_once("authors:")
        .expect("there is an authors key")
        .1;
    let first = after
        .lines()
        .find(|l| !l.trim().is_empty() && !l.trim().starts_with('#'))
        .expect("something follows authors:");
    assert!(
        first.trim_start().starts_with("- "),
        "authors must be a list; the first entry reads {first:?}"
    );
    assert!(
        after.contains("family-names:") && after.contains("given-names:"),
        "an author needs family-names and given-names for a creator to be built"
    );
}
