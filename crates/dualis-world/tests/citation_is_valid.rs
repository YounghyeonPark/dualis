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

/// **The licence is ONE identifier — not an expression, and not a list.**
///
/// This line has failed a deposition twice and the second time was partly this test's fault.
///
/// v0.13.0 put the SPDX **expression** `MIT OR Apache-2.0` here. CFF takes an identifier or a list of
/// them and an expression is neither, so Zenodo rejected it. The fix was a two-element **list**, which
/// validates clean against the CFF 1.2.0 schema — `cffconvert --validate` says so — and Zenodo rejected
/// that too, because Zenodo stores exactly one licence per record and a list is not one.
///
/// And this test asserted the list. It was written to lock in the first fix and it locked in the
/// specification of the second failure instead, which is the shape this workspace has named before:
/// **a test that encodes what the implementation does rather than what the consumer needs.** Validating
/// against the schema was never the requirement; being acceptable to Zenodo was, and the two are
/// different.
///
/// So: one identifier, on the same line as the key. The dual licence lives in `Cargo.toml`, in the two
/// LICENSE files, and in prose in `message` — checked below, because a record that names one of two
/// licences should say that is what it is doing.
#[test]
fn the_licence_is_one_identifier_and_not_an_expression_or_a_list() {
    let text = citation();
    assert!(
        !text.contains("license: MIT OR Apache-2.0"),
        "an SPDX expression is what failed the v0.13.0 deposition"
    );

    let key = text
        .lines()
        .find(|l| l.starts_with("license:"))
        .expect("there is a license key");
    let value = key
        .strip_prefix("license:")
        .expect("just matched it")
        .trim()
        .trim_matches('"');
    assert!(
        !value.is_empty(),
        "`license:` with nothing beside it opens a list, and a list is what failed the v0.14.0          deposition — Zenodo stores one licence"
    );
    assert!(
        ["MIT", "Apache-2.0"].contains(&value),
        "the licence should be one of the project's two, is {value:?}"
    );
    // No list entries anywhere near it, which is how the previous shape would come back.
    assert!(
        !text.contains(
            "
  - MIT"
        ) && !text.contains(
            "
  - Apache-2.0"
        ),
        "a licence list is what Zenodo rejected at v0.14.0"
    );

    // And the dual licence is stated in prose, so naming one is a documented choice.
    assert!(
        text.contains("MIT OR Apache-2.0"),
        "`message` should say what the real licensing is, since `license` can only name one"
    );
}

/// **`.zenodo.json` exists, is the metadata Zenodo asks for, and agrees with `CITATION.cff`.**
///
/// Zenodo prefers `.zenodo.json` over `CITATION.cff` when both are present, so it is the file that
/// decides — and it is written in Zenodo's own vocabulary rather than in CFF's, which means the licence
/// identifier is **lowercase**: `apache-2.0`, from `zenodo.org/api/vocabularies/licenses`. That was asked
/// of the API rather than guessed, because guessing it wrong would have been the third failed deposition
/// on the same field.
///
/// Both files are kept correct rather than relying on the precedence, since "prefers" is documented
/// behaviour rather than a guarantee and the cost of keeping both right is one test.
#[test]
fn the_zenodo_metadata_is_present_and_agrees_with_the_citation() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent().map(|p| p.to_path_buf()))
        .expect("two ancestors");
    let path = root.join(".zenodo.json");
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("reading {}: {e}", path.display()));
    let json: serde_json::Value =
        serde_json::from_str(&text).expect(".zenodo.json is what Zenodo parses; it must be JSON");

    // The licence is a plain string, which is the whole point of this file existing.
    let licence = json["license"]
        .as_str()
        .expect("`license` must be a string, not a list or an object — that is what failed twice");
    assert_eq!(
        licence,
        licence.to_lowercase(),
        "Zenodo's licence vocabulary is lowercase; {licence:?} is not"
    );
    assert_eq!(
        licence, "apache-2.0",
        "the identifier from Zenodo's own vocabulary"
    );

    assert_eq!(
        json["upload_type"].as_str(),
        Some("software"),
        "without this the record is a dataset or a paper"
    );
    assert_eq!(
        json["version"].as_str(),
        Some(env!("CARGO_PKG_VERSION")),
        "the deposition's version must be the crate's"
    );
    assert_eq!(json["access_right"].as_str(), Some("open"));

    let creators = json["creators"]
        .as_array()
        .expect("a record with no creator is worse than a failure, because it succeeds");
    assert!(!creators.is_empty(), "at least one creator");
    for c in creators {
        assert!(
            c["name"].as_str().is_some_and(|n| n.contains(',')),
            "Zenodo wants `Family, Given`, got {:?}",
            c["name"]
        );
        let orcid = c["orcid"]
            .as_str()
            .expect("the ORCID is what disambiguates the author");
        assert!(
            !orcid.starts_with("http"),
            "`.zenodo.json` takes the bare identifier, not the URL: {orcid:?}"
        );
    }

    // And it says the same thing the citation file does, case aside.
    let cff = citation();
    let cff_licence = cff
        .lines()
        .find(|l| l.starts_with("license:"))
        .and_then(|l| l.strip_prefix("license:"))
        .map(|v| v.trim().trim_matches('"').to_lowercase())
        .expect("the citation names a licence");
    assert_eq!(
        cff_licence, licence,
        "the two files must not disagree about the licence, because either one may be the one read"
    );
    assert!(
        text.contains("MIT OR Apache-2.0"),
        "the description should say what the real licensing is, as `CITATION.cff`'s message does"
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
