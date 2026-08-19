//! `FRICTION.md` counts its own findings correctly.
//!
//! That file opens by saying how many findings there are and how many are fixed, and both
//! numbers are written in words by whoever last edited it. They have been wrong twice: once
//! claiming ten fixed when eleven were, and once claiming twenty-one findings when the file
//! contained eighteen — that second time because an edit's anchor silently did not match, so the
//! summary moved and the content did not.
//!
//! The second failure is the instructive one. A count that disagrees with the file is not a
//! typo; it is the only visible symptom of an edit that did not happen. So this checks the
//! summary against the headings, which is the thing a reader cannot do at a glance and the
//! author evidently cannot do reliably either.

/// Words for the numbers the summary is likely to use. Written out because the file is prose and
/// prose does not say "15".
const WORDS: [&str; 41] = [
    "zero",
    "one",
    "two",
    "three",
    "four",
    "five",
    "six",
    "seven",
    "eight",
    "nine",
    "ten",
    "eleven",
    "twelve",
    "thirteen",
    "fourteen",
    "fifteen",
    "sixteen",
    "seventeen",
    "eighteen",
    "nineteen",
    "twenty",
    "twenty-one",
    "twenty-two",
    "twenty-three",
    "twenty-four",
    "twenty-five",
    "twenty-six",
    "twenty-seven",
    "twenty-eight",
    "twenty-nine",
    "thirty",
    "thirty-one",
    "thirty-two",
    "thirty-three",
    "thirty-four",
    "thirty-five",
    "thirty-six",
    "thirty-seven",
    "thirty-eight",
    "thirty-nine",
    "forty",
];

fn friction() -> Option<String> {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("FRICTION.md");
    std::fs::read_to_string(path).ok()
}

/// **The summary line's two numbers are the file's two numbers.**
///
/// Findings are `## <n>.` headings; a fixed one opens its resolution with `**Fixed`. The summary
/// reads "**X of the Y are fixed**, and Z are recorded rather than actioned", so all three have
/// to agree with each other as well as with the file.
#[test]
fn the_summary_counts_what_the_file_contains() {
    let Some(text) = friction() else {
        return; // not packaged; nothing to check
    };

    let findings = text
        .lines()
        .filter(|l| {
            l.strip_prefix("## ")
                .and_then(|r| r.split_once('.'))
                .is_some_and(|(n, _)| !n.is_empty() && n.chars().all(|c| c.is_ascii_digit()))
        })
        .count();
    let fixed = text.lines().filter(|l| l.starts_with("**Fixed")).count();

    assert!(
        findings > 0,
        "no `## <n>.` findings found — was the file reformatted?"
    );
    assert!(
        fixed <= findings,
        "{fixed} fixed against {findings} findings"
    );

    // Numbered consecutively from 1, which catches a finding added with the wrong number as well
    // as one whose content never landed.
    let numbers: Vec<usize> = text
        .lines()
        .filter_map(|l| l.strip_prefix("## ")?.split_once('.')?.0.parse().ok())
        .collect();
    let expected: Vec<usize> = (1..=findings).collect();
    assert_eq!(
        numbers, expected,
        "findings are not numbered 1..={findings}"
    );

    let want_total = WORDS.get(findings).copied().unwrap_or("<too many>");
    let want_fixed = WORDS.get(fixed).copied().unwrap_or("<too many>");
    let want_open = WORDS.get(findings - fixed).copied().unwrap_or("<too many>");

    // The three numbers are *parsed* out of the sentence rather than searched for in it.
    //
    // Substring matching is wrong here and quietly so: `"twenty-one".contains("twenty")` is
    // true, so deleting the last finding — which is exactly what a silently-failed edit looks
    // like — left the summary saying "twenty-one" over twenty findings and the check passed.
    // Found by mutating the file, which is the only way a hole of that shape shows up.
    let lower = text.to_lowercase();
    let summary = lower
        .lines()
        .find(|l| l.contains("are fixed**"))
        .expect("the summary sentence is `**X of the Y are fixed**`");

    let stated_fixed = summary
        .split_once("**")
        .and_then(|(_, rest)| rest.split_once(" of the "))
        .map(|(n, _)| n.trim())
        .expect("the summary names a fixed count");
    let stated_total = summary
        .split_once(" of the ")
        .and_then(|(_, rest)| rest.split_once(" are fixed"))
        .map(|(n, _)| n.trim())
        .expect("the summary names a total");
    let stated_open = summary
        .split_once(", and ")
        .and_then(|(_, rest)| rest.split_once(" are recorded"))
        .map(|(n, _)| n.trim())
        .expect("the summary says how many are recorded rather than actioned");

    assert_eq!(
        stated_fixed, want_fixed,
        "{fixed} findings are marked fixed, so the summary should open with {want_fixed:?}"
    );
    assert_eq!(
        stated_total, want_total,
        "the file contains {findings} findings, so the summary should say {want_total:?}"
    );
    assert_eq!(
        stated_open,
        want_open,
        "{findings} findings with {fixed} fixed leaves {}, so the summary should say \
         {want_open:?}",
        findings - fixed
    );
}
