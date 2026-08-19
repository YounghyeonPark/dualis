//! Counts spelled out in prose, against the thing they count.
//!
//! `friction_counts.rs` checks `FRICTION.md`'s own summary sentence against its headings, and that count
//! has not been wrong since. What it does not check is the **six other places** the same total is
//! restated, and a `prose-auditor` pass before the 0.14.0 release found four of them saying twenty-three
//! when the file said twenty-four. Three more documents said fourteen or eighteen scenes when there were
//! twenty-one.
//!
//! Six of that audit's twenty findings were this one shape: a number of scenes or findings written into a
//! document nothing compares against anything. That is a *class*, and a correction closes one instance of
//! it while a test closes the class.
//!
//! # Exact phrases, and why the first design was worse
//!
//! The first version scanned for any number-word near a subject word — "N scenes", "N findings" — and
//! compared every one against the total. It fired three times immediately and **all three were its own
//! fault**: `scenes/README.md` legitimately says "the first of three scenes with nothing to draw", and
//! `AGENTS.md` legitimately says "twelve worked problems" about fifteen example files, because three of
//! them are a quickstart, a benchmark and a README checker rather than worked problems. Subset counts and
//! differently-defined counts are normal prose and a blunt scan cannot tell them from staleness.
//!
//! I had also predicted the wrong failure mode for that design — written down as "the risk is a false
//! pass" — and what it produced was three false failures.
//!
//! So each claim is registered as a **template** with the number left as `{}`, and the test asserts two
//! things: the template filled with the right word is present, and filled with any *other* number-word it
//! is not. That has no false positives, and rewording a sentence makes the test fail with "not found"
//! rather than pass silently — which is the property that matters, because a vacuous pass is how this
//! class of staleness survives in the first place.

use std::path::{Path, PathBuf};

/// Words for the numbers these documents use. Prose does not say "21".
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

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(Path::parent)
        .expect("crates/dualis-world has two ancestors")
        .to_path_buf()
}

/// Assert that `template` filled with `want` appears in `relative`, and that no other number-word fills
/// it.
///
/// A missing file is skipped: this crate is `publish = false` but its tests run from a checkout, and a
/// packaging arrangement that hid a document is not this test's business. A missing *phrase* is a
/// failure, which is the whole point.
fn phrase(relative: &str, template: &str, want: usize) {
    let Ok(text) = std::fs::read_to_string(root().join(relative)) else {
        return;
    };
    assert!(
        template.contains("{}"),
        "the template needs a hole for the number"
    );
    // Case-insensitively, because a claim at the start of a sentence is capitalised and that is a
    // property of English rather than of the count. Two of these phrases are, and the first run of this
    // test failed on both.
    let text = text.to_lowercase();
    let filled = |word: &str| template.to_lowercase().replace("{}", word);
    let wanted = filled(WORDS[want]);
    assert!(
        text.contains(&wanted),
        "{relative} no longer contains {wanted:?}. If the sentence was reworded, update the \
         template here — a phrase that cannot be found is the failure this test is shaped to give \
         instead of a silent pass"
    );
    for (n, word) in WORDS.iter().enumerate() {
        if n == want {
            continue;
        }
        let wrong = filled(word);
        // A shorter number-word is a *substring* of a longer one: "twenty-one worlds described as data"
        // contains "one worlds described as data", and "twenty-four findings" contains "four findings".
        // So a hit only counts when what precedes it is not part of a word — which is one more false
        // positive this design produced before it stopped producing them.
        let spurious = |at: usize| {
            text[..at]
                .chars()
                .next_back()
                .is_some_and(|c| c == '-' || c.is_alphanumeric())
        };
        let mut from = 0;
        while let Some(at) = text[from..].find(&wrong) {
            let hit = from + at;
            assert!(
                spurious(hit),
                "{relative} contains {wrong:?}; the count is {want}"
            );
            from = hit + 1;
        }
    }
}

/// How many findings `FRICTION.md` holds, and how many are fixed — computed the same way
/// `friction_counts.rs` computes them, so the two tests cannot disagree about the number.
fn friction_totals() -> Option<(usize, usize)> {
    let text = std::fs::read_to_string(root().join("crates/dualis-world/FRICTION.md")).ok()?;
    let findings = text
        .lines()
        .filter(|l| {
            l.strip_prefix("## ")
                .and_then(|r| r.split_once('.'))
                .is_some_and(|(n, _)| !n.is_empty() && n.chars().all(|c| c.is_ascii_digit()))
        })
        .count();
    let fixed = text.lines().filter(|l| l.starts_with("**Fixed")).count();
    (findings > 0 && fixed > 0).then_some((findings, fixed))
}

/// **The findings total is the same in all seven places it is written.**
///
/// Four of these were stale at once before the 0.14.0 release — `CLAUDE.md`, two agent files and
/// `FRICTION.md`'s own closing section — while `friction_counts.rs` passed, because it checks one
/// sentence in one file. One file under test is not enough when seven restate the number.
#[test]
fn the_findings_total_agrees_everywhere_it_is_written() {
    let Some((findings, fixed)) = friction_totals() else {
        return;
    };
    println!("  {findings} findings, {fixed} fixed");

    phrase(
        "crates/dualis-world/FRICTION.md",
        "**{} of the twenty-four are fixed**",
        fixed,
    );
    phrase(
        "crates/dualis-world/FRICTION.md",
        "{} findings, and the source has shifted",
        findings,
    );
    phrase(
        "CLAUDE.md",
        "{} findings from using the SDK as a stranger",
        findings,
    );
    phrase("README.md", "**{} findings, eighteen fixed", findings);
    phrase("README.md", "findings, {} fixed and six argued down", fixed);
    phrase(
        "CHANGELOG.md",
        "has already found {} places it is awkward",
        findings,
    );
    phrase(
        "CHANGELOG.md",
        "places it is awkward, {} of which have been changed",
        fixed,
    );
    phrase(
        "crates/dualis-world/src/lib.rs",
        "beside this crate. {} of the twenty-four are",
        fixed,
    );
    phrase(
        ".claude/agents/README.md",
        "by a distance: {} findings",
        findings,
    );
    phrase(".claude/agents/README.md", "findings, {} fixed", fixed);
    phrase(
        ".claude/agents/consumer-advocate.md",
        "Six of the {} findings were recorded",
        findings,
    );
}

/// **The scene count is the same in all five places it is written.**
///
/// `scene.rs` reads the directory rather than a list, so adding a scene never breaks anything — which is
/// exactly why the sentences about the count drift. Two agent files said fourteen and `scenes/README.md`
/// said eighteen in one place while saying twenty-one in four others.
///
/// The subset claims in the same file — "eleven of these twenty-one have a domain with no field",
/// "twelve of the twenty-one have geometry to export" — are covered here only for the *total* half of
/// each sentence. Their own numerators are not counted by anything and are left alone rather than
/// guarded badly.
#[test]
fn the_scene_count_agrees_everywhere_it_is_written() {
    let dir = root().join("crates/dualis-world/scenes");
    let Ok(entries) = std::fs::read_dir(&dir) else {
        return;
    };
    let scenes = entries
        .filter_map(Result::ok)
        .filter(|e| e.path().extension().is_some_and(|x| x == "json"))
        .count();
    if scenes == 0 {
        return;
    }
    println!("  {scenes} scenes");

    let readme = "crates/dualis-world/scenes/README.md";
    phrase(readme, "{} worlds described as data", scenes);
    phrase(readme, "which is what all {} here", scenes);
    // These two sentences also carry a *numerator* — how many scenes have a fieldless domain, how many
    // export geometry — and the templates deliberately cover only the total. A numerator is a different
    // count with a different source, and folding it in here would make the guard fail every time one
    // moved, which is what it did when scene 22 arrived: the template held "Eleven of these {} scenes" and
    // twelve became correct.
    //
    // Both were verified against the binary rather than reasoned about — 12 scenes print a "not drawn"
    // domain and 13 export a `.gltf` — but nothing counts them continuously, and a guard that demanded a
    // template edit for every such move would be a guard somebody rewrites rather than reads.
    phrase(readme, "of these {} scenes", scenes);
    phrase(readme, "of the {} scenes have geometry", scenes);
    phrase(readme, "runs all {} on every commit", scenes);
    phrase(
        "README.md",
        "with {} scenes across all eleven domains",
        scenes,
    );
    phrase(
        ".claude/agents/domain-builder.md",
        "{} ship, all run by CI",
        scenes,
    );
}

/// **The example count in `RELEASING.md`'s own counting command is the number of example files.**
///
/// Not a prose claim — a shell snippet, which is worse, because it looks like a measurement. The command
/// there globbed `crates/dualis/examples/*.rs` and missed the one example that lives in `dualis-optics`,
/// so a release following it would have counted fourteen of fifteen and had no way to know.
///
/// **`AGENTS.md`'s "twelve worked problems" is deliberately not checked against this.** Three of the
/// fifteen files are a quickstart, a benchmark and a README checker, and calling them worked problems
/// would be the wrong sentence rather than the wrong number. A test that insisted they agree would be
/// asserting a definition, and it fired on exactly that when this file's first draft tried.
#[test]
fn the_releasing_example_command_counts_every_example() {
    let mut examples = 0;
    let mut stack = vec![root().join("crates")];
    while let Some(dir) = stack.pop() {
        let Ok(entries) = std::fs::read_dir(&dir) else {
            continue;
        };
        for entry in entries.filter_map(Result::ok) {
            let path = entry.path();
            if path.is_dir() {
                stack.push(path);
            } else if path.extension().is_some_and(|x| x == "rs")
                && path.parent().is_some_and(|p| p.ends_with("examples"))
            {
                examples += 1;
            }
        }
    }
    if examples == 0 {
        return;
    }
    println!("  {examples} example files across the workspace");
    let Ok(text) = std::fs::read_to_string(root().join("RELEASING.md")) else {
        return;
    };
    let want = format!("# examples: {examples}");
    assert!(
        text.contains(&want),
        "RELEASING.md's counting command should say {want:?}; there are {examples} example files \
         and the answer belongs beside the command that finds them"
    );
    assert!(
        !text.contains("ls crates/dualis/examples/*.rs | grep -vc common"),
        "the old glob missed the example in dualis-optics"
    );
}
