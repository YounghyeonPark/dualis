//! A grid had no void, so every cell a part did not occupy was some other substance.
//!
//! `ARCHITECTURE.md` names it plainly: "a part in a box is surrounded by another material rather
//! than by nothing. Insulating it is a substance with a low conductivity, which is not the same
//! thing." A low conductivity still conducts, still stores heat and still sets a stability limit.
//! Nothing does none of those, and until this existed an assembly of two parts in air was two
//! parts buried in whatever the block was made of — which is the shape that makes assembly wrong
//! rather than merely coarse.

use pantometry_core::conserved::quantity;
use pantometry_core::units::{Length, Temperature, Time};
use pantometry_core::{Domain, Exchange, Substance};
use pantometry_thermal::{Solid3D, HEAT};

/// A bar of `n` cells along z, hot at one end.
fn bar(n: usize, hot_cells: usize) -> Solid3D {
    let mut block = Solid3D::new(
        "bar",
        Substance::copper(),
        (1, 1, n),
        Length::mm(5.0),
        Temperature::celsius(20.0),
    );
    for k in 0..hot_cells {
        block.set_temperature(0, 0, k, Temperature::celsius(300.0));
    }
    block
}

fn march(block: &mut Solid3D, steps: usize) {
    let dt = block.max_stable_dt(Time::ZERO);
    let mut bus = Exchange::new();
    for _ in 0..steps {
        block.step(Time::ZERO, dt, &mut bus).expect("a stable step");
    }
}

/// **Void breaks the conduction path, and a low-conductivity substance does not.**
///
/// Three bars, identical but for what sits in the middle cell: copper, the poorest insulator the
/// catalogue has, and nothing. The middle case is the one that matters — it is what somebody
/// does when the grid has no void, and it is why "use a bad conductor" is not the same answer.
///
/// **The gap is not zero, and this test learned that by failing.** When void arrived it carried
/// nothing at all and this asserted the far end had not moved by a bit. Then gaps started
/// radiating, which is real physics a vacuum clearance always has, and the assertion became
/// false — correctly. What it pins now is the **ordering**, which is the stronger statement:
/// metal carries the most, a poor conductor carries a fraction of that, and a gap carries only
/// what radiation gives, which is less again. The conductances say so too — for these cells
/// PLA's `kA/dx` is 6.5e-4 W/K against radiation's `4σT³A/(1/ε₁+1/ε₂−1)` of 1.25e-4.
#[test]
fn void_breaks_the_path_and_a_poor_conductor_only_slows_it() {
    let steps = 4000;

    let mut solid = bar(5, 2);
    march(&mut solid, steps);
    let through_copper = solid.temperature_at(0, 0, 4).to_si();

    let mut insulated = bar(5, 2);
    insulated.fill(
        Substance::from_name("pla").expect("the catalogue has it"),
        |_, _, k| k == 2,
    );
    march(&mut insulated, steps);
    let through_plastic = insulated.temperature_at(0, 0, 4).to_si();

    let mut gapped = bar(5, 2).empty(|_, _, k| k == 2);
    march(&mut gapped, steps);
    let through_nothing = gapped.temperature_at(0, 0, 4).to_si();

    let start = Temperature::celsius(20.0).to_si();
    assert!(
        through_copper > start + 50.0,
        "copper should carry the heat: {through_copper:.2} K"
    );
    assert!(
        through_plastic > start + 1.0,
        "a poor conductor still conducts, which is the whole point: {through_plastic:.2} K"
    );
    assert!(
        through_nothing > start,
        "a gap radiates, so the far end is not untouched: {through_nothing:.4} K"
    );
    assert!(
        through_plastic < through_copper && through_nothing < through_plastic,
        "metal, then a poor conductor, then only radiation: {through_copper:.2} /          {through_plastic:.2} / {through_nothing:.2} K"
    );
}

/// **A gap with a non-radiating surface really does carry nothing.**
///
/// The control for the test above: set the emissivity to zero and the parallel-plate series has
/// no conductance in it, so the two parts are as uncoupled as void alone used to make them. It
/// is also the check that the gap's exchange is the *radiative* one and not some residue of the
/// conduction stencil leaking across a cell it should not see.
#[test]
fn a_gap_between_mirrors_carries_nothing() {
    let mut mirror = Substance::copper();
    if let Some(thermal) = mirror.thermal.as_mut() {
        thermal.emissivity = 0.0;
    }
    let mut block = Solid3D::new(
        "bar",
        mirror,
        (1, 1, 5),
        Length::mm(5.0),
        Temperature::celsius(20.0),
    );
    for k in 0..2 {
        block.set_temperature(0, 0, k, Temperature::celsius(300.0));
    }
    let mut gapped = block.empty(|_, _, k| k == 2);
    march(&mut gapped, 4000);

    let start = Temperature::celsius(20.0).to_si();
    assert!(
        (gapped.temperature_at(0, 0, 4).to_si() - start).abs() < 1e-9,
        "two perfect mirrors exchange nothing: {:.6} K",
        gapped.temperature_at(0, 0, 4).to_si()
    );
}

/// **A void cell has no temperature, and says so rather than returning one.**
///
/// A zero or an ambient here is a value somebody would plot, average or believe. Every field
/// reader in this workspace already skips a non-finite sample, so `NaN` is the honest answer that
/// costs nothing downstream — and the averages leave it out rather than being dragged by it.
#[test]
fn a_void_cell_has_no_temperature_and_no_vote() {
    let mut block = Solid3D::new(
        "part",
        Substance::copper(),
        (4, 1, 1),
        Length::mm(5.0),
        Temperature::celsius(100.0),
    );
    block.set_temperature(0, 0, 0, Temperature::celsius(300.0));
    let solid_mean = block.mean_temperature().to_si();

    let gapped = block.empty(|i, _, _| i == 3);
    assert!(gapped.is_void(3, 0, 0) && !gapped.is_void(0, 0, 0));
    assert_eq!(gapped.void_cells(), 1);
    assert!(
        gapped.temperature_at(3, 0, 0).to_si().is_nan(),
        "there is nothing there to have a temperature"
    );

    // The mean is over what is there: dropping one 100 C cell of four raises it.
    let voided_mean = gapped.mean_temperature().to_si();
    assert!(
        voided_mean > solid_mean,
        "the average should be over the part, not the box: {solid_mean:.2} then {voided_mean:.2}"
    );
    // And the extremes ignore it too.
    assert!(gapped.peak_temperature().to_si().is_finite());
    assert!(gapped.coldest_temperature().to_si().is_finite());
}

/// **Void holds no heat, so it takes no share of what arrives on the bus** — and the books still
/// balance, because what the void does not take is exactly what the rest does.
///
/// The failure this prevents is quiet: joules spread over cells that cannot hold them are joules
/// the ledger counts and the block does not have.
#[test]
fn void_takes_no_share_of_the_bus() {
    let solid_then_gapped = |void: bool| {
        let mut block = Solid3D::new(
            "part",
            Substance::copper(),
            (4, 1, 1),
            Length::mm(5.0),
            Temperature::celsius(20.0),
        );
        if void {
            block = block.empty(|i, _, _| i >= 2);
        }
        let opening = block.ledger().get(quantity::ENERGY).unwrap_or(0.0);
        let mut bus = Exchange::new();
        bus.publish(HEAT, 1000.0);
        block
            .step(Time::ZERO, Time::from_si(1e-6), &mut bus)
            .expect("a tiny step");
        let closing = block.ledger().get(quantity::ENERGY).unwrap_or(0.0);
        (block, closing - opening)
    };

    let (whole, gained_whole) = solid_then_gapped(false);
    let (half, gained_half) = solid_then_gapped(true);

    // Every joule arrives in both cases — the bus does not know about void.
    assert!((gained_whole - 1000.0).abs() < 1e-6, "{gained_whole}");
    assert!((gained_half - 1000.0).abs() < 1e-6, "{gained_half}");

    // But half the cells hold them, so the rise is twice as big.
    let rise_whole = whole.mean_temperature().to_si() - Temperature::celsius(20.0).to_si();
    let rise_half = half.mean_temperature().to_si() - Temperature::celsius(20.0).to_si();
    assert!(
        (rise_half / rise_whole - 2.0).abs() < 1e-9,
        "two cells should warm twice as fast as four: {rise_whole:.4} K then {rise_half:.4} K"
    );
}

/// **A block that is entirely void is stable and holds nothing**, rather than dividing by a
/// capacity of zero somewhere. The degenerate case, checked because a caller building an
/// assembly cell by cell passes through it.
#[test]
fn a_block_of_nothing_is_not_a_division_by_zero() {
    let mut block = Solid3D::new(
        "nothing",
        Substance::copper(),
        (2, 2, 2),
        Length::mm(5.0),
        Temperature::celsius(20.0),
    )
    .empty(|_, _, _| true);

    assert_eq!(block.void_cells(), 8);
    assert!(block.max_stable_dt(Time::ZERO).to_si().is_infinite());
    assert_eq!(block.ledger().get(quantity::ENERGY).unwrap_or(0.0), 0.0);

    let mut bus = Exchange::new();
    block
        .step(Time::ZERO, Time::from_si(1.0), &mut bus)
        .expect("a block of nothing steps");
    assert!(block.mean_temperature().to_si().is_nan());
}
