"""The binding's own tests, run by CI the way a consumer would run them.

Every number here is compared against something computed in this file rather than read off the
simulation, which is the rule the Rust side follows and the only kind of check worth having: a
value read back from the thing under test agrees with itself no matter what it does.

Run with `python -m pytest`, or as a script.
"""

import dualis


def test_units_crossed_intact():
    assert dualis.one_joule() == 1.0
    assert dualis.__version__ == "0.2.0"


def coupled():
    """A heater feeding a bar, run to four seconds."""
    sim = dualis.Simulation(schedule="multirate", conservation_tolerance=1e-9)
    sim.add_heater("element", watts=2.0, reserve_j=6.0)
    sim.add_bar("bar", length_m=0.020, cells=41, area_m2=1e-4, initial_k=293.15)
    substeps = None
    for _ in range(8):
        substeps = sim.advance(0.5)
    return sim, substeps


def test_every_joule_arrives_and_the_rise_matches_the_closed_form():
    sim, substeps = coupled()

    # 6 J into 20 mm x 1 cm^2 of aluminium. The constant comes from the module, the arithmetic
    # from here, so this is not the simulation checked against itself.
    capacity = dualis.aluminium_heat_capacity_j_per_k(0.020 * 1e-4)
    want = 6.0 / capacity
    got = sim.temperature("bar") - 293.15
    assert abs(got / want - 1.0) < 1e-9, f"wanted {want} K, got {got}"

    # The tank empties and all of it crosses.
    assert sim.reserve_j("element") == 0.0
    assert abs(sim.absorbed_j("bar") - 6.0) < 1e-12

    # The bar subcycles hard and the heater does not, which is what `multirate` is for.
    counts = dict(substeps)
    assert counts["element"] == 1
    assert counts["bar"] > 100


def test_the_fed_end_is_warmer_because_lumped_heat_has_no_place():
    sim, _ = coupled()
    p = sim.profile("bar")
    assert len(p) == 41
    # Heat arriving on a plain channel goes into the first cell, and four seconds of conduction
    # have not levelled it. That gradient is the physics.
    assert p[0] > p[-1]
    # And the profile's mean is the mean, since these are the cells and not a sampled field.
    assert abs(sum(p) / len(p) - sim.temperature("bar")) < 1e-12


def test_the_audit_is_an_exception_and_the_clock_does_not_move():
    # A heater with nobody to take its joules: they leave one domain and arrive nowhere.
    sim = dualis.Simulation(schedule="staggered")
    sim.add_heater("lonely", watts=5.0, reserve_j=100.0)
    try:
        sim.advance(1.0)
        raise AssertionError("the audit should have refused this")
    except dualis.Violation as v:
        assert v.quantity == "energy"
        assert "not consumed" in v.site
        assert v.before == 5.0
        assert v.after == 0.0
        assert v.tolerance > 0.0
    # A refused step leaves the simulation exactly where it was.
    assert sim.time == 0.0


def test_the_mistakes_a_caller_will_make_are_named():
    sim, _ = coupled()
    for call, fragment in [
        (lambda: dualis.Simulation(schedule="magic"), "unknown schedule"),
        (lambda: sim.add_bar("bar", length_m=0.02, cells=41, area_m2=1e-4), "already a domain"),
        (lambda: sim.temperature("nope"), "no domain called"),
        (lambda: sim.add_bar("b", length_m=0.02, cells=1, area_m2=1e-4), "at least two cells"),
        (lambda: sim.add_bar("c", length_m=-1.0, cells=4, area_m2=1e-4), "must be positive"),
    ]:
        try:
            call()
            raise AssertionError(f"{fragment!r} should have been refused")
        except ValueError as e:
            assert fragment in str(e), f"{fragment!r} not in {e}"


def test_the_ledger_is_readable():
    sim, _ = coupled()
    book = dict(sim.ledger())
    # The heater's tank is empty and the bar holds what crossed, so the total is the 6 J that
    # started in the tank.
    assert abs(book["energy"] - 6.0) < 1e-9
    assert sim.domains == ["element", "bar"]


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_"):
            fn()
            print(f"ok  {name}")
    print(chr(10) + "all python tests pass")
