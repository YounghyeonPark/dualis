"""The binding's own tests, run by CI the way a consumer would run them.

Every number here is compared against something computed in this file rather than read off the
simulation, which is the rule the Rust side follows and the only kind of check worth having: a
value read back from the thing under test agrees with itself no matter what it does.

Run with `python -m pytest`, or as a script.
"""

import dualis


def test_units_crossed_intact():
    assert dualis.one_joule() == 1.0

    # `__version__` comes from Rust's CARGO_PKG_VERSION; the wheel's metadata comes from
    # pyproject.toml. They are two hand-maintained files that nothing forces to agree, and a
    # release where they disagree ships a module that reports a version it is not. Comparing
    # them checks that, and — unlike the hardcoded "0.2.0" this replaces — does not need
    # editing every time the version moves, which is the thing that made it a stale literal
    # rather than a check.
    from importlib.metadata import version

    assert dualis.__version__ == version("dualis"), (
        f"Cargo.toml says {dualis.__version__}, pyproject.toml says {version('dualis')}"
    )


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


def motor(watts=6.0, seconds=900):
    """A winding inside a case, run to `seconds`. The shape a lump cannot express."""
    sim = dualis.Simulation(schedule="staggered", conservation_tolerance=1e-9)
    sim.add_heater("losses", watts=watts, reserve_j=watts * seconds + 1.0)
    sim.add_network(
        "motor",
        nodes=[
            {"name": "winding", "material": "copper", "volume_m3": 18e-6,
             "thickness_m": 2e-3, "initial_k": 298.15},
            {"name": "case", "material": "aluminium", "volume_m3": 220e-6,
             "thickness_m": 4e-3, "initial_k": 298.15,
             "ambient_k": 298.15, "area_m2": 0.042},
        ],
        links=[{"from": "winding", "to": "case", "w_per_k": 0.9}],
        absorbing="winding",
    )
    for _ in range(seconds):
        sim.advance(1.0)
    return sim


def test_a_network_carries_the_drop_across_the_joint():
    sim = motor()
    temps = dict(sim.node_temperatures("motor"))
    assert set(temps) == {"winding", "case"}

    drop = temps["winding"] - temps["case"]
    # 6 W across 0.9 W/K is 6.67 K at steady state. Computed here, from the numbers this test
    # passed in -- not read back off the simulation, which would be checking it against itself.
    steady = 6.0 / 0.9
    assert 5.5 < drop < steady, f"drop {drop:.3f} K against a {steady:.3f} K ceiling"

    # The same number by the other route, and they must agree exactly: one reads a handle it
    # resolved itself, the other walks every node.
    assert sim.node_temperature("motor", "winding") == temps["winding"]

    # And the flux is the drop times the conductance, which is the definition of a conductance
    # and the one place the binding could have transposed a pair.
    q = sim.heat_flow_w("motor", "winding", "case")
    assert abs(q - 0.9 * drop) < 1e-9, f"{q:.6f} W against {0.9 * drop:.6f} W"
    # Antisymmetric, so the sign convention is stated rather than assumed.
    assert abs(sim.heat_flow_w("motor", "case", "winding") + q) < 1e-12


def test_a_network_of_one_node_matches_the_lump_it_reduces_to():
    """The reduction the Rust side checks bit-for-bit, verified across the binding boundary.

    If a unit conversion were wrong in `add_network` but right in `add_lump` -- a millimetre
    read as a metre, say -- both would still run and audit green. Only the comparison catches
    it, and it has to be the *same* physical body described through the two calls.
    """
    volume, thickness, area = 1.2e-4, 3e-3, 0.0072
    a = dualis.Simulation(schedule="staggered", conservation_tolerance=1e-9)
    a.add_lump("plate", volume_m3=volume, thickness_m=thickness, area_m2=area,
               initial_k=358.15, ambient_k=293.15)

    b = dualis.Simulation(schedule="staggered", conservation_tolerance=1e-9)
    b.add_network(
        "plate",
        nodes=[{"name": "only", "material": "aluminium", "volume_m3": volume,
                "thickness_m": thickness, "initial_k": 358.15,
                "ambient_k": 293.15, "area_m2": area}],
        links=[],
        absorbing="only",
    )

    # 2000 s against a time constant of about 5800 s (290 J/K over 0.050 W/K of convection
    # plus radiation), so roughly a third of the decay -- enough trajectory that a unit error
    # in one call and not the other could not stay hidden behind a barely-moving value.
    for _ in range(4000):
        a.advance(0.5)
        b.advance(0.5)

    lump = a.temperature("plate")
    node = b.node_temperature("plate", "only")
    assert lump == node, f"lump {lump!r} against network {node!r}"
    # It has to have cooled, or this compared two identical starting values.
    cooled = 358.15 - node
    assert 10.0 < cooled < 65.0, f"the plate cooled {cooled:.3f} K"


def test_the_networks_mistakes_are_refused_by_name():
    sim = dualis.Simulation()
    good = [{"name": "a", "material": "copper", "volume_m3": 1e-5,
             "thickness_m": 1e-3, "initial_k": 300.0}]

    def build(name, nodes, links, absorbing="a"):
        return lambda: sim.add_network(name, nodes=nodes, links=links, absorbing=absorbing)

    half_cooled = [dict(good[0], ambient_k=300.0)]
    bad_material = [dict(good[0], material="unobtainium")]
    missing_key = [{"name": "a", "material": "copper", "thickness_m": 1e-3, "initial_k": 300.0}]
    wrong_type = [dict(good[0], volume_m3="big")]

    for call, fragment in [
        (build("n1", [], []), "at least one node"),
        (build("n2", bad_material, []), "unknown material"),
        (build("n3", missing_key, []), 'missing \"volume_m3\"'),
        (build("n4", wrong_type, []), "volume_m3 should be a number"),
        (build("n5", half_cooled, []), "has ambient_k but not area_m2"),
        (build("n6", good, [{"from": "a", "to": "ghost", "w_per_k": 1.0}]), 'no node named \"ghost\"'),
        (build("n7", good, [], absorbing="ghost"), "which is not a node"),
        (build("n8", good, [{"from": "a", "to": "a", "w_per_k": 1.0}]), "itself"),
    ]:
        try:
            call()
            raise AssertionError(f"{fragment!r} should have been refused")
        except ValueError as e:
            assert fragment in str(e), f"{fragment!r} not in {e}"

    # A network declines to average its own nodes, and says which call to use instead. The
    # whole reason to build one is that its nodes differ, so a mean describes no part of it.
    sim.add_network("ok", nodes=good, links=[], absorbing="a")
    try:
        sim.temperature("ok")
        raise AssertionError("a network should not answer temperature()")
    except ValueError as e:
        assert "node_temperatures" in str(e), str(e)


def test_a_network_that_loses_its_heat_is_still_audited():
    """The audit is the reason to use this from Python rather than numpy, so it has to be live
    on the new domain too. A network that takes joules and books fewer would be caught here --
    and the run below crosses real energy, so there is something to be wrong about."""
    sim = motor(watts=20.0, seconds=300)
    absorbed = sim.absorbed_j("motor")
    assert abs(absorbed - 20.0 * 300) < 1e-6, absorbed
    book = dict(sim.ledger())
    # The heater's remaining tank plus what the network stored and shed. Started at 6001 J.
    assert abs(book["energy"] - (20.0 * 300 + 1.0)) < 1e-6, book


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_"):
            fn()
            print(f"ok  {name}")
    print(chr(10) + "all python tests pass")
