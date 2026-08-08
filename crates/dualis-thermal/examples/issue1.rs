use dualis_core::substance::ThermalProps;
use dualis_core::{Domain, Exchange, Substance};
use dualis_thermal::{Environment, LumpedMass, HEAT};
use dualis_units::{
    Area, Density, Length, Power, SpecificHeat, Temperature, ThermalConductivity, ThermalExpansion,
    Time, Volume,
};

fn part(e: f64) -> LumpedMass {
    let (area, vol) = (0.030_24, 3.456e-4);
    LumpedMass::new(
        "box",
        Substance {
            name: "g".into(),
            density: Density::kg_per_m3(1.122 / vol),
            thermal: Some(ThermalProps {
                conductivity: ThermalConductivity::w_per_m_k(150.0),
                specific_heat: SpecificHeat::j_per_kg_k(600.0),
                expansion: ThermalExpansion::ppm_per_k(23.0),
                emissivity: e,
            }),
            mechanical: None,
            acoustic: None,
        },
        Volume::from_si(vol),
        Length::from_si(vol / area),
        Temperature::celsius(25.0),
        Environment {
            ambient: Temperature::celsius(25.0),
            convection_w_per_m2_k: 7.0,
            area: Area::from_si(area),
        },
    )
}

fn main() {
    println!(
        "{:>5} {:>9} {:>9} {:>7} | {:>9} {:>9} {:>9} {:>7} | max_stable_dt",
        "eps", "quoted dT", "settled", "ratio", "tau cold", "tau hot", "measured", "cold/m"
    );
    for e in [0.05, 0.09, 0.3, 0.5, 0.9, 1.0] {
        let quoted = part(e).equilibrium_rise(Power::w(21.0)).to_si();
        let mut m = part(e);
        let mut bus = Exchange::new();
        for k in 0..600_000 {
            bus.publish(HEAT, 21.0);
            m.step(Time::s(k as f64), Time::s(1.0), &mut bus).unwrap();
        }
        let settled = m.rise().to_si();
        // tau at rest, which is what a caller gets before stepping
        let tau_q = part(e).time_constant().to_si();
        let mut m2 = part(e);
        let mut bus2 = Exchange::new();
        let mut t63 = f64::NAN;
        for k in 0..600_000 {
            bus2.publish(HEAT, 21.0);
            m2.step(Time::s(k as f64), Time::s(1.0), &mut bus2).unwrap();
            if t63.is_nan() && m2.rise().to_si() >= settled * (1.0 - 1.0 / std::f64::consts::E) {
                t63 = k as f64;
            }
        }
        // tau once it is hot: the same call on the settled body.
        let tau_hot = m.time_constant().to_si();
        let dt_cold = part(e).max_stable_dt(Time::ZERO).to_si();
        let dt_hot = m.max_stable_dt(Time::ZERO).to_si();
        println!("{e:>5.2} {quoted:>8.1}K {settled:>8.1}K {:>6.3}x | {:>7.1}min {:>7.1}min {:>7.1}min {:>6.2}x | dt {dt_cold:>5.0}s->{dt_hot:>5.0}s",
                 quoted/settled, tau_q/60.0, tau_hot/60.0, t63/60.0, tau_q/t63);
    }
}
