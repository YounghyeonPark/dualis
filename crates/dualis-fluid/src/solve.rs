//! The two halves of a projection step: advance, then make it divergence-free.

use crate::channel::Channel;

impl Channel {
    /// Advection, diffusion and the body force, all at once, into a provisional velocity.
    ///
    /// # Flux form, and what it buys
    ///
    /// The advection is written as `∇·(uu)` rather than `u·∇u`, which are the same thing for a
    /// divergence-free field and are not the same discretisation. In flux form every face's
    /// contribution appears twice with opposite signs, so **total momentum changes only by what
    /// the boundaries and the force do** — exactly, not to a tolerance. In the convective form it
    /// changes by whatever the truncation error happens to be.
    ///
    /// Central differences, which are second order and add no dissipation of their own. The price
    /// is [`CELL_REYNOLDS_LIMIT`](crate::channel::CELL_REYNOLDS_LIMIT): with nothing damping what
    /// advection sharpens, a mesh too coarse for the viscosity goes unstable, and no time step
    /// rescues it. Upwinding would trade that for a numerical viscosity that is often larger than
    /// the real one — which is how a scheme comes to report a Reynolds number it is not running at.
    pub(crate) fn advance(&mut self, dt: f64) {
        let (nx, ny, nz) = self.counts();
        let h = self.dx();
        let nu = self.fluid().kinematic_viscosity.to_si();
        let force = self.force();
        let (mut du, mut dv, mut dw) = (
            vec![0.0; self.u_ref().len()],
            vec![0.0; self.v_ref().len()],
            vec![0.0; self.w_ref().len()],
        );

        // --- u, on the x faces -------------------------------------------------------------
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let uc = self.u_at_wrapped(i, j as isize, k);
                    // d(uu)/dx: the product lives at the cell centres either side.
                    let ue = 0.5 * (uc + self.u_at_wrapped((i + 1) % nx, j as isize, k));
                    let uw = 0.5 * (self.u_at_wrapped((i + nx - 1) % nx, j as isize, k) + uc);
                    let duudx = (ue * ue - uw * uw) / h;

                    // d(uv)/dy: `v` interpolated in x onto the u column, `u` interpolated in y.
                    let v_up = 0.5 * (self.v_of(i, j + 1, k) + self.v_of((i + nx - 1) % nx, j + 1, k));
                    let v_dn = 0.5 * (self.v_of(i, j, k) + self.v_of((i + nx - 1) % nx, j, k));
                    let u_up = 0.5 * (uc + self.u_at_wrapped(i, j as isize + 1, k));
                    let u_dn = 0.5 * (self.u_at_wrapped(i, j as isize - 1, k) + uc);
                    let duvdy = (u_up * v_up - u_dn * v_dn) / h;

                    // d(uw)/dz, the same again with z.
                    let w_up = 0.5
                        * (self.w_of(i, j, (k + 1) % nz) + self.w_of((i + nx - 1) % nx, j, (k + 1) % nz));
                    let w_dn = 0.5 * (self.w_of(i, j, k) + self.w_of((i + nx - 1) % nx, j, k));
                    let u_f = 0.5 * (uc + self.u_at_wrapped(i, j as isize, (k + 1) % nz));
                    let u_b = 0.5 * (self.u_at_wrapped(i, j as isize, (k + nz - 1) % nz) + uc);
                    let duwdz = (u_f * w_up - u_b * w_dn) / h;

                    let lap = (self.u_at_wrapped((i + 1) % nx, j as isize, k)
                        + self.u_at_wrapped((i + nx - 1) % nx, j as isize, k)
                        + self.u_at_wrapped(i, j as isize + 1, k)
                        + self.u_at_wrapped(i, j as isize - 1, k)
                        + self.u_at_wrapped(i, j as isize, (k + 1) % nz)
                        + self.u_at_wrapped(i, j as isize, (k + nz - 1) % nz)
                        - 6.0 * uc)
                        / (h * h);

                    du[self.iu_pub(i, j, k)] =
                        dt * (-(duudx + duvdy + duwdz) + nu * lap + force.x);
                }
            }
        }

        // --- w, on the z faces, by the same construction rotated -----------------------------
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let wc = self.w_of(i, j, k);
                    let wf = 0.5 * (wc + self.w_of(i, j, (k + 1) % nz));
                    let wb = 0.5 * (self.w_of(i, j, (k + nz - 1) % nz) + wc);
                    let dwwdz = (wf * wf - wb * wb) / h;

                    let u_e = 0.5
                        * (self.u_at_wrapped((i + 1) % nx, j as isize, k)
                            + self.u_at_wrapped((i + 1) % nx, j as isize, (k + nz - 1) % nz));
                    let u_w = 0.5
                        * (self.u_at_wrapped(i, j as isize, k)
                            + self.u_at_wrapped(i, j as isize, (k + nz - 1) % nz));
                    let w_e = 0.5 * (wc + self.w_of((i + 1) % nx, j, k));
                    let w_w = 0.5 * (self.w_of((i + nx - 1) % nx, j, k) + wc);
                    let dwudx = (w_e * u_e - w_w * u_w) / h;

                    let v_up = 0.5 * (self.v_of(i, j + 1, k) + self.v_of(i, j + 1, (k + nz - 1) % nz));
                    let v_dn = 0.5 * (self.v_of(i, j, k) + self.v_of(i, j, (k + nz - 1) % nz));
                    let w_up = 0.5 * (wc + self.w_at_wrapped(i, j as isize + 1, k));
                    let w_dn = 0.5 * (self.w_at_wrapped(i, j as isize - 1, k) + wc);
                    let dwvdy = (w_up * v_up - w_dn * v_dn) / h;

                    let lap = (self.w_of((i + 1) % nx, j, k)
                        + self.w_of((i + nx - 1) % nx, j, k)
                        + self.w_at_wrapped(i, j as isize + 1, k)
                        + self.w_at_wrapped(i, j as isize - 1, k)
                        + self.w_of(i, j, (k + 1) % nz)
                        + self.w_of(i, j, (k + nz - 1) % nz)
                        - 6.0 * wc)
                        / (h * h);

                    dw[self.iw_pub(i, j, k)] =
                        dt * (-(dwwdz + dwudx + dwvdy) + nu * lap + force.z);
                }
            }
        }

        // --- v, on the y faces. Only the interior ones move when there are walls. -------------
        let (lo, hi) = self.v_interior();
        for k in 0..nz {
            for j in lo..hi {
                for i in 0..nx {
                    let vc = self.v_of(i, j, k);
                    let vu = 0.5 * (vc + self.v_wrapped(i, j as isize + 1, k));
                    let vd = 0.5 * (self.v_wrapped(i, j as isize - 1, k) + vc);
                    let dvvdy = (vu * vu - vd * vd) / h;

                    let u_e = 0.5
                        * (self.u_at_wrapped((i + 1) % nx, j as isize, k)
                            + self.u_at_wrapped((i + 1) % nx, j as isize - 1, k));
                    let u_w = 0.5
                        * (self.u_at_wrapped(i, j as isize, k)
                            + self.u_at_wrapped(i, j as isize - 1, k));
                    let v_e = 0.5 * (vc + self.v_of((i + 1) % nx, j, k));
                    let v_w = 0.5 * (self.v_of((i + nx - 1) % nx, j, k) + vc);
                    let dvudx = (v_e * u_e - v_w * u_w) / h;

                    let w_f = 0.5
                        * (self.w_of(i, j, (k + 1) % nz)
                            + self.w_at_wrapped(i, j as isize - 1, (k + 1) % nz));
                    let w_b = 0.5 * (self.w_of(i, j, k) + self.w_at_wrapped(i, j as isize - 1, k));
                    let v_f = 0.5 * (vc + self.v_of(i, j, (k + 1) % nz));
                    let v_b = 0.5 * (self.v_of(i, j, (k + nz - 1) % nz) + vc);
                    let dvwdz = (v_f * w_f - v_b * w_b) / h;

                    let lap = (self.v_of((i + 1) % nx, j, k)
                        + self.v_of((i + nx - 1) % nx, j, k)
                        + self.v_wrapped(i, j as isize + 1, k)
                        + self.v_wrapped(i, j as isize - 1, k)
                        + self.v_of(i, j, (k + 1) % nz)
                        + self.v_of(i, j, (k + nz - 1) % nz)
                        - 6.0 * vc)
                        / (h * h);

                    dv[self.iv_pub(i, j, k)] =
                        dt * (-(dvvdy + dvudx + dvwdz) + nu * lap + force.y);
                }
            }
        }

        self.add_increments(&du, &dv, &dw);
    }
}
