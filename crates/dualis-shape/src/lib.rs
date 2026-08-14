//! Designed geometry as simulation input: a mesh read, measured, and turned into cells.
//!
//! Every domain in this workspace takes a **structured grid** — a box of cubes, `counts` and a cell
//! size. Every 3D file a person designs is a **surface**: a bag of triangles, or a boundary
//! representation that was tessellated into one. This crate is the bridge, and the bridge is where the
//! physics gets decided.
//!
//! # The cell size decides three things at once, and only one of them is obvious
//!
//! A caller picking `cell_mm` is picking, in one number:
//!
//! - **how much of the shape survives.** A 0.5 mm rib voxelised at 2 mm is not a thin rib, it is *gone*,
//!   and the simulation runs perfectly well without it;
//! - **the stability limit**, and therefore the step count, as `dx²`;
//! - **the discretisation error** of whatever physics runs on it.
//!
//! The first is the one that has no symptom. A missing feature does not make a solver fail, produce a
//! `NaN`, or trip the conservation audit — it produces a smooth, plausible answer about a different
//! object. That is the failure this workspace is organised around not having, so this crate's job is not
//! only to rasterise but to **say what the rasterisation lost**: see [`Loss`].
//!
//! # What this is not
//!
//! **Not a domain.** It depends on `dualis-units` and nothing else in the workspace, and no domain
//! depends on it. It produces a predicate; a domain's `fill` consumes one. That is the whole coupling,
//! and it means adding geometry cost the ten domains nothing.
//!
//! **Not a mesh library.** No refinement, no repair, no boolean operations, no simplification. It reads
//! what a CAD tool exported, measures it, and rasterises it.
//!
//! **STL only, and that is a dependency decision rather than a preference.** STEP needs a B-rep kernel
//! and glTF needs a JSON and binary parser; both are external crates, and this workspace gates every one
//! of its twelve through `deny.toml`. STL is the format every CAD tool exports and it can be read in two
//! hundred lines with nothing added. When a second format earns its dependency it goes beside this one.
//!
//! # Where it sits
//!
//! ```text
//!   a designed file ──► Mesh ──► Voxels ──► `|i, j, k| voxels.contains(i, j, k)`
//!                         │         │                        │
//!                    volume,     Loss:                  Solid3D::fill
//!                    closed?     what the grid          Block::fill
//!                                could not hold         Waves::fill
//! ```
//!
//! The predicate on the right is exactly the signature those three already take, which is why this
//! needed no change to any of them.

#![deny(missing_docs)]

mod mesh;
mod voxels;

pub use mesh::{Mesh, Triangle};
pub use voxels::{Loss, Voxels};
