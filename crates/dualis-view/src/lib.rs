//! Drawing a dualis run, without the run knowing it is being drawn.
//!
//! The top layer. [`dualis-scene`](https://docs.rs/dualis-scene) says what one instant of a run
//! looks like — a [`Frame`](dualis_scene::Frame) of panels and readings — and this turns that
//! into something a person can look at.
//!
//! # The view is chosen by the shape of the data
//!
//! Not by the name of the domain. That rule is what makes this layer finishable at all: a
//! physics that arrives tomorrow gets a correct picture without a line here changing.
//!
//! | what the data is | what it becomes |
//! | --- | --- |
//! | scalars over time | a line chart, one series per reading |
//! | a 1D field | a profile that animates, over a faint ghost of the whole run |
//! | a 2D field | a heatmap that animates, on one colour scale throughout |
//! | a 3D field | every z-slice as a montage, on one colour scale, animating together |
//! | points in space | a rotatable 3D scene, depth-sorted, that animates |
//!
//! # The scale is fixed across a run, everywhere, and that is not a detail
//!
//! Every view here holds one scale for the whole run. A picture that renormalises per frame makes
//! a quantity *look* constant while it changes by orders of magnitude, and makes a decay look
//! like a steady state — which is the one thing a picture of a simulation must never do. It is
//! also the easiest thing in the world to do by accident, because per-frame normalisation is what
//! you get if you do not think about it.
//!
//! # Why this is a crate and not a module in an application
//!
//! It was a module in an application, and that application is `publish = false`. So a consumer
//! who could state a simulation and run it could not draw it: the shape of the answer and every
//! view of it were both locked inside a binary nobody can depend on. That is the researcher this
//! layer is for — someone who can write down the physics and does not want to write a plotting
//! stack to see it.
//!
//! # No dependencies
//!
//! SVG and HTML are text, so a `format!` and a file write is the whole renderer: no image
//! encoder, no font, no plotting library, and no network at run time. The HTML report inlines its
//! own viewer for the same reason — the promise is *open the file and it works*, and a script tag
//! pointing at a CDN breaks that on the first machine without a network.

#![deny(missing_docs)]
#![forbid(unsafe_code)]

pub mod data;
pub mod filmstrip;
pub mod report;

pub use data::{readings_csv, to_json};
pub use filmstrip::svg;
pub use report::html;
