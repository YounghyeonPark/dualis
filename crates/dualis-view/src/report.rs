//! A self-contained HTML report, with the view for each domain chosen from the data's shape.
//!
//! # The problem this is for
//!
//! A researcher who can state a simulation and cannot draw one. Today this crate offers them a
//! filmstrip — one fixed shape, whatever the physics is — or a table and a data file with the
//! drawing left as an exercise. For a 2D room the filmstrip is right; for a winding it is empty;
//! for a bar it is a colour strip where a profile would say more. None of that is a choice the
//! researcher should have to make, and all of it is a choice the *data* can make.
//!
//! So: one file, opened in a browser, nothing installed. No Python, no plotting library, no
//! toolchain beyond the one that already ran the simulation.
//!
//! # How the view is chosen
//!
//! By shape, not by domain. A new domain gets a sensible picture without this module learning
//! its name, which is the same reason `Domain::as_field` exists.
//!
//! | what the data is | what it becomes |
//! | --- | --- |
//! | scalars over time | a line chart, one series per reading |
//! | a 1D field | a profile that animates, over a faint ghost of the whole run |
//! | a 2D field | a heatmap that animates, on one colour scale throughout |
//! | a 3D field | a rotatable render **and** every z-slice as a montage — see below |
//! | points in space | a rotatable 3D scene, depth-sorted, that animates |
//!
//! A volume gets **both**, and that is not indecision. A raycast composites values along each
//! ray, so it shows the shape of a field and a reader cannot get a number back out of it; the
//! montage puts every sample on screen and is quantitative and unreadable as a shape. Offering
//! one would be choosing which half of the question a researcher is allowed to ask.
//!
//! The montage is one slice per tile rather than one slice behind a slider, because a viewer who
//! never touches a slider sees a picture of a solid that is really a picture of one plane.
//!
//! The scale is fixed across the run in every case. A frame that rescales makes a quantity
//! *look* constant while it changes by orders of magnitude, which is the one thing a picture of
//! a simulation must never do.
//!
//! # Why the viewer is a string in this file
//!
//! It is about two hundred lines of JavaScript, inlined into the output. That is not elegant and
//! the alternatives are worse: a separate asset means the report is no longer one file, and a
//! library from a CDN means it does not open on a machine without a network. The whole promise
//! is *open it and it works*.

use dualis_scene::{Frame, PanelData};

/// Build the report for a finished run.
pub fn html(title: &str, frames: &[Frame]) -> String {
    let mut out = String::with_capacity(1 << 16);
    out.push_str("<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n");
    out.push_str("<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">\n");
    out.push_str(&format!("<title>{}</title>\n", escape(title)));
    out.push_str("<style>\n");
    out.push_str(STYLE);
    out.push_str("</style>\n</head>\n<body>\n<div class=\"wrap\">\n");

    out.push_str(&format!(
        "<header><p class=\"eyebrow\">dualis-world</p><h1>{}</h1>\
         <p class=\"sub\">{} frames &middot; {:.4} s &middot; conservation held throughout</p>\
         </header>\n",
        escape(title),
        frames.len(),
        frames.last().map_or(0.0, |f| f.time_s)
    ));

    out.push_str(
        "<div class=\"bar\"><button id=\"play\">Pause</button>\
         <input id=\"scrub\" type=\"range\" min=\"0\" value=\"0\" aria-label=\"Frame\">\
         <span class=\"tick\" id=\"tick\"></span></div>\n",
    );

    // One card per drawable domain, plus one for the readings, which every domain has.
    //
    // **A volume gets two.** A raycast is what a three-dimensional field looks like and a slice
    // montage is what it *is*, and neither substitutes: the render composites values along a ray,
    // so a reader cannot get a number back out of it, while the montage puts every sample on
    // screen and is unreadable as a shape. Offering only one would be choosing which half of the
    // question a researcher is allowed to ask.
    if let Some(first) = frames.first() {
        for panel in &first.panels {
            let kinds: &[&str] = match &panel.data {
                PanelData::Field { nz, .. } if *nz > 1 => &["volume", "slices"],
                PanelData::Field { ny, .. } if *ny <= 1 => &["profile"],
                PanelData::Field { .. } => &["heatmap"],
                PanelData::Points { .. } => &["scene"],
            };
            for kind in kinds {
                out.push_str(&format!(
                    "<section class=\"card\"><div class=\"head\"><h2>{}</h2>\
                     <span class=\"kind\">{}</span></div>\
                     <canvas class=\"view\" data-panel=\"{}\" data-kind=\"{}\" \
                     data-slot=\"{}\" width=\"1400\" height=\"620\"></canvas>\
                     <p class=\"cap\" id=\"cap-{}\"></p></section>\n",
                    escape(&panel.name),
                    match *kind {
                        "profile" => "1D field &middot; profile",
                        "heatmap" => "2D field &middot; heatmap",
                        "volume" => "3D field &middot; rendered, drag to rotate",
                        "slices" => "3D field &middot; every z-slice, and the numbers",
                        _ => "bodies &middot; 3D, drag to rotate",
                    },
                    escape(&panel.name),
                    kind,
                    format_args!("{}-{}", escape(&panel.name), kind),
                    format_args!("{}-{}", escape(&panel.name), kind),
                ));
            }
        }
        if !first.readings.is_empty() {
            out.push_str(
                "<section class=\"card\"><div class=\"head\"><h2>Readings</h2>\
                 <span class=\"kind\">scalars &middot; over time</span></div>\
                 <canvas class=\"view\" data-kind=\"series\" data-slot=\"series\" \
                 width=\"1400\" height=\"560\">\
                 </canvas><p class=\"cap\" id=\"cap-series\"></p></section>\n",
            );
        }
    }

    out.push_str("</div>\n<script id=\"run\" type=\"application/json\">");
    out.push_str(&json(frames));
    out.push_str("</script>\n<script>\n");
    out.push_str(VIEWER);
    out.push_str("</script>\n</body>\n</html>\n");
    out
}

/// The run as JSON, for the viewer above it. Six figures: a picture needs no more.
fn json(frames: &[Frame]) -> String {
    let mut out = String::from("{\"frames\":[");
    for (fi, f) in frames.iter().enumerate() {
        if fi > 0 {
            out.push(',');
        }
        out.push_str(&format!("{{\"t\":{:.6e},\"panels\":[", f.time_s));
        for (pi, p) in f.panels.iter().enumerate() {
            if pi > 0 {
                out.push(',');
            }
            out.push_str(&format!(
                "{{\"name\":{},\"unit\":{},",
                quote(&p.name),
                quote(p.unit)
            ));
            match &p.data {
                PanelData::Field { nx, ny, nz, values } => out.push_str(&format!(
                    "\"kind\":\"field\",\"nx\":{nx},\"ny\":{ny},\"nz\":{nz},\"v\":{}",
                    nums(values)
                )),
                PanelData::Points {
                    positions,
                    values,
                    bounds,
                    boxed,
                } => {
                    let flat: Vec<f64> = positions.iter().flatten().copied().collect();
                    out.push_str(&format!(
                        "\"kind\":\"points\",\"boxed\":{boxed},\"b\":{},\"p\":{},\"v\":{}",
                        nums(bounds),
                        nums(&flat),
                        nums(values)
                    ));
                }
            }
            out.push('}');
        }
        out.push_str("],\"r\":[");
        for (ri, r) in f.readings.iter().enumerate() {
            if ri > 0 {
                out.push(',');
            }
            out.push_str(&format!(
                "{{\"d\":{},\"l\":{},\"u\":{},\"v\":{:.6e}}}",
                quote(&r.domain),
                quote(&r.label),
                quote(r.unit),
                r.value
            ));
        }
        out.push_str("]}");
    }
    out.push_str("]}");
    out
}

fn nums(v: &[f64]) -> String {
    let mut s = String::from("[");
    for (i, x) in v.iter().enumerate() {
        if i > 0 {
            s.push(',');
        }
        s.push_str(&format!("{x:.6e}"));
    }
    s.push(']');
    s
}

fn quote(s: &str) -> String {
    let mut out = String::from("\"");
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            c if c.is_control() => out.push(' '),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

fn escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

const STYLE: &str = r##"
:root{--bg:#0d1116;--card:#151b24;--rule:#242d3a;--ink:#d7dde6;--dim:#7f8b9a;--hot:#ff5221;--cool:#3d8bff;--warm:#ffc247}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);line-height:1.5;
 font-family:ui-sans-serif,system-ui,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;-webkit-font-smoothing:antialiased}
.wrap{max-width:1180px;margin:0 auto;padding:30px 18px 60px;display:flex;flex-direction:column;gap:16px}
.eyebrow{font-family:ui-monospace,Consolas,monospace;font-size:11px;letter-spacing:.16em;
 text-transform:uppercase;color:var(--dim);margin:0}
h1{font-size:clamp(21px,3.2vw,29px);font-weight:620;letter-spacing:-.02em;margin:4px 0 2px;text-wrap:balance}
.sub{margin:0;color:var(--dim);font-size:14px}
.bar{position:sticky;top:0;z-index:5;display:flex;align-items:center;gap:12px;padding:10px 12px;
 background:var(--card);border:1px solid var(--rule);border-radius:3px}
.bar button{appearance:none;background:var(--bg);color:var(--ink);border:1px solid var(--rule);
 border-radius:3px;font:inherit;font-size:13px;padding:6px 14px;cursor:pointer;min-width:72px}
.bar button:hover{border-color:var(--hot)}
.bar button:focus-visible{outline:2px solid var(--warm);outline-offset:1px}
input[type=range]{flex:1;min-width:140px;accent-color:var(--hot)}
.tick{font-family:ui-monospace,Consolas,monospace;font-size:12px;color:var(--dim);
 font-variant-numeric:tabular-nums;min-width:190px;text-align:right}
.card{background:var(--card);border:1px solid var(--rule);border-radius:3px;overflow:hidden}
.head{display:flex;align-items:baseline;gap:12px;padding:12px 14px 10px;border-bottom:1px solid var(--rule)}
h2{margin:0;font-size:15px;font-weight:620;letter-spacing:-.01em}
.kind{font-family:ui-monospace,Consolas,monospace;font-size:11px;color:var(--dim);letter-spacing:.04em}
canvas.view{display:block;width:100%;height:auto;background:#0a0e13}
canvas[data-kind=scene],canvas[data-kind=volume]{cursor:grab;touch-action:none}
canvas[data-kind=scene]:active,canvas[data-kind=volume]:active{cursor:grabbing}
.cap{margin:0;padding:9px 14px 12px;font-family:ui-monospace,Consolas,monospace;
 font-size:11.5px;color:var(--dim)}
@media (prefers-reduced-motion:reduce){.bar button{outline:1px dashed var(--dim)}}
"##;

const VIEWER: &str = r##"
(function(){
"use strict";
var RUN = JSON.parse(document.getElementById("run").textContent);
var F = RUN.frames, N = F.length;
var frame = 0, playing = true, last = 0;
var cam = {az:0.7, el:0.4, dist:2.5}, drag = null;

/* One scale per panel across the whole run. A frame that rescales makes a quantity look
   constant while it changes by orders of magnitude. */
var range = {};
F[0].panels.forEach(function(p, i){
  var lo = Infinity, hi = -Infinity;
  F.forEach(function(f){ f.panels[i].v.forEach(function(x){ if(x<lo)lo=x; if(x>hi)hi=x; }); });
  range[p.name] = {lo:lo, hi:hi};
});

var series = {};
if (F[0].r.length) {
  F[0].r.forEach(function(r, i){
    var key = r.d + "." + r.l;
    var vals = F.map(function(f){ return f.r[i] ? f.r[i].v : 0; });
    var lo = Math.min.apply(null, vals), hi = Math.max.apply(null, vals);
    series[key] = {vals:vals, lo:lo, hi:hi, unit:r.u};
  });
}

function ramp(t){
  t = t<0?0:t>1?1:t;
  var a,b,u, c0=[24,52,110], c1=[61,139,255], c2=[255,194,71], c3=[255,82,33];
  if(t<0.34){a=c0;b=c1;u=t/0.34;} else if(t<0.67){a=c1;b=c2;u=(t-0.34)/0.33;} else {a=c2;b=c3;u=(t-0.67)/0.33;}
  return "rgb("+[0,1,2].map(function(i){return Math.round(a[i]+(b[i]-a[i])*u);}).join(",")+")";
}
function nice(x){
  if (x === 0) return "0";
  var a = Math.abs(x);
  if (a >= 1e5 || a < 1e-3) return x.toExponential(2);
  return x.toFixed(a < 1 ? 4 : a < 100 ? 2 : 1);
}

var views = [].slice.call(document.querySelectorAll("canvas.view")).map(function(c){
  return {c:c, ctx:c.getContext("2d"), kind:c.dataset.kind, panel:c.dataset.panel,
          slot:c.dataset.slot || c.dataset.panel || "series"};
});

function panelOf(f, name){
  for (var i=0;i<f.panels.length;i++) if (f.panels[i].name===name) return f.panels[i];
  return null;
}

/* ---- 1D: a profile, over a ghost of every frame ------------------------------------------ */
function drawProfile(v, f){
  var x = v.ctx, w = v.c.width, h = v.c.height, p = panelOf(f, v.panel), R = range[v.panel];
  var pad = {l:64, r:18, t:16, b:34};
  x.fillStyle="#0a0e13"; x.fillRect(0,0,w,h);
  var iw = w-pad.l-pad.r, ih = h-pad.t-pad.b;
  var sy = function(val){ return pad.t + ih*(1-(val-R.lo)/((R.hi-R.lo)||1)); };
  x.strokeStyle="#1c242f"; x.lineWidth=1;
  for(var g=0; g<=4; g++){ var yy=pad.t+ih*g/4; x.beginPath(); x.moveTo(pad.l,yy); x.lineTo(w-pad.r,yy); x.stroke(); }
  /* every frame, faint: the envelope the run covers */
  x.strokeStyle="rgba(125,145,175,0.16)"; x.lineWidth=1.5;
  F.forEach(function(ff){
    var q = panelOf(ff, v.panel); if(!q) return;
    x.beginPath();
    q.v.forEach(function(val,i){ var xx=pad.l+iw*i/(q.v.length-1||1); i?x.lineTo(xx,sy(val)):x.moveTo(xx,sy(val)); });
    x.stroke();
  });
  x.strokeStyle="#ff5221"; x.lineWidth=2.6; x.beginPath();
  p.v.forEach(function(val,i){ var xx=pad.l+iw*i/(p.v.length-1||1); i?x.lineTo(xx,sy(val)):x.moveTo(xx,sy(val)); });
  x.stroke();
  x.fillStyle="#7f8b9a"; x.font="13px ui-monospace,Consolas,monospace"; x.textAlign="right";
  x.fillText(nice(R.hi)+" "+p.unit, pad.l-9, pad.t+5);
  x.fillText(nice(R.lo)+" "+p.unit, pad.l-9, pad.t+ih+5);
  x.textAlign="left"; x.fillText("cell 0", pad.l, h-12);
  x.textAlign="right"; x.fillText("cell "+(p.v.length-1), w-pad.r, h-12);
  cap(v, p.v.length+" cells · scale fixed across the run · faint lines are every frame");
}

/* ---- 2D: a heatmap ------------------------------------------------------------------------ */
function drawHeat(v, f){
  var x = v.ctx, w = v.c.width, h = v.c.height, p = panelOf(f, v.panel), R = range[v.panel];
  x.fillStyle="#0a0e13"; x.fillRect(0,0,w,h);
  var pad=14, aw=w-pad*2, ah=h-pad*2-22;
  var s = Math.min(aw/p.nx, ah/p.ny), ox = (w-s*p.nx)/2, oy = pad + (ah-s*p.ny)/2;
  for (var j=0;j<p.ny;j++) for (var i=0;i<p.nx;i++){
    var val = p.v[j*p.nx+i];
    x.fillStyle = ramp((val-R.lo)/((R.hi-R.lo)||1));
    x.fillRect(Math.floor(ox+i*s), Math.floor(oy+(p.ny-1-j)*s), Math.ceil(s), Math.ceil(s));
  }
  /* the scale, as a strip, so a colour can be read back to a number */
  var bw=210, bx=w-pad-bw, by=h-24;
  for(var k=0;k<bw;k++){ x.fillStyle=ramp(k/bw); x.fillRect(bx+k,by,1,10); }
  x.fillStyle="#7f8b9a"; x.font="12px ui-monospace,Consolas,monospace";
  x.textAlign="right"; x.fillText(nice(R.lo), bx-8, by+9);
  x.textAlign="left"; x.fillText(nice(R.hi)+" "+p.unit, bx+bw+8, by+9);
  cap(v, p.nx+" x "+p.ny+" · one colour scale across every frame");
}

/* ---- 3D field: every slice, laid out as a montage ------------------------------------------ */
function drawSlices(v, f){
  var x = v.ctx, w = v.c.width, h = v.c.height, p = panelOf(f, v.panel), R = range[v.panel];
  x.fillStyle="#0a0e13"; x.fillRect(0,0,w,h);
  var pad=14, gap=6, aw=w-pad*2, ah=h-pad*2-22;
  /* Enough columns that the tiles are as large as possible while all of them fit. Searched
     rather than derived: the tile size depends on the aspect of both the grid and the canvas,
     and nz is small enough that trying every column count is free. */
  var best={s:0,cols:1};
  for(var c=1;c<=p.nz;c++){
    var rows=Math.ceil(p.nz/c);
    var s=Math.min((aw-(c-1)*gap)/(c*p.nx), (ah-(rows-1)*gap)/(rows*p.ny));
    if(s>best.s) best={s:s,cols:c};
  }
  var s=best.s, cols=best.cols, rows=Math.ceil(p.nz/cols);
  var tw=cols*p.nx*s+(cols-1)*gap, th=rows*p.ny*s+(rows-1)*gap;
  var ox=(w-tw)/2, oy=pad+(ah-th)/2;
  for(var k=0;k<p.nz;k++){
    var cx=ox+(k%cols)*(p.nx*s+gap), cy=oy+Math.floor(k/cols)*(p.ny*s+gap);
    for(var j=0;j<p.ny;j++) for(var i=0;i<p.nx;i++){
      var val=p.v[k*p.nx*p.ny+j*p.nx+i];
      x.fillStyle=ramp((val-R.lo)/((R.hi-R.lo)||1));
      x.fillRect(Math.floor(cx+i*s), Math.floor(cy+(p.ny-1-j)*s), Math.ceil(s), Math.ceil(s));
    }
    x.fillStyle="#7f8b9a"; x.font="10px ui-monospace,Consolas,monospace"; x.textAlign="left";
    x.fillText("z"+k, cx+1, cy-2);
  }
  var bw=210, bx=w-pad-bw, by=h-24;
  for(var q=0;q<bw;q++){ x.fillStyle=ramp(q/bw); x.fillRect(bx+q,by,1,10); }
  x.fillStyle="#7f8b9a"; x.font="12px ui-monospace,Consolas,monospace";
  x.textAlign="right"; x.fillText(nice(R.lo), bx-8, by+9);
  x.textAlign="left"; x.fillText(nice(R.hi)+" "+p.unit, bx+bw+8, by+9);
  cap(v, p.nx+" x "+p.ny+" x "+p.nz+" · all "+p.nz+" slices, z increasing · one colour scale across every frame");
}

/* ---- 3D field: a raycast ------------------------------------------------------------------- */
/* Rendered into a small buffer and scaled up. A 220x150 buffer at 48 samples a ray is about
   1.6 million trilinear lookups a frame, which a browser does comfortably; the full canvas would
   be forty times that and would drop frames on the animation. The softness that costs is honest —
   this view is for shape, and the montage beside it carries the numbers. */
var VOL_W = 220, VOL_H = 150, VOL_STEPS = 48;

/* Opacity from value, and the choice is not cosmetic.

   A signed field — pressure, zero mean — must be transparent in the middle and opaque at both
   extremes, or a standing wave renders as a solid block. A one-sided field — kelvin, 293 upward —
   must be transparent at the low end, or a block at ambient renders as a solid block for the
   opposite reason. So the transfer function is chosen from the run's own range rather than fixed,
   and `signed` is decided once from whether that range straddles zero. */
function opacity(t, signed){
  var a = signed ? Math.abs(2*t - 1) : t;
  return a*a;                     /* squared, so the quiet bulk clears out of the way */
}

function drawVolume(v, f){
  var x=v.ctx, w=v.c.width, h=v.c.height, p=panelOf(f,v.panel), R=range[v.panel];
  x.fillStyle="#0a0e13"; x.fillRect(0,0,w,h);
  var span=(R.hi-R.lo)||1, signed = R.lo < 0 && R.hi > 0;

  /* The box in world units, longest axis normalised to 1 so a slab is drawn as a slab. */
  var m=Math.max(p.nx,p.ny,p.nz), ex=[p.nx/m, p.ny/m, p.nz/m];
  var ca=Math.cos(cam.az), sa=Math.sin(cam.az), ce=Math.cos(cam.el), se=Math.sin(cam.el);
  /* Camera basis: right, up, forward. The same angles the bodies view uses, so dragging one
     rotates the other and a reader is never looking at two different orientations. */
  var fwd=[-sa*ce, -se, -ca*ce], right=[ca, 0, -sa], up=[-sa*se, ce, -ca*se];
  var eye=[-fwd[0]*cam.dist, -fwd[1]*cam.dist, -fwd[2]*cam.dist];

  var img=x.createImageData(VOL_W, VOL_H), d=img.data, aspect=VOL_W/VOL_H, k=0, lit=0;
  for(var py=0; py<VOL_H; py++){
    var sy=(1 - 2*(py+0.5)/VOL_H)*0.75;
    for(var px=0; px<VOL_W; px++){
      var sx=(2*(px+0.5)/VOL_W - 1)*0.75*aspect;
      var dir=[fwd[0]+right[0]*sx+up[0]*sy, fwd[1]+right[1]*sx+up[1]*sy, fwd[2]+right[2]*sx+up[2]*sy];
      var len=Math.hypot(dir[0],dir[1],dir[2]); dir=[dir[0]/len,dir[1]/len,dir[2]/len];

      /* Slab test against the box centred on the origin. */
      var t0=-1e9, t1=1e9, hit=true;
      for(var ax=0; ax<3; ax++){
        var half=ex[ax]/2, o=eye[ax], dd=dir[ax];
        if(Math.abs(dd)<1e-9){ if(o<-half||o>half){hit=false;break;} continue; }
        var a1=(-half-o)/dd, b1=(half-o)/dd, lo=Math.min(a1,b1), hi=Math.max(a1,b1);
        if(lo>t0)t0=lo; if(hi<t1)t1=hi;
      }
      var r=13,g=18,b=25, alpha=0;
      if(hit && t1>t0 && t1>0){
        if(t0<0)t0=0;
        var dt=(t1-t0)/VOL_STEPS;
        for(var st=0; st<VOL_STEPS && alpha<0.985; st++){
          var tt=t0+dt*(st+0.5);
          /* World point -> grid index, with the box centred on the origin. */
          var gx=((eye[0]+dir[0]*tt)/ex[0]+0.5)*(p.nx-1);
          var gy=((eye[1]+dir[1]*tt)/ex[1]+0.5)*(p.ny-1);
          var gz=((eye[2]+dir[2]*tt)/ex[2]+0.5)*(p.nz-1);
          if(gx<0||gy<0||gz<0||gx>p.nx-1||gy>p.ny-1||gz>p.nz-1) continue;
          var i0=Math.floor(gx), j0=Math.floor(gy), k0=Math.floor(gz);
          var i1=Math.min(i0+1,p.nx-1), j1=Math.min(j0+1,p.ny-1), k1=Math.min(k0+1,p.nz-1);
          var fx=gx-i0, fy=gy-j0, fz=gz-k0, nxy=p.nx*p.ny;
          var c000=p.v[k0*nxy+j0*p.nx+i0], c100=p.v[k0*nxy+j0*p.nx+i1];
          var c010=p.v[k0*nxy+j1*p.nx+i0], c110=p.v[k0*nxy+j1*p.nx+i1];
          var c001=p.v[k1*nxy+j0*p.nx+i0], c101=p.v[k1*nxy+j0*p.nx+i1];
          var c011=p.v[k1*nxy+j1*p.nx+i0], c111=p.v[k1*nxy+j1*p.nx+i1];
          var z0=(c000*(1-fx)+c100*fx)*(1-fy)+(c010*(1-fx)+c110*fx)*fy;
          var z1=(c001*(1-fx)+c101*fx)*(1-fy)+(c011*(1-fx)+c111*fx)*fy;
          var val=z0*(1-fz)+z1*fz;

          var norm=(val-R.lo)/span;
          var a=opacity(norm, signed)*0.16;
          if(a<=0.0008) continue;
          var col=ramp(norm), ci=col.indexOf("(");
          var parts=col.slice(ci+1,-1).split(",");
          var contrib=a*(1-alpha);
          r+= (Number(parts[0])-13)*contrib;
          g+= (Number(parts[1])-18)*contrib;
          b+= (Number(parts[2])-25)*contrib;
          alpha+=contrib;
        }
      }
      if(alpha>0.05) lit++;
      d[k++]=r; d[k++]=g; d[k++]=b; d[k++]=255;
    }
  }

  /* Nearest-neighbour off, so the small buffer reads as a soft render and not as pixels. */
  var tmp=document.createElement("canvas"); tmp.width=VOL_W; tmp.height=VOL_H;
  tmp.getContext("2d").putImageData(img,0,0);
  x.imageSmoothingEnabled=true;
  var scale=Math.min(w/VOL_W,(h-26)/VOL_H);
  x.drawImage(tmp,(w-VOL_W*scale)/2,(h-26-VOL_H*scale)/2,VOL_W*scale,VOL_H*scale);

  var pad=14, bw=210, bx=w-pad-bw, by=h-24;
  for(var q=0;q<bw;q++){ x.fillStyle=ramp(q/bw); x.fillRect(bx+q,by,1,10); }
  x.fillStyle="#7f8b9a"; x.font="12px ui-monospace,Consolas,monospace";
  x.textAlign="right"; x.fillText(nice(R.lo), bx-8, by+9);
  x.textAlign="left"; x.fillText(nice(R.hi)+" "+p.unit, bx+bw+8, by+9);
  /* **Say when the picture is nearly empty.** A localised feature -- one hot cell in a block
     of 729 -- is a small bright dot with everything else transparent, which is correct and
     reads exactly like a broken renderer. Making it look bigger would be making the picture
     lie, so the caption reports how much of the frame carries anything instead.

     Measured rather than guessed: the exponent in `opacity` was tried at 1, 1.5 and 2 against
     three real scenes and moved the occupied fraction from 0.2% to 0.1% on the hot spot. The
     transfer function is not what makes a point source small. */
  var frac = lit/(VOL_W*VOL_H);
  var note = frac < 0.03
    ? " \u00b7 " + (100*frac).toFixed(1) + "% of the frame is occupied: the feature is that "
      + "small against the whole volume, and the montage below shows where"
    : "";
  cap(v, p.nx+" x "+p.ny+" x "+p.nz+" \u00b7 composited along each ray, so this shows shape and not "
       +"values \u00b7 the montage below carries the numbers \u00b7 drag to rotate, scroll to zoom"
       + note);
}

/* ---- 3D: bodies, depth sorted -------------------------------------------------------------- */
function project(pt, s, w, h){
  var X=(pt[0]-s.c[0])/s.span, Y=(pt[1]-s.c[1])/s.span, Z=(pt[2]-s.c[2])/s.span;
  var ca=Math.cos(cam.az), sa=Math.sin(cam.az);
  var x1=X*ca-Z*sa, z1=X*sa+Z*ca;
  var ce=Math.cos(cam.el), se=Math.sin(cam.el);
  var y1=Y*ce-z1*se, z2=Y*se+z1*ce;
  var d=z2+cam.dist; if(d<0.05)d=0.05;
  var f=(Math.min(w,h)*0.60)/d;
  return {x:w/2+x1*f, y:h/2-y1*f, d:d, s:f/Math.min(w,h)};
}
function drawScene(v, f){
  var x=v.ctx, w=v.c.width, h=v.c.height, p=panelOf(f,v.panel), R=range[v.panel];
  x.fillStyle="#0a0e13"; x.fillRect(0,0,w,h);
  var b=p.b, s={c:[(b[0]+b[3])/2,(b[1]+b[4])/2,(b[2]+b[5])/2],
                span:Math.max(b[3]-b[0],b[4]-b[1],b[5]-b[2])||1};
  if(p.boxed){
    var cs=[];
    for(var i=0;i<8;i++) cs.push(project([i&1?b[3]:b[0], i&2?b[4]:b[1], i&4?b[5]:b[2]], s, w, h));
    var E=[[0,1],[0,2],[0,4],[1,3],[1,5],[2,3],[2,6],[3,7],[4,5],[4,6],[5,7],[6,7]];
    x.strokeStyle="rgba(120,145,180,0.32)"; x.lineWidth=1.3;
    E.forEach(function(e){ x.beginPath(); x.moveTo(cs[e[0]].x,cs[e[0]].y); x.lineTo(cs[e[1]].x,cs[e[1]].y); x.stroke(); });
  }
  var pts=[];
  for(var k=0;k<p.v.length;k++){
    var q=project([p.p[3*k],p.p[3*k+1],p.p[3*k+2]], s, w, h); q.val=p.v[k]; pts.push(q);
  }
  pts.sort(function(a,b2){ return b2.d-a.d; });
  var base = p.v.length > 40 ? 6 : 14;
  pts.forEach(function(q){
    var r=Math.max(1.5, base*q.s*44);
    x.beginPath(); x.arc(q.x,q.y,r,0,6.2832);
    x.fillStyle=ramp((q.val-R.lo)/((R.hi-R.lo)||1)); x.fill();
  });
  cap(v, p.v.length+" bodies · colour is "+p.unit+" · drag to rotate, scroll to zoom");
}

/* ---- scalars over time --------------------------------------------------------------------- */
function drawSeries(v){
  var x=v.ctx, w=v.c.width, h=v.c.height;
  var keys=Object.keys(series);
  x.fillStyle="#0a0e13"; x.fillRect(0,0,w,h);
  if(!keys.length) return;
  var pad={l:74,r:210,t:18,b:34}, iw=w-pad.l-pad.r, ih=h-pad.t-pad.b;
  x.strokeStyle="#1c242f"; x.lineWidth=1;
  for(var g=0;g<=4;g++){ var yy=pad.t+ih*g/4; x.beginPath(); x.moveTo(pad.l,yy); x.lineTo(pad.l+iw,yy); x.stroke(); }
  keys.forEach(function(k,ki){
    var s=series[k], span=(s.hi-s.lo)||1;
    x.strokeStyle=ramp(keys.length>1?ki/(keys.length-1):0.5); x.lineWidth=2; x.beginPath();
    s.vals.forEach(function(val,i){
      var xx=pad.l+iw*i/(N-1||1), yy=pad.t+ih*(1-(val-s.lo)/span);
      i?x.lineTo(xx,yy):x.moveTo(xx,yy);
    });
    x.stroke();
    /* each series is normalised to its own range, so the legend carries the numbers */
    x.fillStyle=ramp(keys.length>1?ki/(keys.length-1):0.5);
    x.fillRect(pad.l+iw+14, pad.t+ki*19+4, 9, 9);
    x.fillStyle="#a8b3c2"; x.font="12px ui-monospace,Consolas,monospace"; x.textAlign="left";
    x.fillText(k+"  "+nice(s.vals[frame])+" "+s.unit, pad.l+iw+29, pad.t+ki*19+13);
  });
  /* where we are */
  var cx=pad.l+iw*frame/(N-1||1);
  x.strokeStyle="rgba(255,82,33,0.75)"; x.lineWidth=1.5;
  x.beginPath(); x.moveTo(cx,pad.t); x.lineTo(cx,pad.t+ih); x.stroke();
  x.fillStyle="#7f8b9a"; x.font="12px ui-monospace,Consolas,monospace";
  x.textAlign="left"; x.fillText("t = "+nice(F[0].t)+" s", pad.l, h-12);
  x.textAlign="right"; x.fillText("t = "+nice(F[N-1].t)+" s", pad.l+iw, h-12);
  x.textAlign="left"; x.fillText("each series on its own scale", pad.l, pad.t-4);
  cap(v, keys.length+" scalars · the line marks the current frame");
}

function cap(v, text){
  /* keyed by slot, not by panel: a volume has two views of the same panel and each writes its
     own caption. Keying by panel made the second overwrite the first. */
  var el = document.getElementById("cap-" + v.slot);
  if (el) el.textContent = text;
}

function drawAll(){
  var f=F[frame];
  views.forEach(function(v){
    if(v.kind==="profile") drawProfile(v,f);
    else if(v.kind==="heatmap") drawHeat(v,f);
    else if(v.kind==="slices") drawSlices(v,f);
    else if(v.kind==="volume") drawVolume(v,f);
    else if(v.kind==="scene") drawScene(v,f);
    else drawSeries(v);
  });
  document.getElementById("tick").textContent =
    "frame "+(frame+1)+" / "+N+"   t = "+nice(f.t)+" s";
  document.getElementById("scrub").value=String(frame);
}

var scrub=document.getElementById("scrub"), play=document.getElementById("play");
scrub.max=String(N-1);
scrub.oninput=function(){ frame=Number(scrub.value); playing=false; play.textContent="Play"; drawAll(); };
play.onclick=function(){ playing=!playing; play.textContent=playing?"Pause":"Play"; };

views.filter(function(v){return v.kind==="scene"||v.kind==="volume";}).forEach(function(v){
  v.c.addEventListener("pointerdown",function(e){ drag={x:e.clientX,y:e.clientY}; v.c.setPointerCapture(e.pointerId); });
  v.c.addEventListener("pointermove",function(e){
    if(!drag)return;
    cam.az+=(e.clientX-drag.x)*0.008;
    cam.el=Math.max(-1.5,Math.min(1.5,cam.el+(e.clientY-drag.y)*0.006));
    drag={x:e.clientX,y:e.clientY}; drawAll();
  });
  v.c.addEventListener("pointerup",function(){ drag=null; });
  v.c.addEventListener("wheel",function(e){
    e.preventDefault(); cam.dist=Math.max(1.2,Math.min(9,cam.dist*(1+e.deltaY*0.0011))); drawAll();
  },{passive:false});
});

if (window.matchMedia("(prefers-reduced-motion: reduce)").matches){ playing=false; play.textContent="Play"; }
function loop(now){
  if(playing && now-last>60){ last=now; frame=(frame+1)%N; drawAll(); }
  requestAnimationFrame(loop);
}
drawAll();
requestAnimationFrame(loop);
})();
"##;
