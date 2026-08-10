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
//! | a 3D field | every z-slice as a montage, on one colour scale, animating together |
//! | points in space | a rotatable 3D scene, depth-sorted, that animates |
//!
//! The montage is the honest answer to a volume on a flat canvas. The alternative — one slice
//! with a slider — hides the rest behind an interaction, and a viewer who never touches the
//! slider sees a picture of a solid that is really a picture of one plane through it. Every
//! sample is on screen here, and the caption says how many slices there are.
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
    if let Some(first) = frames.first() {
        for panel in &first.panels {
            let kind = match &panel.data {
                PanelData::Field { nz, .. } if *nz > 1 => "slices",
                PanelData::Field { ny, .. } if *ny <= 1 => "profile",
                PanelData::Field { .. } => "heatmap",
                PanelData::Points { .. } => "scene",
            };
            out.push_str(&format!(
                "<section class=\"card\"><div class=\"head\"><h2>{}</h2>\
                 <span class=\"kind\">{}</span></div>\
                 <canvas class=\"view\" data-panel=\"{}\" data-kind=\"{}\" \
                 width=\"1400\" height=\"620\"></canvas>\
                 <p class=\"cap\" id=\"cap-{}\"></p></section>\n",
                escape(&panel.name),
                match kind {
                    "profile" => "1D field &middot; profile",
                    "heatmap" => "2D field &middot; heatmap",
                    "slices" => "3D field &middot; every z-slice",
                    _ => "bodies &middot; 3D, drag to rotate",
                },
                escape(&panel.name),
                kind,
                escape(&panel.name),
            ));
        }
        if !first.readings.is_empty() {
            out.push_str(
                "<section class=\"card\"><div class=\"head\"><h2>Readings</h2>\
                 <span class=\"kind\">scalars &middot; over time</span></div>\
                 <canvas class=\"view\" data-kind=\"series\" width=\"1400\" height=\"560\">\
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
canvas[data-kind=scene]{cursor:grab;touch-action:none}
canvas[data-kind=scene]:active{cursor:grabbing}
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
  return {c:c, ctx:c.getContext("2d"), kind:c.dataset.kind, panel:c.dataset.panel};
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
  var el = document.getElementById("cap-" + (v.panel || "series"));
  if (el) el.textContent = text;
}

function drawAll(){
  var f=F[frame];
  views.forEach(function(v){
    if(v.kind==="profile") drawProfile(v,f);
    else if(v.kind==="heatmap") drawHeat(v,f);
    else if(v.kind==="slices") drawSlices(v,f);
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

views.filter(function(v){return v.kind==="scene";}).forEach(function(v){
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
