// The browser path, exercised without a browser.
//
// `runtime/viewer` learned this once already: a window nobody can photograph proves the program
// did not panic, which is a much weaker claim than it looks. A page nobody clicks is the same
// claim again. This module imports nothing — no wasm-bindgen, no environment — so any host can
// instantiate it, and Node is a host. What runs below is exactly what the page runs: the same
// bytes, the same exports, the same length-prefixed contract.
//
//   cargo build --release -p editor-wasm --target wasm32-unknown-unknown
//   node editor-wasm/selftest.mjs
//
// Exits non-zero on the first failed claim, so it can be a gate step rather than a thing to
// read.

import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';

const here = dirname(fileURLToPath(import.meta.url));
const wasmPath = join(here, '..', 'target', 'wasm32-unknown-unknown', 'release', 'editor_wasm.wasm');

const { instance } = await WebAssembly.instantiate(readFileSync(wasmPath), {});
const w = instance.exports;
const mem = w.memory;
const enc = new TextEncoder(), dec = new TextDecoder();

function take(ptr) {
  const len = new DataView(mem.buffer).getUint32(ptr, true);
  const body = dec.decode(new Uint8Array(mem.buffer, ptr + 4, len));
  w.dualis_free(ptr, len + 4);
  return JSON.parse(body);
}
function call(fn, str) {
  const bytes = enc.encode(str);
  const ptr = w.dualis_alloc(bytes.length);
  new Uint8Array(mem.buffer, ptr, bytes.length).set(bytes);
  const out = take(fn(ptr, bytes.length));
  w.dualis_free(ptr, bytes.length);
  return out;
}

let failures = 0;
function ok(claim, condition, detail = '') {
  if (condition) { console.log(`  ok    ${claim}`); }
  else { console.log(`  FAIL  ${claim}${detail ? ' — ' + detail : ''}`); failures++; }
}

const HOT = JSON.stringify({
  title: 'a block hot enough to glow',
  duration_s: 2.0, frames: 4,
  domains: [{
    kind: 'block', name: 'block', cells: [6, 6, 6], cell_mm: 5.0,
    initial_c: 1200.0, hot_spot: { at: [3, 3, 3], above_k: 400.0 },
  }],
});

console.log('the browser path, in a host that is not a browser');

// A scene checks, and its geometry comes back for the page to wireframe.
const checked = call(w.dualis_check, HOT);
ok('a valid scene checks clean', !checked.error, checked.error);
ok('the summary names the scene', (checked.summary || '').includes('hot enough to glow'));
ok('one placed box comes back', checked.boxes.length === 1);
ok('the box has eight corners', checked.boxes[0].corners.length === 8);
ok('the bounds are finite', Array.isArray(checked.bounds) && checked.bounds.every(Number.isFinite));

// A bad scene is an error with a position, not a crash.
const bad = call(w.dualis_check, '{ "title": ');
ok('a truncated scene reports line:column', /^1:/.test(bad.error || ''), bad.error);

// A scene that did not check must not be runnable off the back of the stored text: the module
// refuses with the check's own message rather than deriving a second, worse one. This ordering
// -- check something bad, then press run -- is what a page with a missing guard would do, and
// it is how this claim came to be here.
const refused = take(w.dualis_run());
ok('run refuses a scene that did not check', /^1:/.test(refused.error || ''), refused.error);

// The run is the CLI's run.
call(w.dualis_check, HOT);
const ran = take(w.dualis_run());
ok('the scene runs', !ran.error, ran.error);
ok('it produced its frames', ran.frames === 5, `got ${ran.frames}`);

// The draw call hands back primitives with every colour resolved.
const drawn = call(w.dualis_draw, JSON.stringify({
  azimuth: 0.7, elevation: 0.4, distance: 2.5, scale: 1.0,
  aspect: 1.6, frame: 0, fit: true,
}));
ok('the wireframe has twelve edges per box', drawn.lines.length === 12);
ok('the field became cells to paint', drawn.dots.length > 100, `${drawn.dots.length} dots`);
ok('every coordinate is finite', drawn.dots.every(d => Number.isFinite(d[0]) && Number.isFinite(d[1])));
ok('the camera came back fitted', Number.isFinite(drawn.camera.scale) && drawn.camera.scale > 0);

// And the colour is Planck's, not a palette: a 1473 K block is orange, so across the cells the
// red channel dominates the blue. This is the claim the whole colour module exists for, checked
// through the boundary the page actually uses.
ok('the canvas says the colour is computed',
   (drawn.notes || []).some(n => n.includes("Planck's")), JSON.stringify(drawn.notes));
const channels = drawn.dots
  .map(d => /rgba\((\d+),(\d+),(\d+)/.exec(d[3]))
  .filter(Boolean)
  .map(m => [+m[1], +m[2], +m[3]]);
ok('the cells carry rgba colours', channels.length > 100);
const reddest = channels.filter(c => c[0] > 8);
ok('a glowing block runs red over blue',
   reddest.length > 0 && reddest.every(c => c[0] >= c[2]),
   `${reddest.length} lit cells`);

// A cool scene must NOT claim a computed colour — the fallback is the honest half.
const COOL = HOT.replace('"initial_c":1200', '"initial_c":20').replace('"initial_c": 1200', '"initial_c": 20');
call(w.dualis_check, COOL);
take(w.dualis_run());
const cool = call(w.dualis_draw, JSON.stringify({ aspect: 1.6, frame: 0, fit: true }));
ok('a cool field says it is false colour',
   (cool.notes || []).some(n => n.includes('false colour')), JSON.stringify(cool.notes));

// The battery runs in the page too, and says the same things.
call(w.dualis_check, HOT);
const verified = take(w.dualis_verify(0));
ok('verify returns a report', typeof verified.report === 'string' && verified.report.length > 50);
ok('the report carries the determinism line', (verified.report || '').includes('determinism'));
ok('a clean scene has no findings', verified.findings === 0, `${verified.findings} findings`);

console.log(failures ? `\n${failures} failed` : '\nall claims held');
process.exit(failures ? 1 : 0);
