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


// --- CAD, which is the whole reason this page is an IDE and not a viewer ----------------------
// A scene's `parts` names a file. On a machine that is a path; here there is no filesystem, so
// the page hands the bytes over and the name is a label. The claims below run that path end to
// end -- a real binary STL, built here from the format's own layout, through the boundary, into
// the voxeliser -- in a host that is not a browser.

function putBytes(u8) {
  const ptr = w.dualis_alloc(u8.length);
  new Uint8Array(mem.buffer, ptr, u8.length).set(u8);
  return [ptr, u8.length];
}
function sendPart(name, u8) {
  const nb = enc.encode(name);
  const np = w.dualis_alloc(nb.length);
  new Uint8Array(mem.buffer, np, nb.length).set(nb);
  const [dp, dl] = putBytes(u8);
  const out = take(w.dualis_part(np, nb.length, dp, dl));
  w.dualis_free(np, nb.length);
  w.dualis_free(dp, dl);
  return out;
}
function forgetPart(name) {
  const nb = enc.encode(name);
  const np = w.dualis_alloc(Math.max(1, nb.length));
  if (nb.length) new Uint8Array(mem.buffer, np, nb.length).set(nb);
  const out = take(w.dualis_forget_part(np, nb.length));
  w.dualis_free(np, Math.max(1, nb.length));
  return out;
}

// A binary STL of an axis-aligned box, written from the format: 80 bytes of header, a
// little-endian triangle count, then 50 bytes per triangle. Millimetres, which is what
// `Mesh::from_stl` reads -- a fixture in metres voxelises to nothing, which is a mistake this
// repository has already made once.
function boxStl(lo, hi) {
  const c = [
    [lo[0], lo[1], lo[2]], [hi[0], lo[1], lo[2]], [hi[0], hi[1], lo[2]], [lo[0], hi[1], lo[2]],
    [lo[0], lo[1], hi[2]], [hi[0], lo[1], hi[2]], [hi[0], hi[1], hi[2]], [lo[0], hi[1], hi[2]],
  ];
  const quads = [[0,1,2,3],[5,4,7,6],[4,5,1,0],[3,2,6,7],[4,0,3,7],[1,5,6,2]];
  const tris = [];
  for (const [a, b, cc, d] of quads) { tris.push([c[a], c[b], c[cc]]); tris.push([c[a], c[cc], c[d]]); }
  const buf = new ArrayBuffer(84 + tris.length * 50);
  const dv = new DataView(buf);
  dv.setUint32(80, tris.length, true);
  tris.forEach((t, i) => {
    let o = 84 + i * 50 + 12;               // the normal is left at zero; the reader recomputes
    for (const v of t) for (const x of v) { dv.setFloat32(o, x, true); o += 4; }
  });
  return new Uint8Array(buf);
}

const ASSEMBLY = JSON.stringify({
  title: 'one uploaded part',
  duration_s: 0.1, frames: 2,
  domains: [{
    kind: 'block', name: 'assembly', cells: [4, 4, 4], cell_mm: 5.0, initial_c: 20.0,
    parts: [{ stl: 'bracket.stl', material: 'copper' }],
  }],
});

// **A scene naming a file nobody uploaded is refused, and the refusal says what is here.**
// The failure this prevents is the one that wastes an afternoon: a misspelt name and a message
// that says only "not found", with three files sitting in the tab under other spellings.
forgetPart('');
const missing = call(w.dualis_check, ASSEMBLY);
ok('an un-uploaded part is refused', !!missing.error && missing.error.includes('bracket.stl'),
   JSON.stringify(missing.error));
ok('and the refusal says nothing has been uploaded',
   !!missing.error && missing.error.includes('nothing has been uploaded'), missing.error);

// **Uploading reports what the module now holds**, which is what the page's list shows.
const stl = boxStl([0, 0, 0], [20, 20, 20]);
const added = sendPart('bracket.stl', stl);
ok('an upload is acknowledged by name', !added.error && added.names.length === 1
   && added.names[0] === 'bracket.stl', JSON.stringify(added));
ok('and by size', added.bytes === stl.length, `${added.bytes} against ${stl.length}`);

// **And now the same scene builds** -- same bytes, same builder, same voxeliser as the CLI.
const built = call(w.dualis_check, ASSEMBLY);
ok('the scene builds once its part is here', !built.error, JSON.stringify(built.error));
ok('and the voxeliser reports what it cost',
   built.notes.some(n => n.includes('bracket.stl') && n.includes('filled')),
   JSON.stringify(built.notes));
ok('and the block is drawn', built.boxes.length === 1, String(built.boxes.length));

// **A 20 mm box on a 20 mm grid fills every cell**, the claim that says the bytes were read as
// millimetres and landed where the scene put them rather than merely parsing.
ok('a 20 mm part on a 20 mm block fills all 64 cells',
   built.notes.some(n => n.includes('filled 64 cells')), JSON.stringify(built.notes));

// **The run is the CLI's run.** A check that builds and a run that does not is the failure worth
// separating, because the page shows the first and the user asked for the second.
const ranOut = take(w.dualis_run());
ok('an uploaded assembly runs', !ranOut.error, JSON.stringify(ranOut.error));

// **Forgetting is real**, or a page that lets somebody remove a file would keep running the old
// bytes and show a stale answer as a fresh one.
const cleared = forgetPart('bracket.stl');
ok('forgetting empties the list', !cleared.error && cleared.names.length === 0, JSON.stringify(cleared));
const gone = call(w.dualis_check, ASSEMBLY);
ok('and the scene is refused again', !!gone.error && gone.error.includes('bracket.stl'), gone.error);

// **An empty upload is refused where it can be named.** A zero-byte file is what a failed read in
// the page looks like; `Mesh::from_stl` would call it a malformed mesh, which is true and about
// the wrong thing.
const empty = sendPart('nothing.stl', new Uint8Array(0));
ok('an empty upload is refused as an upload', !!empty.error && empty.error.includes('empty'),
   JSON.stringify(empty.error));


// --- fitting a grid, which is the step between "I dropped a file on this" and "here is a scene"
// A user who uploads CAD has no way to know what `cells` and `cell_mm` hold it, and getting it
// wrong does not fail -- a part finer than the grid rasterises to nothing and the run is well
// behaved about a different object.

function fitGrid(budget, material) {
  const mb = enc.encode(material);
  const mp = w.dualis_alloc(Math.max(1, mb.length));
  if (mb.length) new Uint8Array(mem.buffer, mp, mb.length).set(mb);
  const out = take(w.dualis_fit(budget, mp, mb.length));
  w.dualis_free(mp, Math.max(1, mb.length));
  return out;
}

forgetPart('');
sendPart('bracket.stl', boxStl([0, 0, 0], [45, 30, 4]));
sendPart('boss.stl', boxStl([18, 12, 4], [27, 18, 12]));

const fitted = fitGrid(400000, 'aluminium');
ok('a fit comes back with a table', !fitted.error && fitted.table.includes('cell_mm'),
   JSON.stringify(fitted.error || fitted.table.slice(0, 60)));

// **The extent is the union of the parts**, which is the first thing that would be wrong if the
// two files were being fitted separately rather than onto one grid.
ok('the assembly is the union of its parts',
   fitted.table.includes('45.0 x 30.0 x 12.0 mm'), fitted.table.split('\n')[0]);

// **The ladder resolves the thinnest feature by powers of two.** The boss is 4 mm through, so the
// rows are 4, 2, 1, 0.5 mm -- each one a sentence about the assembly rather than a round number.
for (const mm of ['4.000', '2.000', '1.000', '0.500']) {
  ok(`the ladder has a ${mm} mm row`, fitted.table.includes(mm), fitted.table);
}

// **A coarse row loses the small part and the table says so**, which is the failure this whole
// step exists to prevent: at 4 mm the boss fills four cells and is 40% short by volume.
const coarse = fitted.table.split('\n').find(l => l.includes('4.000'));
ok('the coarsest row is charged 100% boundary', coarse.includes('100.0'), coarse);

// **The recommendation is a scene fragment that names every uploaded part.**
ok('a fragment comes back', !!fitted.fragment, JSON.stringify(fitted.fragment));
ok('and it names both parts',
   fitted.fragment.includes('bracket.stl') && fitted.fragment.includes('boss.stl'),
   fitted.fragment);
ok('and carries the grid it measured',
   fitted.fragment.includes(`"cells": [${fitted.cells.join(', ')}]`), fitted.fragment);

// **And the fragment builds.** A table nobody can act on is a table; this is the round trip from
// two dropped files to a running scene with nothing in between to get wrong.
const assembled = call(w.dualis_check, `{ "title": "from the fitter", "duration_s": 1.0, "frames": 2,
  "domains": [{ "kind": "block", "name": "assembly", "initial_c": 20.0,
${fitted.fragment} }] }`);
ok('the suggested grid builds', !assembled.error, JSON.stringify(assembled.error));
ok('with both parts in it',
   assembled.notes.filter(n => n.includes('filled')).length === 2, JSON.stringify(assembled.notes));

// **A budget too small to hold anything is an answer, not a crash.** `fragment` is null and the
// table still says what each row cost, which is what a user needs in order to raise it.
const starved = fitGrid(50, 'aluminium');
ok('an impossible budget still returns a table', !starved.error && !!starved.table,
   JSON.stringify(starved.error));
ok('and recommends nothing rather than guessing', starved.fragment === null,
   JSON.stringify(starved.fragment));

forgetPart('');
const nothing = fitGrid(400000, 'aluminium');
ok('fitting nothing is refused by name', !!nothing.error && nothing.error.includes('at least one'),
   JSON.stringify(nothing.error));

console.log(failures ? `\n${failures} failed` : '\nall claims held');
process.exit(failures ? 1 : 0);
