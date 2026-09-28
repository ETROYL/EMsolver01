/**
 * TinyEM 2D Yee-grid FDTD engine. Framework-free; only typed arrays in the time loop.
 *
 * Index k = i*ny + j. Node positions (Yee staggering):
 *  TMz: Ez(i,j), Hx(i,j+1/2), Hy(i+1/2,j)
 *  TEz: Ex(i+1/2,j), Ey(i,j+1/2), Hz(i+1/2,j+1/2)
 * H components beyond the last staggered node are never updated and stay 0 (acts as PMC);
 * PEC / Mur are applied to tangential E on each side after the E update.
 * Time levels: E at n*dt, H at (n+1/2)*dt.
 */
import { rasterize } from './geometry';
import { sampleWaveform } from './waveform';
import { C0, EPS0, MU0, modeComponents, type BoundaryType, type Component, type Material, type SimConfig } from './types';

type FieldMap = Partial<Record<Component, Float64Array>>;

export interface PortState {
  name: string;
  comp: Component;
  k: number;
  len: number;
  R: number;
  cs: number;
  eOld: number;
  vs: Float64Array; // source voltage at (n+1/2)dt
  v: Float64Array; // port voltage at (n+1/2)dt (average of E^n, E^{n+1})
  i: Float64Array; // port current (Vs - V)/R at (n+1/2)dt
}

export interface ProbeState {
  name: string;
  components: Component[];
  every: number;
  idx: Int32Array;
  coords: [number, number][];
  nSamples: number;
  data: Float64Array; // [sample][point][component]
}

interface SourceState { arr: Float64Array; k: number; w: Float64Array; }

interface EdgeState {
  type: BoundaryType;
  arr: Float64Array;
  b: Int32Array;
  inner: Int32Array;
  old0: Float64Array;
  old1: Float64Array;
  kMur: number;
}

export class Fdtd2D {
  readonly cfg: SimConfig;
  readonly nx: number; readonly ny: number; readonly dx: number; readonly dy: number;
  readonly dt: number; readonly dtMax: number; readonly steps: number; readonly depth: number;
  readonly mode: SimConfig['simulation']['mode'];
  readonly materials: Material[];
  readonly matMap: Uint8Array;
  readonly fields: FieldMap = {};
  readonly warnings: string[] = [];
  readonly ports: PortState[] = [];
  readonly probes: ProbeState[] = [];
  readonly memoryBytes: number;
  step = 0;

  private readonly ca: FieldMap = {};
  private readonly cb: FieldMap = {};
  private readonly db: FieldMap = {};
  private readonly eSources: SourceState[] = [];
  private readonly hSources: SourceState[] = [];
  private readonly edges: EdgeState[] = [];

  constructor(cfg: SimConfig) {
    this.cfg = cfg;
    const s = cfg.simulation;
    if (s.mode !== 'TMz' && s.mode !== 'TEz') throw new Error('Unsupported solver mode (use "TMz" or "TEz").');
    if (!(s.nx > 2 && s.ny > 2)) throw new Error('Grid dimensions must be greater than 2.');
    if (!(s.dx > 0 && s.dy > 0)) throw new Error('Cell size must be greater than zero.');
    if (!(s.steps > 0)) throw new Error('Number of timesteps must be greater than zero.');
    this.mode = s.mode;
    this.nx = Math.floor(s.nx); this.ny = Math.floor(s.ny);
    this.dx = s.dx; this.dy = s.dy; this.steps = Math.floor(s.steps);
    this.depth = s.depth ?? s.dx;
    this.dtMax = 1 / (C0 * Math.sqrt(1 / (s.dx * s.dx) + 1 / (s.dy * s.dy)));
    this.dt = s.dt && s.dt > 0 ? s.dt : (s.courant ?? 0.99) * this.dtMax;
    if (this.dt > this.dtMax) this.warnings.push(`CFL stability condition violated: dt=${this.dt.toExponential(3)} > dtMax=${this.dtMax.toExponential(3)}.`);

    const N = this.nx * this.ny;
    const comps = modeComponents(this.mode);
    this.memoryBytes = N * 8 * (comps.e.length * 3 + comps.h.length * 2) + N;
    if (this.memoryBytes > 1.5e9) throw new Error(`Simulation too large: ~${(this.memoryBytes / 1e9).toFixed(2)} GB required.`);

    // Materials and voxelisation
    this.materials = cfg.materials;
    if (this.materials.length > 255) throw new Error('At most 255 materials are supported.');
    const matIndex = new Map(this.materials.map((m, i) => [m.name, i] as [string, number]));
    const bg = matIndex.get(cfg.background);
    if (bg === undefined) throw new Error(`Invalid background material "${cfg.background}".`);
    this.matMap = new Uint8Array(N).fill(bg);
    rasterize(cfg.geometry, matIndex, this.matMap, this.nx, this.ny, this.dx, this.dy);

    // Update coefficients (semi-implicit conductivity averaging)
    for (const c of comps.e) {
      this.fields[c] = new Float64Array(N);
      const ca = new Float64Array(N), cb = new Float64Array(N);
      for (let k = 0; k < N; k++) {
        const m = this.materials[this.matMap[k]];
        if (m.pec) continue;
        const eps = EPS0 * m.epsR;
        const a = (m.sigma * this.dt) / (2 * eps);
        ca[k] = (1 - a) / (1 + a);
        cb[k] = this.dt / eps / (1 + a);
      }
      this.ca[c] = ca; this.cb[c] = cb;
    }
    for (const c of comps.h) {
      this.fields[c] = new Float64Array(N);
      const db = new Float64Array(N);
      for (let k = 0; k < N; k++) db[k] = this.dt / (MU0 * this.materials[this.matMap[k]].muR);
      this.db[c] = db;
    }

    this.initPorts();
    this.initSources();
    this.initProbes();
    this.initBoundaries();
  }

  private node(x: number, y: number, what: string): number {
    const i = Math.round(x / this.dx), j = Math.round(y / this.dy);
    if (i < 0 || j < 0 || i >= this.nx || j >= this.ny) throw new Error(`${what} lies outside computational domain.`);
    return i * this.ny + j;
  }

  private initPorts(): void {
    for (const p of this.cfg.ports) {
      const comp = ('E' + p.orientation) as Component;
      if (!this.fields[comp]) throw new Error(`Port "${p.name}": orientation "${p.orientation}" not available in ${this.mode} mode.`);
      if (!(p.resistance > 0)) throw new Error(`Port "${p.name}": resistance must be > 0.`);
      const k = this.node(p.x, p.y, `Port "${p.name}"`);
      const len = p.orientation === 'x' ? this.dx : p.orientation === 'y' ? this.dy : this.depth;
      const area = p.orientation === 'x' ? this.dy * this.depth : p.orientation === 'y' ? this.dx * this.depth : this.dx * this.dy;
      const m = this.materials[this.matMap[k]];
      const eps = EPS0 * m.epsR;
      const sigma = m.pec ? 0 : m.sigma;
      // Resistive voltage source (Taflove lumped element): beta = dt*len/(2*R*eps*area)
      const a = (sigma * this.dt) / (2 * eps) + (this.dt * len) / (2 * p.resistance * eps * area);
      this.ca[comp]![k] = (1 - a) / (1 + a);
      this.cb[comp]![k] = this.dt / eps / (1 + a);
      this.ports.push({
        name: p.name, comp, k, len, R: p.resistance, eOld: 0,
        cs: this.dt / (p.resistance * eps * area) / (1 + a),
        vs: sampleWaveform(p.waveform, this.steps, this.dt, 0.5),
        v: new Float64Array(this.steps), i: new Float64Array(this.steps),
      });
    }
  }

  private initSources(): void {
    for (const s of this.cfg.sources) {
      const arr = this.fields[s.component];
      if (!arr) throw new Error(`Source component ${s.component} not available in ${this.mode} mode.`);
      const isH = s.component[0] === 'H';
      const st = { arr, k: this.node(s.x, s.y, 'Source'), w: sampleWaveform(s.waveform, this.steps, this.dt, isH ? 1 : 0.5) };
      (isH ? this.hSources : this.eSources).push(st);
    }
  }

  private initProbes(): void {
    for (const p of this.cfg.probes) {
      for (const c of p.components) if (!this.fields[c]) throw new Error(`Probe "${p.name}": component ${c} not available in ${this.mode} mode.`);
      const n = p.type === 'line' ? Math.max(2, Math.floor(p.n ?? 50)) : 1;
      const coords: [number, number][] = [];
      const idx = new Int32Array(n);
      for (let q = 0; q < n; q++) {
        const t = n === 1 ? 0 : q / (n - 1);
        const x = p.x + t * ((p.x1 ?? p.x) - p.x), y = p.y + t * ((p.y1 ?? p.y) - p.y);
        idx[q] = this.node(x, y, `Probe "${p.name}"`);
        coords.push([x, y]);
      }
      const every = Math.max(1, Math.floor(p.every ?? 1));
      const nSamples = Math.floor(this.steps / every);
      this.probes.push({ name: p.name, components: p.components, every, idx, coords, nSamples, data: new Float64Array(nSamples * n * p.components.length) });
    }
  }

  private initBoundaries(): void {
    const { nx, ny } = this;
    const sides: ['xmin' | 'xmax' | 'ymin' | 'ymax', Component][] =
      this.mode === 'TMz'
        ? [['xmin', 'Ez'], ['xmax', 'Ez'], ['ymin', 'Ez'], ['ymax', 'Ez']]
        : [['xmin', 'Ey'], ['xmax', 'Ey'], ['ymin', 'Ex'], ['ymax', 'Ex']];
    for (const [side, comp] of sides) {
      const type = this.cfg.boundaries[side];
      if (type === 'pmc') continue;
      if (type !== 'pec' && type !== 'mur') throw new Error(`Unknown boundary "${type}" on ${side}.`);
      const isX = side[0] === 'x';
      const len = isX ? ny : nx;
      const b = new Int32Array(len), inner = new Int32Array(len);
      for (let q = 0; q < len; q++) {
        if (side === 'xmin') { b[q] = q; inner[q] = ny + q; }
        else if (side === 'xmax') { b[q] = (nx - 1) * ny + q; inner[q] = (nx - 2) * ny + q; }
        else if (side === 'ymin') { b[q] = q * ny; inner[q] = q * ny + 1; }
        else { b[q] = q * ny + ny - 1; inner[q] = q * ny + ny - 2; }
      }
      const h = isX ? this.dx : this.dy;
      this.edges.push({
        type, arr: this.fields[comp]!, b, inner,
        old0: new Float64Array(len), old1: new Float64Array(len),
        kMur: (C0 * this.dt - h) / (C0 * this.dt + h),
      });
    }
  }

  /** Advances up to `count` timesteps. Returns true when finished. */
  advance(count: number): boolean {
    const end = Math.min(this.steps, this.step + count);
    while (this.step < end) this.stepOnce();
    return this.step >= this.steps;
  }

  private stepOnce(): void {
    const n = this.step;
    if (this.mode === 'TMz') this.updateH_TM(); else this.updateH_TE();
    for (const s of this.hSources) s.arr[s.k] += s.w[n];

    for (const p of this.ports) p.eOld = this.fields[p.comp]![p.k];
    for (const e of this.edges) {
      if (e.type !== 'mur') continue;
      for (let q = 0; q < e.b.length; q++) { e.old0[q] = e.arr[e.b[q]]; e.old1[q] = e.arr[e.inner[q]]; }
    }

    if (this.mode === 'TMz') this.updateE_TM(); else this.updateE_TE();
    for (const s of this.eSources) s.arr[s.k] += s.w[n];

    for (const p of this.ports) {
      const f = this.fields[p.comp]!;
      f[p.k] += p.cs * p.vs[n];
      const v = (0.5 * (p.eOld + f[p.k])) * p.len;
      p.v[n] = v;
      p.i[n] = (p.vs[n] - v) / p.R;
    }

    for (const e of this.edges) {
      const a = e.arr, b = e.b;
      if (e.type === 'pec') for (let q = 0; q < b.length; q++) a[b[q]] = 0;
      else for (let q = 0; q < b.length; q++) a[b[q]] = e.old1[q] + e.kMur * (a[e.inner[q]] - e.old0[q]);
    }

    this.step = n + 1;
    this.recordProbes();
  }

  private recordProbes(): void {
    for (const p of this.probes) {
      if (this.step % p.every !== 0) continue;
      const s = this.step / p.every - 1;
      if (s >= p.nSamples) continue;
      const nc = p.components.length, np = p.idx.length;
      let o = s * np * nc;
      for (let q = 0; q < np; q++) for (let c = 0; c < nc; c++) p.data[o++] = this.fields[p.components[c]]![p.idx[q]];
    }
  }

  // ---- TMz: Ez, Hx, Hy ----
  private updateH_TM(): void {
    const { nx, ny, dx, dy } = this;
    const ez = this.fields.Ez!, hx = this.fields.Hx!, hy = this.fields.Hy!;
    const dbx = this.db.Hx!, dby = this.db.Hy!;
    for (let i = 0; i < nx; i++) {
      const r = i * ny;
      for (let j = 0; j < ny - 1; j++) { const k = r + j; hx[k] -= (dbx[k] * (ez[k + 1] - ez[k])) / dy; }
      if (i < nx - 1) for (let j = 0; j < ny; j++) { const k = r + j; hy[k] += (dby[k] * (ez[k + ny] - ez[k])) / dx; }
    }
  }

  private updateE_TM(): void {
    const { nx, ny, dx, dy } = this;
    const ez = this.fields.Ez!, hx = this.fields.Hx!, hy = this.fields.Hy!;
    const ca = this.ca.Ez!, cb = this.cb.Ez!;
    const idx = 1 / dx, idy = 1 / dy;
    for (let i = 0; i < nx; i++) {
      const r = i * ny;
      for (let j = 0; j < ny; j++) {
        const k = r + j;
        const hyW = i > 0 ? hy[k - ny] : 0;
        const hxS = j > 0 ? hx[k - 1] : 0;
        ez[k] = ca[k] * ez[k] + cb[k] * ((hy[k] - hyW) * idx - (hx[k] - hxS) * idy);
      }
    }
  }

  // ---- TEz: Ex, Ey, Hz ----
  private updateH_TE(): void {
    const { nx, ny, dx, dy } = this;
    const ex = this.fields.Ex!, ey = this.fields.Ey!, hz = this.fields.Hz!, db = this.db.Hz!;
    const idx = 1 / dx, idy = 1 / dy;
    for (let i = 0; i < nx - 1; i++) {
      const r = i * ny;
      for (let j = 0; j < ny - 1; j++) {
        const k = r + j;
        hz[k] -= db[k] * ((ey[k + ny] - ey[k]) * idx - (ex[k + 1] - ex[k]) * idy);
      }
    }
  }

  private updateE_TE(): void {
    const { nx, ny, dx, dy } = this;
    const ex = this.fields.Ex!, ey = this.fields.Ey!, hz = this.fields.Hz!;
    const cax = this.ca.Ex!, cbx = this.cb.Ex!, cay = this.ca.Ey!, cby = this.cb.Ey!;
    const idx = 1 / dx, idy = 1 / dy;
    for (let i = 0; i < nx; i++) {
      const r = i * ny;
      for (let j = 0; j < ny; j++) {
        const k = r + j;
        if (i < nx - 1) ex[k] = cax[k] * ex[k] + cbx[k] * (hz[k] - (j > 0 ? hz[k - 1] : 0)) * idy;
        if (j < ny - 1) ey[k] = cay[k] * ey[k] - cby[k] * (hz[k] - (i > 0 ? hz[k - ny] : 0)) * idx;
      }
    }
  }
}
 
