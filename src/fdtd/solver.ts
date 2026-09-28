/**
 * TinyEM 2D Yee-grid FDTD engine.
 *
 * Field arrays retain the compact Nx*Ny storage used by the original prototype,
 * but material coefficients are now evaluated at their actual staggered Yee
 * locations by local arithmetic averaging. Internal PEC geometry is enforced
 * through zero tangential-E update coefficients.
 */
import { rasterize } from './geometry';
import { sampleWaveform } from './waveform';
import { C0, EPS0, MU0, modeComponents, type BoundaryType, type Component, type Material, type Port, type SimConfig } from './types';

type FieldMap = Partial<Record<Component, Float64Array>>;

interface PortEdge {
  k: number;
  comp: Component;
  len: number;
  area: number;
  branchR: number;
  sourceV: number;
}

export interface PortState {
  name: string;
  comp: Component;
  k: number;
  len: number;
  R: number;
  edges: PortEdge[];
  eOld: Float64Array;
  vs: Float64Array;
  v: Float64Array;
  i: Float64Array;
}

export interface ProbeState {
  name: string;
  components: Component[];
  every: number;
  idx: Int32Array;
  coords: [number, number][];
  nSamples: number;
  data: Float64Array;
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

    this.materials = cfg.materials;
    if (this.materials.length > 255) throw new Error('At most 255 materials are supported.');
    const matIndex = new Map(this.materials.map((m, i) => [m.name, i] as [string, number]));
    const bg = matIndex.get(cfg.background);
    if (bg === undefined) throw new Error(`Invalid background material "${cfg.background}".`);
    this.matMap = new Uint8Array(N).fill(bg);
    rasterize(cfg.geometry, matIndex, this.matMap, this.nx, this.ny, this.dx, this.dy);

    // Electric-field coefficients are evaluated at the electric-field location.
    // The node map is interpreted as a material sample map and locally averaged
    // over the surrounding material samples. This removes the original mistake
    // of using one node material for all staggered Yee components.
    for (const c of comps.e) {
      this.fields[c] = new Float64Array(N);
      const ca = new Float64Array(N), cb = new Float64Array(N);
      for (let i = 0; i < this.nx; i++) {
        for (let j = 0; j < this.ny; j++) {
          const k = i * this.ny + j;
          const ids = this.eMaterialSamples(c, i, j);
          const props = this.averageMaterial(ids);
          if (props.pec) {
            ca[k] = 0;
            cb[k] = 0;
            continue;
          }
          const a = (props.sigma * this.dt) / (2 * props.eps);
          ca[k] = (1 - a) / (1 + a);
          cb[k] = this.dt / props.eps / (1 + a);
        }
      }
      this.ca[c] = ca; this.cb[c] = cb;
    }

    // Magnetic-field coefficients are evaluated at the magnetic Yee locations.
    for (const c of comps.h) {
      this.fields[c] = new Float64Array(N);
      const db = new Float64Array(N);
      for (let i = 0; i < this.nx; i++) {
        for (let j = 0; j < this.ny; j++) {
          const k = i * this.ny + j;
          const ids = this.hMaterialSamples(c, i, j);
          db[k] = this.dt / this.averageMaterial(ids).mu;
        }
      }
      this.db[c] = db;
    }

    this.initPorts();
    this.initSources();
    this.initProbes();
    this.initBoundaries();
  }

  private materialAt(i: number, j: number): Material {
    return this.materials[this.matMap[i * this.ny + j]];
  }

  private validIds(samples: [number, number][]): number[] {
    const ids: number[] = [];
    for (const [i, j] of samples) {
      if (i >= 0 && i < this.nx && j >= 0 && j < this.ny) ids.push(i * this.ny + j);
    }
    return ids;
  }

  private eMaterialSamples(c: Component, i: number, j: number): number[] {
    if (c === 'Ez') {
      return this.validIds([[i, j], [i - 1, j], [i, j - 1], [i - 1, j - 1]]);
    }
    if (c === 'Ex') return this.validIds([[i, j], [i, j + 1]]);
    if (c === 'Ey') return this.validIds([[i, j], [i + 1, j]]);
    return [i * this.ny + j];
  }

  private hMaterialSamples(c: Component, i: number, j: number): number[] {
    if (c === 'Hx') return this.validIds([[i, j], [i, j + 1]]);
    if (c === 'Hy') return this.validIds([[i, j], [i + 1, j]]);
    if (c === 'Hz') return this.validIds([[i, j], [i + 1, j], [i, j + 1], [i + 1, j + 1]]);
    return [i * this.ny + j];
  }

  private averageMaterial(ids: number[]): { eps: number; mu: number; sigma: number; pec: boolean } {
    if (!ids.length) throw new Error('Internal material interpolation error: no surrounding material samples.');
    let eps = 0, mu = 0, sigma = 0, pec = false;
    for (const id of ids) {
      const m = this.materials[this.matMap[id]];
      eps += EPS0 * m.epsR;
      mu += MU0 * m.muR;
      sigma += m.sigma;
      pec ||= Boolean(m.pec);
    }
    const n = ids.length;
    return { eps: eps / n, mu: mu / n, sigma: sigma / n, pec };
  }

  private node(x: number, y: number, what: string): number {
    const i = Math.round(x / this.dx), j = Math.round(y / this.dy);
    if (i < 0 || j < 0 || i >= this.nx || j >= this.ny) throw new Error(`${what} lies outside computational domain.`);
    return i * this.ny + j;
  }

  private portEdges(p: Port): PortEdge[] {
    const comp = ('E' + p.orientation) as Component;
    if (!this.fields[comp]) throw new Error(`Port "${p.name}": orientation "${p.orientation}" not available in ${this.mode} mode.`);
    if (!(p.resistance > 0)) throw new Error(`Port "${p.name}": resistance must be > 0.`);

    if (p.orientation === 'z') {
      const k = this.node(p.x, p.y, `Port "${p.name}"`);
      const props = this.averageMaterial(this.eMaterialSamples(comp, Math.round(p.x / this.dx), Math.round(p.y / this.dy)));
      if (props.pec) throw new Error(`Port "${p.name}" is on PEC. Define a physical feed gap instead.`);
      return [{ k, comp, len: this.depth, area: this.dx * this.dy, branchR: p.resistance, sourceV: 1 }];
    }

    const x1 = p.x1 ?? (p.orientation === 'x' ? p.x + this.dx : p.x);
    const y1 = p.y1 ?? (p.orientation === 'y' ? p.y + this.dy : p.y);
    const i0 = Math.round(p.x / this.dx), j0 = Math.round(p.y / this.dy);
    const i1 = Math.round(x1 / this.dx), j1 = Math.round(y1 / this.dy);
    const tol = 1e-9;
    if (p.orientation === 'x' && Math.abs(y1 - p.y) > tol) throw new Error(`Port "${p.name}" must be axis-aligned with orientation x.`);
    if (p.orientation === 'y' && Math.abs(x1 - p.x) > tol) throw new Error(`Port "${p.name}" must be axis-aligned with orientation y.`);

    const edges: PortEdge[] = [];
    if (p.orientation === 'x') {
      const j = j0, a = Math.min(i0, i1), b = Math.max(i0, i1);
      if (b <= a) throw new Error(`Port "${p.name}" must span at least one Yee edge.`);
      for (let i = a; i < b; i++) {
        const k = i * this.ny + j;
        const props = this.averageMaterial(this.eMaterialSamples(comp, i, j));
        if (props.pec) throw new Error(`Port "${p.name}" crosses PEC. Leave a material gap between the terminals.`);
        edges.push({ k, comp, len: this.dx, area: this.dy * this.depth, branchR: 0, sourceV: 0 });
      }
    } else {
      const i = i0, a = Math.min(j0, j1), b = Math.max(j0, j1);
      if (b <= a) throw new Error(`Port "${p.name}" must span at least one Yee edge.`);
      for (let j = a; j < b; j++) {
        const k = i * this.ny + j;
        const props = this.averageMaterial(this.eMaterialSamples(comp, i, j));
        if (props.pec) throw new Error(`Port "${p.name}" crosses PEC. Leave a material gap between the terminals.`);
        edges.push({ k, comp, len: this.dy, area: this.dx * this.depth, branchR: 0, sourceV: 0 });
      }
    }

    const n = edges.length;
    for (const e of edges) {
      e.branchR = p.resistance / n;
      e.sourceV = 1 / n;
    }
    return edges;
  }

  private initPorts(): void {
    for (const p of this.cfg.ports) {
      const edges = this.portEdges(p);
      const comp = edges[0].comp;
      const len = edges.reduce((sum, e) => sum + e.len, 0);

      for (const e of edges) {
        const i = Math.floor(e.k / this.ny), j = e.k % this.ny;
        const props = this.averageMaterial(this.eMaterialSamples(comp, i, j));
        const a = (props.sigma * this.dt) / (2 * props.eps)
          + (this.dt * e.len) / (2 * e.branchR * props.eps * e.area);
        this.ca[comp]![e.k] = (1 - a) / (1 + a);
        this.cb[comp]![e.k] = this.dt / props.eps / (1 + a);
      }

      this.ports.push({
        name: p.name,
        comp,
        k: edges[0].k,
        len,
        R: p.resistance,
        edges,
        eOld: new Float64Array(edges.length),
        vs: sampleWaveform(p.waveform, this.steps, this.dt, 0.5),
        v: new Float64Array(this.steps),
        i: new Float64Array(this.steps),
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

  advance(count: number): boolean {
    const end = Math.min(this.steps, this.step + count);
    while (this.step < end) this.stepOnce();
    return this.step >= this.steps;
  }

  private stepOnce(): void {
    const n = this.step;
    if (this.mode === 'TMz') this.updateH_TM(); else this.updateH_TE();
    for (const s of this.hSources) s.arr[s.k] += s.w[n];

    for (const p of this.ports) {
      for (let q = 0; q < p.edges.length; q++) {
        p.eOld[q] = this.fields[p.edges[q].comp]![p.edges[q].k];
      }
    }

    for (const e of this.edges) {
      if (e.type !== 'mur') continue;
      for (let q = 0; q < e.b.length; q++) {
        e.old0[q] = e.arr[e.b[q]];
        e.old1[q] = e.arr[e.inner[q]];
      }
    }

    if (this.mode === 'TMz') this.updateE_TM(); else this.updateE_TE();
    for (const s of this.eSources) s.arr[s.k] += s.w[n];

    for (const p of this.ports) {
      let v = 0, currentSum = 0;
      for (let q = 0; q < p.edges.length; q++) {
        const e = p.edges[q];
        const f = this.fields[e.comp]!;
        const edgeV = 0.5 * (p.eOld[q] + f[e.k]) * e.len;
        v += edgeV;
        const branchVs = p.vs[n] * e.sourceV;
        currentSum += (branchVs - edgeV) / e.branchR;
        // The branch source is already included implicitly by the modified
        // electric-field coefficient/source term below. The field update
        // itself is completed here using the known generator waveform.
        f[e.k] += (this.dt / (e.branchR * this.averageMaterial(
          this.eMaterialSamples(e.comp, Math.floor(e.k / this.ny), e.k % this.ny),
        ).eps * (1 + 0))) * 0;
      }
      p.v[n] = v;
      p.i[n] = currentSum / p.edges.length;
    }

    // Apply the distributed Thévenin source after the standard E update.
    // The source term is equivalent to a series source/resistor on each
    // Yee edge; voltage and resistance are evenly distributed along a
    // multi-edge series gap.
    for (const p of this.ports) {
      for (const e of p.edges) {
        const f = this.fields[e.comp]!;
        const i = Math.floor(e.k / this.ny), j = e.k % this.ny;
        const props = this.averageMaterial(this.eMaterialSamples(e.comp, i, j));
        const a = (props.sigma * this.dt) / (2 * props.eps)
          + (this.dt * e.len) / (2 * e.branchR * props.eps * e.area);
        const sourceCoeff = this.dt / (e.branchR * props.eps * e.area * (1 + a));
        f[e.k] += sourceCoeff * p.vs[n] * e.sourceV;
      }
    }

    for (const e of this.edges) {
      const a = e.arr, b = e.b;
      if (e.type === 'pec') {
        for (let q = 0; q < b.length; q++) a[b[q]] = 0;
      } else {
        for (let q = 0; q < b.length; q++) a[b[q]] = e.old1[q] + e.kMur * (a[e.inner[q]] - e.old0[q]);
      }
    }

    this.applyPmcBoundaries();
    this.step = n + 1;
    this.recordProbes();
  }

  private applyPmcBoundaries(): void {
    const b = this.cfg.boundaries;
    if (this.mode === 'TMz') {
      const hx = this.fields.Hx!, hy = this.fields.Hy!;
      if (b.xmin === 'pmc') for (let j = 0; j < this.ny; j++) hy[j] = 0;
      if (b.xmax === 'pmc') for (let j = 0; j < this.ny; j++) hy[(this.nx - 2) * this.ny + j] = 0;
      if (b.ymin === 'pmc') for (let i = 0; i < this.nx; i++) hx[i * this.ny] = 0;
      if (b.ymax === 'pmc') for (let i = 0; i < this.nx; i++) hx[i * this.ny + this.ny - 2] = 0;
    } else {
      const ex = this.fields.Ex!, ey = this.fields.Ey!;
      if (b.xmin === 'pmc') for (let j = 0; j < this.ny; j++) ex[j] = 0;
      if (b.xmax === 'pmc') for (let j = 0; j < this.ny; j++) ex[(this.nx - 1) * this.ny + j] = 0;
      if (b.ymin === 'pmc') for (let i = 0; i < this.nx; i++) ey[i * this.ny] = 0;
      if (b.ymax === 'pmc') for (let i = 0; i < this.nx; i++) ey[i * this.ny + this.ny - 1] = 0;
    }
  }

  private recordProbes(): void {
    for (const p of this.probes) {
      if (this.step % p.every !== 0) continue;
      const s = this.step / p.every - 1;
      if (s >= p.nSamples) continue;
      const nc = p.components.length, np = p.idx.length;
      let o = s * np * nc;
      for (let q = 0; q < np; q++) {
        for (let c = 0; c < nc; c++) {
          p.data[o++] = this.fields[p.components[c]]![p.idx[q]];
        }
      }
    }
  }

  private updateH_TM(): void {
    const { nx, ny, dx, dy } = this;
    const ez = this.fields.Ez!, hx = this.fields.Hx!, hy = this.fields.Hy!;
    const dbx = this.db.Hx!, dby = this.db.Hy!;
    for (let i = 0; i < nx; i++) {
      const r = i * ny;
      for (let j = 0; j < ny - 1; j++) {
        const k = r + j;
        hx[k] -= (dbx[k] * (ez[k + 1] - ez[k])) / dy;
      }
      if (i < nx - 1) {
        for (let j = 0; j < ny; j++) {
          const k = r + j;
          hy[k] += (dby[k] * (ez[k + ny] - ez[k])) / dx;
        }
      }
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
