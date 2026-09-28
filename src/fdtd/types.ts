// Configuration schema for the TinyEM 2D FDTD engine. All lengths in metres, times in seconds, SI units.

export type Mode = 'TMz' | 'TEz';
export type BoundaryType = 'pec' | 'pmc' | 'mur';
export type Component = 'Ex' | 'Ey' | 'Ez' | 'Hx' | 'Hy' | 'Hz';

export interface Material {
  name: string;
  epsR: number;
  muR: number;
  sigma: number;
  pec?: boolean;
}

export type Shape =
  | { type: 'box'; material: string; xmin: number; xmax: number; ymin: number; ymax: number }
  | { type: 'circle'; material: string; cx: number; cy: number; r: number }
  | { type: 'polygon'; material: string; points: [number, number][] }
  | { type: 'line'; material: string; x0: number; y0: number; x1: number; y1: number; width?: number }
  | { type: 'cells'; material: string; cells: [number, number][] };

export type Waveform =
  | { type: 'gaussian'; amplitude: number; t0: number; width: number }
  | { type: 'ricker'; amplitude: number; f0: number; delay: number }
  | { type: 'cw'; amplitude: number; frequency: number; phase?: number; ramp?: number };

export interface Source {
  component: Component;
  x: number;
  y: number;
  waveform: Waveform;
}

export interface Port {
  name: string;
  x: number;
  y: number;
  orientation: 'x' | 'y' | 'z';
  resistance: number;
  waveform: Waveform;
}

export interface Probe {
  name: string;
  type: 'point' | 'line';
  x: number;
  y: number;
  x1?: number;
  y1?: number;
  n?: number;
  components: Component[];
  every?: number;
}

export interface SimConfig {
  simulation: {
    mode: Mode;
    nx: number;
    ny: number;
    dx: number;
    dy: number;
    dt?: number | null;
    courant?: number;
    steps: number;
    depth?: number;
  };
  boundaries: { xmin: BoundaryType; xmax: BoundaryType; ymin: BoundaryType; ymax: BoundaryType };
  background: string;
  materials: Material[];
  geometry: Shape[];
  sources: Source[];
  ports: Port[];
  probes: Probe[];
  output: { z0: number; fmin: number; fmax: number; fftPad?: number; snapshotEvery: number };
}

export const C0 = 299792458;
export const MU0 = 4e-7 * Math.PI;
export const EPS0 = 1 / (MU0 * C0 * C0);

export function modeComponents(mode: Mode): { e: Component[]; h: Component[] } {
  return mode === 'TMz' ? { e: ['Ez'], h: ['Hx', 'Hy'] } : { e: ['Ex', 'Ey'], h: ['Hz'] };
}
 
