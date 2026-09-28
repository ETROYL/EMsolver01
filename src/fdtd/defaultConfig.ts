 import type { SimConfig } from './types';

// Example: TEz strip dipole (PEC line, y-directed lumped port at centre) next to a dielectric cylinder.
export const defaultConfig: SimConfig = {
  simulation: { mode: 'TEz', nx: 200, ny: 200, dx: 0.001, dy: 0.001, dt: null, courant: 0.99, steps: 4000, depth: 0.001 },
  boundaries: { xmin: 'mur', xmax: 'mur', ymin: 'mur', ymax: 'mur' },
  background: 'air',
  materials: [
    { name: 'air', epsR: 1, muR: 1, sigma: 0 },
    { name: 'copper', epsR: 1, muR: 1, sigma: 0, pec: true },
    { name: 'dielectric', epsR: 4.4, muR: 1, sigma: 0.002 },
  ],
  geometry: [
    { type: 'line', material: 'copper', x0: 0.1, y0: 0.07, x1: 0.1, y1: 0.13 },
    { type: 'circle', material: 'dielectric', cx: 0.145, cy: 0.1, r: 0.015 },
  ],
  sources: [],
  ports: [
    { name: 'P1', x: 0.1, y: 0.1, orientation: 'y', resistance: 50, waveform: { type: 'gaussian', amplitude: 1, t0: 2.4e-10, width: 6e-11 } },
  ],
  probes: [
    { name: 'pt1', type: 'point', x: 0.05, y: 0.1, components: ['Ey', 'Hz'], every: 1 },
    { name: 'line1', type: 'line', x: 0.02, y: 0.03, x1: 0.18, y1: 0.03, n: 81, components: ['Ex', 'Ey'], every: 10 },
  ],
  output: { z0: 50, fmin: 2e8, fmax: 1e10, fftPad: 4, snapshotEvery: 10 },
};

