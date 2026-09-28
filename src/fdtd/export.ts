import type { Fdtd2D } from './solver';

/** Raw results: config + dt + port V/I/Vs time series + probe samples; enough for independent post-processing. */
export function exportJSON(sim: Fdtd2D): string {
  const n = sim.step;
  return JSON.stringify({
    format: 'tinyem2d-v1',
    notes: 'Port samples at t=(n+0.5)*dt. V = port voltage, I = (Vs-V)/R flowing into structure. Probe sample s at t=(s+1)*every*dt (E) / (s+1.5)*every*dt approx (H at half step). Probe data layout [sample][point][component].',
    config: sim.cfg,
    dt: sim.dt, dtMax: sim.dtMax, stepsCompleted: n,
    ports: sim.ports.map((p) => ({
      name: p.name, component: p.comp, resistance: p.R, gapLength: p.len,
      t: Array.from({ length: n }, (_, k) => (k + 0.5) * sim.dt),
      Vs: Array.from(p.vs.subarray(0, n)), V: Array.from(p.v.subarray(0, n)), I: Array.from(p.i.subarray(0, n)),
    })),
    probes: sim.probes.map((p) => {
      const ns = Math.min(p.nSamples, Math.floor(n / p.every));
      return {
        name: p.name, components: p.components, every: p.every, points: p.coords,
        t: Array.from({ length: ns }, (_, s) => (s + 1) * p.every * sim.dt),
        data: Array.from(p.data.subarray(0, ns * p.idx.length * p.components.length)),
      };
    }),
  });
}

export function exportPortCSV(sim: Fdtd2D): string {
  const n = sim.step;
  const head = ['t', ...sim.ports.flatMap((p) => [`${p.name}_Vs`, `${p.name}_V`, `${p.name}_I`, `${p.name}_P`])];
  const rows = [head.join(',')];
  for (let k = 0; k < n; k++) {
    const r = [((k + 0.5) * sim.dt).toExponential(9)];
    for (const p of sim.ports) r.push(p.vs[k].toExponential(9), p.v[k].toExponential(9), p.i[k].toExponential(9), (p.v[k] * p.i[k]).toExponential(9));
    rows.push(r.join(','));
  }
  return rows.join('\n');
}

export function download(name: string, text: string, type: string): void {
  const url = URL.createObjectURL(new Blob([text], { type }));
  const a = document.createElement('a');
  a.href = url; a.download = name; a.click();
  URL.revokeObjectURL(url);
}
 
