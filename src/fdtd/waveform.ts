import type { Waveform } from './types';

export function evalWaveform(w: Waveform, t: number): number {
  switch (w.type) {
    case 'gaussian': {
      const u = (t - w.t0) / w.width;
      return w.amplitude * Math.exp(-u * u);
    }
    case 'ricker': {
      const a = Math.PI * w.f0 * (t - w.delay);
      const a2 = a * a;
      return w.amplitude * (1 - 2 * a2) * Math.exp(-a2);
    }
    case 'cw': {
      const ramp = w.ramp && w.ramp > 0 ? Math.min(1, t / w.ramp) : 1;
      return w.amplitude * ramp * Math.sin(2 * Math.PI * w.frequency * t + (w.phase ?? 0));
    }
  }
}

/** Samples waveform at t = (n + offset) * dt for n = 0..steps-1 (preallocated before the time loop). */
export function sampleWaveform(w: Waveform, steps: number, dt: number, offset: number): Float64Array {
  const out = new Float64Array(steps);
  for (let n = 0; n < steps; n++) out[n] = evalWaveform(w, (n + offset) * dt);
  return out;
}
 
