/** FFT-based port post-processing: V(f), I(f), Z(f), Gamma(f), S11(f), VSWR. */

export function fft(re: Float64Array, im: Float64Array): void {
  const n = re.length;
  for (let i = 1, j = 0; i < n; i++) {
    let bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) { [re[i], re[j]] = [re[j], re[i]]; [im[i], im[j]] = [im[j], im[i]]; }
  }
  for (let len = 2; len <= n; len <<= 1) {
    const ang = (-2 * Math.PI) / len, wr = Math.cos(ang), wi = Math.sin(ang);
    for (let i = 0; i < n; i += len) {
      let cr = 1, ci = 0;
      for (let k = 0; k < len / 2; k++) {
        const a = i + k, b = a + len / 2;
        const tr = re[b] * cr - im[b] * ci, ti = re[b] * ci + im[b] * cr;
        re[b] = re[a] - tr; im[b] = im[a] - ti; re[a] += tr; im[a] += ti;
        const nr = cr * wr - ci * wi; ci = cr * wi + ci * wr; cr = nr;
      }
    }
  }
}

export interface PortSpectrum {
  f: number[]; vr: number[]; vi: number[]; ir: number[]; ii: number[];
  zr: number[]; zi: number[]; gammaMag: number[]; s11: number[]; vswr: number[];
}

/** Spectrum X(f) = sum x[n] e^{-j2πf n dt} dt, zero-padded to next power of 2 >= len*pad. */
export function spectrum(x: Float64Array, len: number, dt: number, pad: number): { re: Float64Array; im: Float64Array; df: number } {
  let N = 1;
  while (N < len * pad) N <<= 1;
  const re = new Float64Array(N), im = new Float64Array(N);
  for (let n = 0; n < len; n++) re[n] = x[n] * dt;
  fft(re, im);
  return { re, im, df: 1 / (N * dt) };
}

export function portSpectrum(v: Float64Array, i: Float64Array, len: number, dt: number, z0: number, fmin: number, fmax: number, pad = 4): PortSpectrum {
  const V = spectrum(v, len, dt, pad), I = spectrum(i, len, dt, pad);
  const out: PortSpectrum = { f: [], vr: [], vi: [], ir: [], ii: [], zr: [], zi: [], gammaMag: [], s11: [], vswr: [] };
  const k0 = Math.max(1, Math.ceil(fmin / V.df)), k1 = Math.min(V.re.length / 2, Math.floor(fmax / V.df));
  const stride = Math.max(1, Math.floor((k1 - k0) / 600));
  for (let k = k0; k <= k1; k += stride) {
    const vr = V.re[k], vi = V.im[k], ir = I.re[k], ii = I.im[k];
    const d = ir * ir + ii * ii || 1e-300;
    const zr = (vr * ir + vi * ii) / d, zi = (vi * ir - vr * ii) / d;
    // Gamma = (Z - Z0) / (Z + Z0)
    const nr = zr - z0, dr = zr + z0;
    const g = Math.sqrt((nr * nr + zi * zi) / (dr * dr + zi * zi));
    out.f.push(k * V.df); out.vr.push(vr); out.vi.push(vi); out.ir.push(ir); out.ii.push(ii);
    out.zr.push(zr); out.zi.push(zi); out.gammaMag.push(g);
    out.s11.push(20 * Math.log10(Math.max(g, 1e-12)));
    out.vswr.push(g < 1 ? (1 + g) / (1 - g) : Infinity);
  }
  return out;
}
 
