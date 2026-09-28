 import { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import { CartesianGrid, Line, LineChart, ResponsiveContainer, Tooltip, XAxis, YAxis } from 'recharts';
import { Fdtd2D } from '@/fdtd/solver';
import { defaultConfig } from '@/fdtd/defaultConfig';
import { portSpectrum, type PortSpectrum } from '@/fdtd/postprocess';
import { download, exportJSON, exportPortCSV } from '@/fdtd/export';
import type { Component, SimConfig } from '@/fdtd/types';

type View = Component | '|E|' | '|H|' | 'material';
const MAT_COLORS = [[20, 22, 28], [230, 170, 60], [80, 170, 220], [120, 220, 140], [220, 90, 150], [180, 180, 180]];

export default function Index() {
  const [text, setText] = useState(JSON.stringify(defaultConfig, null, 2));
  const [error, setError] = useState('');
  const [sim, setSim] = useState<Fdtd2D | null>(null);
  const [running, setRunning] = useState(false);
  const [step, setStep] = useState(0);
  const [view, setView] = useState<View>('Ey');
  const [spec, setSpec] = useState<PortSpectrum | null>(null);
  const [portIdx, setPortIdx] = useState(0);
  const [perf, setPerf] = useState('');
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const scaleRef = useRef(1e-9);

  const build = useCallback(() => {
    try {
      const cfg = JSON.parse(text) as SimConfig;
      const s = new Fdtd2D(cfg);
      setSim(s); setStep(0); setSpec(null); setError(s.warnings.join('\n')); setRunning(false); setPortIdx(0);
      scaleRef.current = 1e-9;
      const avail: View[] = cfg.simulation.mode === 'TMz' ? ['Ez'] : ['Ey'];
      setView(avail[0]);
      return s;
    } catch (e) {
      setError(`ERROR: ${(e as Error).message}`); setSim(null);
      return null;
    }
  }, [text]);

  useEffect(() => { build(); }, []); // eslint-disable-line react-hooks/exhaustive-deps

  const draw = useCallback((s: Fdtd2D) => {
    const cv = canvasRef.current; if (!cv) return;
    const { nx, ny } = s;
    if (cv.width !== nx || cv.height !== ny) { cv.width = nx; cv.height = ny; }
    const ctx = cv.getContext('2d')!;
    const img = ctx.createImageData(nx, ny);
    const d = img.data, mm = s.matMap, bgIdx = s.materials.findIndex((m) => m.name === s.cfg.background);
    const f = s.fields;
    let vals: Float64Array | null = null, mag = false;
    if (view === '|E|' || view === '|H|') {
      const cs = (view === '|E|' ? ['Ex', 'Ey', 'Ez'] : ['Hx', 'Hy', 'Hz']) as Component[];
      const arrs = cs.map((c) => f[c]).filter(Boolean) as Float64Array[];
      vals = new Float64Array(nx * ny);
      for (const a of arrs) for (let k = 0; k < vals.length; k++) vals[k] += a[k] * a[k];
      for (let k = 0; k < vals.length; k++) vals[k] = Math.sqrt(vals[k]);
      mag = true;
    } else if (view !== 'material') vals = f[view] ?? null;
    let mx = 0;
    if (vals) for (let k = 0; k < vals.length; k++) { const a = Math.abs(vals[k]); if (a > mx) mx = a; }
    scaleRef.current = Math.max(mx, scaleRef.current * 0.97, 1e-30);
    const sc = scaleRef.current;
    for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) {
      const k = i * ny + j, p = ((ny - 1 - j) * nx + i) * 4;
      let r = 12, g = 14, b = 20;
      if (vals) {
        const v = Math.max(-1, Math.min(1, vals[k] / sc));
        if (mag) { r = 255 * Math.min(1, v * 2); g = 255 * Math.max(0, v * 2 - 1) * 0.9; b = 40 * v; }
        else if (v > 0) { r = 12 + 243 * v; g = 14 + 80 * v; b = 20; } else { r = 12; g = 14 + 130 * -v; b = 20 + 235 * -v; }
      }
      const m = mm[k];
      if (m !== bgIdx) {
        const c = MAT_COLORS[(m % (MAT_COLORS.length - 1)) + 1];
        const a = view === 'material' ? 1 : 0.35;
        r = r * (1 - a) + c[0] * a; g = g * (1 - a) + c[1] * a; b = b * (1 - a) + c[2] * a;
      }
      d[p] = r; d[p + 1] = g; d[p + 2] = b; d[p + 3] = 255;
    }
    ctx.putImageData(img, 0, 0);
    const mark = (x: number, y: number, col: string, sz: number) => {
      ctx.fillStyle = col; ctx.fillRect(Math.round(x / s.dx) - sz / 2, ny - 1 - Math.round(y / s.dy) - sz / 2, sz, sz);
    };
    for (const p of s.cfg.ports) mark(p.x, p.y, '#ff3df0', 4);
    for (const pr of s.probes) for (const [x, y] of pr.coords) mark(x, y, '#f5f5f5', 2);
    for (const so of s.cfg.sources) mark(so.x, so.y, '#ffe600', 3);
  }, [view]);

  useEffect(() => { if (sim) draw(sim); }, [sim, view, draw, step]);

  const computeSpectrum = useCallback((s: Fdtd2D, pi: number) => {
    const p = s.ports[pi]; if (!p || s.step < 8) return setSpec(null);
    const o = s.cfg.output;
    setSpec(portSpectrum(p.v, p.i, s.step, s.dt, o.z0, o.fmin, o.fmax, o.fftPad ?? 4));
  }, []);

  useEffect(() => {
    if (!running || !sim) return;
    let raf = 0, total = 0, t0 = performance.now();
    const chunk = Math.max(1, sim.cfg.output.snapshotEvery);
    const loop = () => {
      const a = performance.now();
      const before = sim.step;
      const done = sim.advance(chunk);
      total += sim.step - before;
      setStep(sim.step);
      const el = (performance.now() - t0) / 1000;
      if (el > 0.5) { setPerf(`${((total * sim.nx * sim.ny) / el / 1e6).toFixed(1)} Mcells/s`); t0 = performance.now(); total = 0; }
      void a;
      if (done) { setRunning(false); computeSpectrum(sim, portIdx); return; }
      raf = requestAnimationFrame(loop);
    };
    raf = requestAnimationFrame(loop);
    return () => cancelAnimationFrame(raf);
  }, [running, sim, computeSpectrum, portIdx]);

  const runFast = () => {
    const s = sim ?? build(); if (!s) return;
    setRunning(false);
    const t = performance.now(), before = s.step;
    s.advance(s.steps);
    const el = (performance.now() - t) / 1000;
    setPerf(`${(((s.step - before) * s.nx * s.ny) / Math.max(el, 1e-6) / 1e6).toFixed(1)} Mcells/s (${el.toFixed(2)} s)`);
    setStep(s.step); computeSpectrum(s, portIdx);
  };

  const views = useMemo<View[]>(() => {
    if (!sim) return ['material'];
    return sim.mode === 'TMz' ? ['Ez', 'Hx', 'Hy', '|H|', 'material'] : ['Ex', 'Ey', 'Hz', '|E|', 'material'];
  }, [sim]);

  const chartData = useMemo(() => spec ? spec.f.map((f, k) => ({ f: f / 1e9, s11: spec.s11[k], zr: spec.zr[k], zi: spec.zi[k], vswr: Math.min(spec.vswr[k], 50) })) : [], [spec]);
  const timeData = useMemo(() => {
    const p = sim?.ports[portIdx]; if (!p || !sim) return [];
    const n = sim.step, st = Math.max(1, Math.floor(n / 500)), out = [];
    for (let k = 0; k < n; k += st) out.push({ t: ((k + 0.5) * sim.dt) / 1e-9, v: p.v[k], i: p.i[k] * 50, p: p.v[k] * p.i[k] });
    return out;
  }, [sim, portIdx, step]); // eslint-disable-line react-hooks/exhaustive-deps

  const btn = 'px-3 py-1.5 text-xs font-mono border border-[#2a3140] bg-[#161b24] text-[#d8dee9] hover:bg-[#222a38] disabled:opacity-40 rounded';

  return (
    <div className="min-h-screen bg-[#0d1016] text-[#d8dee9] font-mono text-xs">
      <header className="flex items-center gap-4 px-4 py-2 border-b border-[#1f2530]">
        <span className="text-[#f0a030] font-bold text-sm">TinyEM·2D</span>
        <span className="text-[#6b7688]">Yee FDTD · {sim ? `${sim.mode} ${sim.nx}×${sim.ny} dt=${sim.dt.toExponential(3)}s (CFL max ${sim.dtMax.toExponential(3)}) mem≈${(sim.memoryBytes / 1e6).toFixed(1)}MB` : 'no simulation'}</span>
      </header>
      <div className="grid grid-cols-1 xl:grid-cols-[360px_1fr_1fr] gap-3 p-3">
        <section className="flex flex-col gap-2">
          <div className="text-[#6b7688]">simulation.json</div>
          <textarea value={text} onChange={(e) => setText(e.target.value)} spellCheck={false}
            className="h-[560px] w-full bg-[#10141b] border border-[#1f2530] rounded p-2 text-[11px] leading-4 text-[#c8d0dc] outline-none focus:border-[#f0a030]" />
          <div className="flex flex-wrap gap-2">
            <button className={btn} onClick={build}>Build / Reset</button>
            <button className={btn} disabled={!sim} onClick={() => setRunning((r) => !r)}>{running ? 'Pause' : 'Run (animated)'}</button>
            <button className={btn} disabled={!sim} onClick={runFast}>Run to end (fast)</button>
          </div>
          <div className="flex flex-wrap gap-2">
            <button className={btn} disabled={!sim} onClick={() => sim && download('tinyem_results.json', exportJSON(sim), 'application/json')}>Export JSON</button>
            <button className={btn} disabled={!sim?.ports.length} onClick={() => sim && download('tinyem_ports.csv', exportPortCSV(sim), 'text/csv')}>Export port CSV</button>
            <button className={btn} onClick={() => download('simulation.json', text, 'application/json')}>Save config</button>
          </div>
          {error && <pre className="whitespace-pre-wrap text-[#ff6b6b]">{error}</pre>}
        </section>

        <section className="flex flex-col gap-2">
          <div className="flex flex-wrap gap-1">
            {views.map((v) => (
              <button key={v} onClick={() => setView(v)} className={`${btn} ${view === v ? '!bg-[#f0a030] !text-[#0d1016]' : ''}`}>{v}</button>
            ))}
          </div>
          <canvas ref={canvasRef} className="w-full aspect-square bg-black border border-[#1f2530]" style={{ imageRendering: 'pixelated' }} />
          <div className="h-1.5 bg-[#1f2530] rounded"><div className="h-full bg-[#f0a030] rounded" style={{ width: `${sim ? (100 * step) / sim.steps : 0}%` }} /></div>
          <div className="flex justify-between text-[#6b7688]">
            <span>step {step}/{sim?.steps ?? 0} · t={sim ? (step * sim.dt * 1e9).toFixed(3) : 0} ns</span>
            <span>{perf}</span>
          </div>
          <div className="text-[#6b7688]"><span className="text-[#ff3df0]">■</span> port <span className="text-white">■</span> probe <span className="text-[#ffe600]">■</span> source · materials tinted</div>
        </section>

        <section className="flex flex-col gap-2">
          <div className="flex items-center gap-2">
            <span className="text-[#6b7688]">Port</span>
            <select value={portIdx} onChange={(e) => { const v = +e.target.value; setPortIdx(v); if (sim) computeSpectrum(sim, v); }}
              className="bg-[#161b24] border border-[#2a3140] rounded px-2 py-1">
              {sim?.ports.map((p, k) => <option key={p.name} value={k}>{p.name}</option>)}
            </select>
            <button className={btn} disabled={!sim?.ports.length} onClick={() => sim && computeSpectrum(sim, portIdx)}>Compute FFT now</button>
            <span className="text-[#6b7688]">Z0={sim?.cfg.output.z0 ?? 50}Ω</span>
          </div>
          <Chart title="V(t) [V], I(t)·50 [V], P(t) [W]  vs t [ns]" data={timeData} x="t" lines={[['v', '#f0a030'], ['i', '#50aadc'], ['p', '#8fe39a']]} />
          <Chart title="S11 [dB] vs f [GHz]" data={chartData} x="f" lines={[['s11', '#f0a030']]} />
          <Chart title="Re(Z), Im(Z) [Ω] vs f [GHz]" data={chartData} x="f" lines={[['zr', '#50aadc'], ['zi', '#dc5096']]} />
          <Chart title="VSWR (clipped at 50) vs f [GHz]" data={chartData} x="f" lines={[['vswr', '#8fe39a']]} />
        </section>
      </div>
    </div>
  );
}

function Chart({ title, data, x, lines }: { title: string; data: object[]; x: string; lines: [string, string][] }) {
  return (
    <div className="border border-[#1f2530] rounded p-2 bg-[#10141b]">
      <div className="text-[#6b7688] mb-1">{title}</div>
      <div className="h-[150px]">
        {data.length ? (
          <ResponsiveContainer width="100%" height="100%">
            <LineChart data={data} margin={{ top: 4, right: 8, bottom: 0, left: 0 }}>
              <CartesianGrid stroke="#1f2530" />
              <XAxis dataKey={x} type="number" domain={['dataMin', 'dataMax']} tick={{ fill: '#6b7688', fontSize: 10 }} tickFormatter={(v: number) => v.toFixed(2)} />
              <YAxis tick={{ fill: '#6b7688', fontSize: 10 }} width={52} tickFormatter={(v: number) => (Math.abs(v) >= 1000 || (Math.abs(v) < 0.01 && v !== 0) ? v.toExponential(1) : v.toFixed(2))} />
              <Tooltip contentStyle={{ background: '#161b24', border: '1px solid #2a3140', fontSize: 11 }} />
              {lines.map(([k, c]) => <Line key={k} dataKey={k} stroke={c} dot={false} isAnimationActive={false} strokeWidth={1.2} />)}
            </LineChart>
          </ResponsiveContainer>
        ) : <div className="h-full flex items-center justify-center text-[#3d4656]">run simulation to populate</div>}
      </div>
    </div>
  );
}

