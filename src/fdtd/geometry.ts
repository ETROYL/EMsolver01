import type { Shape } from './types';

/**
 * Voxelises shapes onto the node map (node (i,j) at x=i*dx, y=j*dy). Later shapes overwrite earlier ones.
 */
export function rasterize(
  shapes: Shape[],
  matIndex: Map<string, number>,
  map: Uint8Array,
  nx: number,
  ny: number,
  dx: number,
  dy: number,
): void {
  const lx = (nx - 1) * dx + 1e-12;
  const ly = (ny - 1) * dy + 1e-12;
  shapes.forEach((s, si) => {
    const m = matIndex.get(s.material);
    if (m === undefined) throw new Error(`Invalid material "${s.material}" in geometry #${si}.`);
    if (s.type === 'cells') {
      for (const [i, j] of s.cells) {
        if (i < 0 || j < 0 || i >= nx || j >= ny) throw new Error(`Geometry exceeds computational domain (cell ${i},${j}).`);
        map[i * ny + j] = m;
      }
      return;
    }
    const [x0, x1, y0, y1] = bbox(s);
    if (x0 < -1e-12 || y0 < -1e-12 || x1 > lx || y1 > ly) throw new Error(`Geometry exceeds computational domain (shape #${si}, ${s.type}).`);
    const i0 = Math.max(0, Math.floor(x0 / dx) - 1), i1 = Math.min(nx - 1, Math.ceil(x1 / dx) + 1);
    const j0 = Math.max(0, Math.floor(y0 / dy) - 1), j1 = Math.min(ny - 1, Math.ceil(y1 / dy) + 1);
    const hw = s.type === 'line' ? Math.max((s.width ?? 0) / 2, 0.5 * Math.min(dx, dy)) + 1e-12 : 0;
    for (let i = i0; i <= i1; i++) {
      for (let j = j0; j <= j1; j++) {
        if (inside(s, i * dx, j * dy, hw)) map[i * ny + j] = m;
      }
    }
  });
}

function bbox(s: Exclude<Shape, { type: 'cells' }>): [number, number, number, number] {
  switch (s.type) {
    case 'box': return [s.xmin, s.xmax, s.ymin, s.ymax];
    case 'circle': return [s.cx - s.r, s.cx + s.r, s.cy - s.r, s.cy + s.r];
    case 'line': return [Math.min(s.x0, s.x1), Math.max(s.x0, s.x1), Math.min(s.y0, s.y1), Math.max(s.y0, s.y1)];
    case 'polygon': {
      const xs = s.points.map((p) => p[0]), ys = s.points.map((p) => p[1]);
      return [Math.min(...xs), Math.max(...xs), Math.min(...ys), Math.max(...ys)];
    }
  }
}

function inside(s: Exclude<Shape, { type: 'cells' }>, x: number, y: number, hw: number): boolean {
  const eps = 1e-12;
  switch (s.type) {
    case 'box': return x >= s.xmin - eps && x <= s.xmax + eps && y >= s.ymin - eps && y <= s.ymax + eps;
    case 'circle': return (x - s.cx) ** 2 + (y - s.cy) ** 2 <= s.r * s.r + eps;
    case 'line': {
      const vx = s.x1 - s.x0, vy = s.y1 - s.y0;
      const L2 = vx * vx + vy * vy;
      const t = L2 > 0 ? Math.max(0, Math.min(1, ((x - s.x0) * vx + (y - s.y0) * vy) / L2)) : 0;
      return Math.hypot(x - s.x0 - t * vx, y - s.y0 - t * vy) <= hw;
    }
    case 'polygon': {
      let c = false;
      const p = s.points;
      for (let a = 0, b = p.length - 1; a < p.length; b = a++) {
        if ((p[a][1] > y) !== (p[b][1] > y) && x < ((p[b][0] - p[a][0]) * (y - p[a][1])) / (p[b][1] - p[a][1]) + p[a][0]) c = !c;
      }
      return c;
    }
  }
}
 
