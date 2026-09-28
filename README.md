# TinyEM 2D FDTD Solver

TinyEM is a browser-based 2D electromagnetic simulator implementing TMz and TEz finite-difference time-domain (FDTD) updates on a Yee grid.

> **Status:** research/prototype software. The numerical engine is implemented, but it is not yet validated against analytical/reference cases sufficiently to support quantitative antenna-design claims.

## Features

- TMz: Ez, Hx, Hy
- TEz: Ex, Ey, Hz
- CFL-limited time stepping
- Lossy isotropic materials with relative permittivity, permeability and conductivity
- PEC, PMC and first-order Mur outer boundaries
- Rasterised box, circle, polygon, line and cell geometry
- Gaussian, Ricker and CW excitations
- Explicit multi-edge feed gaps with finite Thévenin resistance
- Point and line probes
- FFT-based impedance, reflection coefficient, S11 and VSWR post-processing
- JSON/CSV result export

## Numerical model

The solver uses compact typed-array storage, while constitutive coefficients are evaluated at the corresponding staggered Yee locations. Electric material parameters are locally averaged from surrounding material samples; magnetic parameters are evaluated at the magnetic-field locations. Internal PEC regions force the corresponding tangential electric-field update coefficients to zero.

Finite-resistance ports are modelled as distributed Thévenin sources over the declared feed gap. For a multi-edge gap, the source voltage and resistance are distributed in series across the edges. Port voltage is obtained from the discrete line integral of the electric field and port current from the source circuit relation.

These are still grid-based approximations. Geometry is staircased/rasterised, and sub-cell conformal modelling is not implemented.

## Boundaries

- **Mur:** first-order absorbing boundary; useful for prototype work but not equivalent to a PML.
- **PEC:** tangential electric field is forced to zero at the selected outer boundary.
- **PMC:** the appropriate tangential-H / normal-E condition is explicitly enforced for the selected 2D mode.

A PML is **not** implemented yet.

## Port and S-parameter interpretation

The default example contains an actual PEC feed gap rather than placing the port inside a continuous PEC object.

The port resistance is the physical Thévenin source resistance. `output.z0` is the wave-reference impedance used when converting the computed port impedance to S11:

```
Z(f) = V(f) / I(f)
Gamma(f) = (Z(f) - Z0) / (Z(f) + Z0)
S11(dB) = 20 log10(|Gamma|)
```

Broadband results are only meaningful when the time record is long enough for the transient to decay, the feed is adequately resolved, and the spatial mesh is sufficiently fine for the highest frequency/material wavelength of interest.

## Time sampling

Electric and magnetic fields are staggered in time by half a timestep. Exported probe data therefore contains separate `t_e` and `t_h` axes. Port voltage/current histories are sampled consistently at the source half-step.

## Development

This repository is the **app project root**. Run development commands from this directory.

The repository contains a `pnpm-lock.yaml`, so pnpm is the reproducible package-manager path:

```bash
pnpm install
pnpm run dev
pnpm run build
pnpm run lint
```

npm can also execute the package scripts, but do not mix package managers in the same working tree unless you intentionally regenerate the lockfile.

For a local browser session, open:

```
http://localhost:3000/
```

## Project structure

```text
prerender/          Blog prerendering and sitemap helpers
public/             Static assets
src/fdtd/           FDTD engine, geometry, waveforms and post-processing
src/pages/          React application pages
src/components/     UI components
src/lib/             Shared configuration/blog utilities
vite.config.ts      Vite configuration
package.json        Dependencies and scripts
pnpm-lock.yaml      Reproducible pnpm dependency lockfile
```

## Validation roadmap

Before using TinyEM for quantitative antenna or GPR work, validate at minimum:

1. homogeneous-medium plane-wave propagation and numerical dispersion;
2. PEC reflection;
3. dielectric-interface Fresnel reflection/transmission;
4. Mur reflection versus distance/incidence angle;
5. feed-gap voltage/current against a known simple load;
6. impedance/S11 against an analytical or trusted reference;
7. mesh/time-step convergence.

The current implementation is deliberately kept 2D. PML, far-field transformation, directivity/gain, and general multiport S-parameter studies remain future work rather than silently being presented as implemented capabilities.
