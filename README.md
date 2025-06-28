# Spherical Harmonics

3D printable models of spherical harmonics.

The `.stl` files are available on [Thingiverse](https://www.thingiverse.com/thing:7077370) and in the [GitHub releases](https://github.com/djiboshin/spherical-harmonics/releases) (`SH_optimized.zip` file).

## Definition
The surface of each model is defined by real spherical harmonics:

$$\text S\text H_{lm+}(\theta, \phi) = \sqrt{2}(-1)^{m}\text{Re}\left[Y_{lm}(\theta, \phi)\right]$$
$$\text S\text H_{lm-}(\theta, \phi) = \sqrt{2}(-1)^{m}\text{Im}\left[Y_{lm}(\theta, \phi)\right]$$
$$\text S\text H_{l0+}(\theta, \phi) = \text{Re}\left[Y_{lm}(\theta, \phi)\right]$$

where $l=0,1,\dots$, $m=0,1,\dots,l$, $Y_{lm}$ is the complex-valued spherical harmonic as [defined in SciPy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sph_harm_y.html#scipy.special.sph_harm_y).
The $+$ and $-$ signs in the $\text{SH}$ index indicate the sign of $m$ in the interval $[-l, l]$.

A small central sphere has been added to each `.stl` model to soften sharp edges and provide additional support for floating "leaves", improving printability.

## Setup

Install all dependencies using [`uv`](https://github.com/astral-sh/uv):
```bash
uv sync
```

## Generation

Run the script to generate the models:
```bash
uv run generating.py
```
This script generates `.stl` files of spherical harmonics on a regular mesh and then optimizes them using `gmsh` library.

# License
Licensed under CC BY-NC-SA 4.0
