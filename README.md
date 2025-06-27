# spherical-harmonics

3D printable models of spherical harmonics.

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
