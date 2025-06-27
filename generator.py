# %% [markdown]
# # SH .stl generation

# %% [markdown]
# Packages import

# %%
import numpy as np
import scipy.special as spf
from stl import mesh  # numpy-stl
from pathlib import Path

# Set folders where raw and optimized `.stl` files will be stored.
SH_DIR = Path("SH")
SH_OPT_DIR = Path("SH_optimized")

# %% [markdown]
# Define functions for `.stl` generation.
# `spfunc()` can be redefined in order to draw any $r(\theta, \phi)$ surface.


# %%
def spfunc(phi, theta, m, l, r, r_smooth: float = 1):
    """Points on real spherical harmonic (SH) surface.
    m, l are integers, l>=0, l>=|m|.
    r is a radius of little supporting sphere in a middle of SH."""
    # Real SH [doi.org/10.1016/S0166-1280(97)00185-1]
    Phi, Theta = np.meshgrid(phi, theta)

    if m > 0:
        HARMONIC = np.sqrt(2) * (-1) ** m * spf.sph_harm_y(l, abs(m), Theta, Phi).real
    elif m < 0:
        HARMONIC = np.sqrt(2) * (-1) ** m * spf.sph_harm_y(l, abs(m), Theta, Phi).imag
    else:
        HARMONIC = spf.sph_harm_y(l, abs(m), Theta, Phi).real

    ABS = np.abs(HARMONIC) * 100  # useful scale
    # add little sphere inside SH

    # soft-max for smoothing transition between SH and little sphere
    r_smooth = 1
    R = r + 1 / r_smooth * np.log(1 + np.exp(r_smooth * (ABS - r)))

    X = R * np.sin(Theta) * np.cos(Phi)
    Y = R * np.sin(Theta) * np.sin(Phi)
    Z = R * np.cos(Theta)

    return np.stack((X, Y, Z), axis=2).reshape(Phi.shape[0] * Phi.shape[1], 3)


def linrange(start, step, num):
    return np.linspace(start, start + step * num, num, endpoint=False)


def faces_generator(phi, theta):
    """
    Generate mesh faces for STL from phi/theta linspaces.
    Returns faces: (N, 3) numpy array of triangle indices for STL mesh
    """
    n_theta = len(theta)
    n_phi = len(phi)
    faces = []

    # Side faces (between latitude bands)
    for t in range(n_theta - 1):
        for p in range(n_phi):
            p_next = (p + 1) % n_phi
            i0 = t * n_phi + p
            i1 = t * n_phi + p_next
            i2 = (t + 1) * n_phi + p
            i3 = (t + 1) * n_phi + p_next

            # Each quad is split into two triangles
            faces.append([i0, i2, i1])
            faces.append([i1, i2, i3])

    # Top pole faces
    for p in range(n_phi):
        p_next = (p + 1) % n_phi
        faces.append([0, p, p_next])

    # Bottom pole faces
    base = (n_theta - 1) * n_phi
    for p in range(n_phi):
        p_next = (p + 1) % n_phi
        faces.append([base + p, base, base + p_next])

    return np.array(faces, dtype=int)


# %% [markdown]
# Here the `.stl` files are generating with `numpy-stl` package.

# %%
# create folder
SH_DIR.mkdir(exist_ok=True, parents=True)

# generate .stl for spherical harmonics
max_l = 3
# little sphere radius
r = 7  # size in mm
# overall scale
SCALE = 1.5
# number of phi and theta samples
n_phi, n_theta = 500, 500


phi = np.linspace(0, 2 * np.pi, n_phi, endpoint=False)
theta = np.linspace(0, np.pi, n_theta)

for l in range(max_l + 1):
    for m in range(-l, l + 1):
        points = spfunc(phi, theta, m=m, l=l, r=r)
        faces = faces_generator(phi, theta)

        my_mesh: mesh.Mesh = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
        my_mesh.vectors = points[faces, :] * 0.5 * SCALE  # 0.5 somehow works
        my_mesh.save(SH_DIR / f'{l}{abs(m)}{"-" if m < 0 else "+"}.stl')

# %% [markdown]
# # .stl optimization
#
# Define function for mesh optimisation with the `gmsh` package.

# %%
import gmsh
import math
import tqdm


def optimize(stl_path: str, result_path: str, target_size):
    gmsh.clear()
    gmsh.merge(stl_path)

    angle = 40
    forceParametrizablePatches = True
    includeBoundary = True
    curveAngle = 180

    gmsh.model.mesh.classifySurfaces(
        angle * math.pi / 180.0,
        includeBoundary,
        forceParametrizablePatches,
        curveAngle * math.pi / 180.0,
    )

    gmsh.model.mesh.createGeometry()

    s = gmsh.model.getEntities(2)
    l = gmsh.model.geo.addSurfaceLoop([e[1] for e in s])

    gmsh.model.geo.addVolume([l])
    gmsh.model.geo.synchronize()

    # Set mesh size to reduce number of faces
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", target_size)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 2 * target_size)

    # Smoothing, optimizing
    gmsh.option.setNumber("Mesh.Smoothing", 5)
    gmsh.option.setNumber("Mesh.Optimize", 1)  # Basic optimization
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)  # Use Netgen optimizer

    # Generate mesh
    gmsh.model.mesh.generate(2)

    # Save as STL (surface mesh)
    gmsh.write(result_path)


# %% [markdown]
# Optimization of all `.stl` files in `SH_DIR` folder.

# %%
# create folder
SH_OPT_DIR.mkdir(exist_ok=True, parents=True)

for file in tqdm.tqdm(sorted(SH_DIR.iterdir())):
    if not file.is_file():
        continue
    try:
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)

        optimize(
            stl_path=str(file), result_path=str(SH_OPT_DIR / file.name), target_size=0.4
        )
    finally:
        gmsh.finalize()

# %%
