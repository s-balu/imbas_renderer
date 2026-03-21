# render_image

A high-performance N-body simulation volume renderer. Reads particle snapshots in Gadget binary or HDF5 format and produces PNG images with full control over colourmap, opacity, density scaling, and SPH kernel smoothing. Supports single-frame renders, zoom sequences, and rotation animations. Parallelised with OpenMP and optionally MPI.

---

## Features

- **CIC and SPH deposit modes** — cloud-in-cell for speed, or kernel-smoothed SPH projection for publication-quality images
- **O(N) smoothing length estimator** (`-fast_smooth`) — replaces ~160s exact KNN search with a 0.9s density-grid estimate; smoothing lengths are cached to disk
- **Direct 2D SPH projection** — projects particles straight onto the image plane using the analytically integrated cubic spline kernel, bypassing the 3D grid entirely
- **10 built-in colourmap palettes** — viridis, magma, inferno, plasma, hot, fire, ice, grayscale, coolwarm, custom
- **6 opacity transfer functions** — flat, linear, sqrt, power, log, threshold
- **Auto-levelling** — percentile-based vmin/vmax with optional frame-to-frame locking for animations
- **Scene presets** — cluster, scattered, filament
- **Animation support** — zoom sequences and rotation sequences with Rodrigues rotation
- **Multi-particle-type support** — gas, dark matter, stars, or any bitmask combination
- **OpenMP parallelism** — deposit and neighbour search scale across all available cores
- **MPI parallelism** — x-slab decomposition for distributed-memory clusters

---

## Dependencies

| Library | Purpose |
|---------|---------|
| libpng | PNG output |
| libhdf5 | HDF5 snapshot reading |
| libz | zlib compression (HDF5 dependency) |
| OpenMP | Shared-memory parallelism (optional) |
| MPI | Distributed-memory parallelism (optional) |

On macOS with Homebrew: `brew install gcc hdf5 libpng`

---

## Building

Edit the `COMPILE_ON_SYSTEM` variable at the top of `makefile` to match your platform, then:

```bash
make clean && make
```

Supported platform targets: `MacBook`, `MacPro`, `Magnus`, `OzSTAR`.

### Compile-time options

Set these in the `OPT` lines in the makefile or pass on the command line:

| Flag | Description |
|------|-------------|
| `-DCIC` | Cloud-in-cell deposit (default, fast) |
| `-DKERNEL_SMOOTHING` | SPH kernel deposit (slower, smoother) |
| `-DNONPERIODIC` | Non-periodic boundary conditions |
| `-DENABLE_OPENMP` | Enable OpenMP threading |
| `-DENABLE_MPI` | Enable MPI for distributed runs |
| `-DLONG_IDS` | 64-bit particle IDs |
| `-DMAX_H_VOXELS=N` | Max kernel radius in pixels at 768px (default 8) |
| `-DMAX_H_DEPOSIT_PX=N` | Hard pixel cap applied during deposit (default 8) |
| `-DMIN_RHO_FRAC=f` | Skip particles below this fraction of mean density (default 0.1) |
| `-DIMAGE_DIMENSIONX=N` | Output image width in pixels (default 768) |
| `-DIMAGE_DIMENSIONY=N` | Output image height in pixels (default 768) |

---

## Usage

```
./render_image.exe -input <snapshot> -output <prefix> [options]
```

### Required

| Flag | Description |
|------|-------------|
| `-input <file>` | Snapshot root path (without `.hdf5` extension) |
| `-output <file>` | Output image filename root |

### Particle selection

| Flag | Description |
|------|-------------|
| `-isHDF5` | HDF5 input format (default: Gadget binary) |
| `-dark_matter` | Render type 1 particles (default) |
| `-gas` | Render type 0 particles |
| `-stars` | Render type 4 particles |
| `-all_types` | Render all particle types |
| `-ptype <N>` | Keep particle type N (repeatable) |

### View volume

| Flag | Description |
|------|-------------|
| `-xc <val>` | X centre of view volume |
| `-yc <val>` | Y centre of view volume |
| `-zc <val>` | Z centre of view volume |
| `-lbox <val>` | Side length of view volume (same units as snapshot) |

### Colour and opacity

| Flag | Description |
|------|-------------|
| `-colormap <name>` | Palette: `viridis` `magma` `inferno` `plasma` `hot` `fire` `ice` `grayscale` `coolwarm` `custom` |
| `-reverse_colormap` | Reverse the palette |
| `-bg_color R,G,B[,A]` | Background colour in [0,1] (default `0,0,0,1`) |
| `-opacity <val>` | Global opacity multiplier 0–1 (default 1) |
| `-opacity_func <name>` | Transfer function: `flat` `linear` `sqrt` `power` `log` `threshold` |
| `-opacity_gamma <val>` | Exponent for `power` function |
| `-opacity_threshold <val>` | Cutoff for `threshold` function |

### Density scaling

| Flag | Description |
|------|-------------|
| `-auto_pct_lo <val>` | Low clip percentile 0–1 (default 0.001) |
| `-auto_pct_hi <val>` | High clip percentile 0–1 (default 0.999) |
| `-no_auto_levels` | Use fixed vmin/vmax instead of auto-levelling |
| `-vmin <val>` | Lower density clip in log10 units (with `-no_auto_levels`) |
| `-vmax <val>` | Upper density clip in log10 units (with `-no_auto_levels`) |
| `-linear_scale` | Use linear density scaling instead of log10 |

### Scene presets

| Flag | Description |
|------|-------------|
| `-scene cluster` | Magma colormap, sqrt opacity, 5% low clip — good for galaxy clusters |
| `-scene scattered` | Plasma colormap, flat opacity, 10% low clip — good for sparse distributions |
| `-scene filament` | Inferno colormap, power γ=0.5 opacity, 2% low clip — good for cosmic web |

### Animation

| Flag | Description |
|------|-------------|
| `-itmax <N>` | Render N frames with fixed view volume |
| `-zoom <N>` | Render N frames, shrinking the box each frame |
| `-zoom_factor <f>` | Multiply box side by f per zoom frame, 0 < f < 1 (default 0.5) |
| `-rot_dangle <deg>` | Rotate view by this many degrees per frame |
| `-rot_axis x,y,z` | Rotation axis (default `0,0,1`) |
| `-lock_levels` | Auto-level on frame 0, lock vmin/vmax for all subsequent frames |

Output frames are named `<prefix>.NNNN.png`.

### SPH kernel smoothing (requires `-DKERNEL_SMOOTHING`)

| Flag | Description |
|------|-------------|
| `-fast_smooth` | O(N) grid-based h estimator (~0.9s) instead of exact KNN (~160s) |
| `-sph_cache <file>` | Cache smoothing lengths; reloaded if particle count and parameters match |
| `-num_ngb <N>` | Target neighbour count for smoothing length (default 32) |
| `-sph_eta <val>` | h = eta × dist_to_Nth_neighbour (default 1.2; larger = smoother) |
| `-ngrid_z <N>` | Depth of 3D density grid for CIC mode (default = image width) |

---

## Examples

### Single cluster image

```bash
./render_image.exe \
    -input snapshot_122 -output cluster \
    -isHDF5 -xc 56.18 -yc 43.35 -zc 49.13 -lbox 0.2 \
    -dark_matter -scene cluster \
    -fast_smooth -sph_cache snap122.hcache -itmax 1
```

### Large volume render

```bash
./render_image.exe \
    -input snapshot_122 -output volume \
    -isHDF5 -xc 56.18 -yc 43.35 -zc 49.13 -lbox 1.0 \
    -dark_matter -colormap plasma \
    -auto_pct_lo 0.1 -auto_pct_hi 0.9999 -bg_color 0,0,0,1 \
    -fast_smooth -sph_cache snap122.hcache -itmax 1
```

### 36-frame rotation animation

```bash
./render_image.exe \
    -input snapshot_122 -output frames \
    -isHDF5 -xc 56.18 -yc 43.35 -zc 49.13 -lbox 0.2 \
    -dark_matter -scene cluster \
    -fast_smooth -sph_cache snap122.hcache \
    -itmax 36 -rot_dangle 10 -rot_axis 0,1,0 -lock_levels

ffmpeg -framerate 24 -i frames.%04d.png -c:v libx264 -pix_fmt yuv420p rotation.mp4
```

### Zoom sequence

```bash
./render_image.exe \
    -input snapshot_122 -output zoom \
    -isHDF5 -xc 56.18 -yc 43.35 -zc 49.13 -lbox 2.0 \
    -dark_matter -scene cluster \
    -fast_smooth -sph_cache snap122.hcache \
    -zoom 20 -zoom_factor 0.85 -lock_levels
```

---

## Source file overview

| File | Purpose |
|------|---------|
| `render_image.c` | Main entry point, CLI parsing, render loop, rotation/zoom logic |
| `io.c` | HDF5 and Gadget binary snapshot reader |
| `find_neighbours.c` | Smoothing length computation — O(N) density-grid estimator and exact KNN; disk caching |
| `flat_kd_tree.h` | Header-only flat-array kd-tree for exact KNN search |
| `deposit_sph_2d.c` | Direct 2D SPH kernel projection (separate TU for compiler inlining) |
| `smooth_to_mesh.c` | Density field construction — SPH 2D deposit or CIC 3D deposit + projection |
| `kernels.c` | Cubic spline kernel and 2D projected kernel; lookup tables |
| `colormap.c` | Colour palettes, opacity transfer functions, `render_config_t` pipeline |
| `colormap.h` | `render_config_t` struct definition and palette/function enums |
| `select_particles.c` | View-volume selection and particle type masking |
| `make_tree.c` | Octree construction for spatial indexing |
| `walk_tree.c` | Octree traversal functions |
| `split_across_tasks.c` | MPI x-slab decomposition |
| `header.c` / `header.h` | Global state, `sim_info` struct, shared declarations |
| `write_to_ppm.c` | PNG output with colormap, opacity, and background compositing |

---

## Performance notes

Typical timings on a 4-core Apple M2 for 14M dark matter particles at 768×768:

| Step | Time |
|------|------|
| HDF5 read + particle selection | ~4s |
| `-fast_smooth` smoothing lengths | ~0.9s |
| Smoothing length cache load | <0.1s |
| SPH 2D deposit | ~15–30s |
| PNG write | <0.1s |

The SPH deposit is memory-bandwidth limited (stride-768 write pattern). The `-fast_smooth` flag is strongly recommended; the exact KNN alternative takes ~160s for the same dataset. Smoothing lengths are stable between frames of a rotation animation so the cache eliminates the cost entirely after the first frame.

For large volumes or exploratory work, omit `-DKERNEL_SMOOTHING` to use CIC deposit (~3s).

---

## Smoothing length tuning (`-fast_smooth`)

The O(N) estimator bins particles onto a 128³ density grid, applies a single box-smooth pass, then estimates h from the local density via trilinear interpolation. Particles below `MIN_RHO_FRAC` (default 0.1) times the mean density are skipped — they contribute negligible signal and avoiding them eliminates isolated-dot artefacts in void regions. Skipped pixels are filled by the zero-hole inpainting pass in the render loop.

To adjust aggressiveness:
- `-DMIN_RHO_FRAC=0.05` — keep more sparse particles (more dots visible in voids)
- `-DMIN_RHO_FRAC=0.2` — skip more sparse particles (cleaner background)

---

## Custom colourmap

Pass a text file with one `R G B` triplet per line (values 0–255, 256 lines) via `-colormap_file <path>` combined with `-colormap custom`.
