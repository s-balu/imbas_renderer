#ifndef RENDER_IMAGE_HEADER_H
#define RENDER_IMAGE_HEADER_H

#include <stdint.h>
#include <png.h>
#include "colormap.h"

#define MAXNODE 8
#define NDIM 3

/* Maximum neighbours returned in a flat neighbour buffer.
   Increase if you have very high local particle counts. */
#define MAX_NGP 200000

/* Number of children per node: 2^NDIM */
#define NCHILDREN (1 << NDIM)

/* Upper bound on particle types used only to size the fixed Gadget-binary
   header block (256 bytes). The actual number of types present in a file
   is stored in sim_info.num_types and used for all runtime loops. */
#define GADGET_MAX_TYPES 6

struct link_list {
  float x[NDIM];
  struct link_list *ptr;
};

struct interaction_list {
  int index;
  struct tree_node_nd *node;
  struct interaction_list *left, *right;
};

struct tree_node {
  float x[MAXNODE];
  float xmin, xmax;
  int num_members;
  int i;
  int split;
  struct tree_node *left, *right;
};

/* N-dimensional tree node.
   children[k] replaces the named p0..p7 pointers.
   k is a bitmask: bit i set => coordinate i is in the upper half-cell. */
struct tree_node_nd {
  float x[NDIM][MAXNODE];
  float xmin[NDIM], xmax[NDIM];
  float xmean[NDIM];
  int num_members;
  int index;
  int split;
  struct tree_node_nd *children[NCHILDREN];
};

struct btree_node {
  float x;
  struct btree_node *left, *right;
};

struct point {
  float pos[NDIM];
  float mass;
  float density;
  float dist_ngb;
};

/* Flat neighbour buffer: replaces malloc'd linked lists in hot paths.
 *
 * Stored on the heap with capacity set at runtime from num_ngb, so the
 * buffer fits in L2 cache rather than being a fixed 3 MB struct.
 * Allocated once per thread by ngb_buf_alloc() and freed by ngb_buf_free().
 *
 * SoA layout (separate x0/x1/x2 arrays) allows the distance loop to be
 * auto-vectorised.  r2[] is filled during the tree walk so no second
 * distance pass is needed.
 */
typedef struct {
  float  *x0;    /* x coordinates of candidate neighbours */
  float  *x1;    /* y coordinates                         */
#if NDIM == 3
  float  *x2;    /* z coordinates                         */
#endif
  float  *r2;    /* squared distances (computed in-place)  */
  int     count;
  int     capacity;
} ngb_buf_t;

/* Allocate a buffer sized for the expected number of candidates.
 * cap = 4 * num_ngb gives a comfortable margin with the sphere filter
 * (which limits actual fill to ~2*num_ngb) and stays well under 100 KB. */
ngb_buf_t *ngb_buf_alloc(int num_ngb);
void       ngb_buf_free(ngb_buf_t *buf);

/* ------------------------------------------------------------------ */
/* Simulation metadata                                                  */
/* ------------------------------------------------------------------ */

/*
 * sim_info holds everything read from the snapshot header.
 *
 * num_types   – number of particle types actually present in the file.
 *               Discovered at header-read time.  All runtime loops over
 *               types use this value instead of the hard-coded 6.
 *
 * npart / nall / massarr are dynamically allocated arrays of length
 * num_types so the code works with any format that reports a different
 * number of species (e.g. SWIFT uses 7, some custom codes use 2).
 *
 * The Gadget binary header block is still 256 bytes and stores exactly
 * GADGET_MAX_TYPES (6) entries; the reader zero-pads or truncates as
 * needed and then sets num_types from the count of non-zero entries.
 */
struct sim_info {
  int    num_types;        /* actual number of types in this snapshot   */
  int   *npart;            /* [num_types] particles in this file        */
  int   *nall;             /* [num_types] particles across all files    */
  double *massarr;         /* [num_types] table masses (0 = per-part)   */
  double time;
  double redshift;
  int    NumFiles;
  double BoxSize;
  double Omega0, OmegaLambda, HubbleParam;
  int    SnapFormat;
};

/* Allocate / free the dynamic arrays inside a sim_info */
void sim_info_alloc(struct sim_info *h, int num_types);
void sim_info_free(struct sim_info *h);

/* ------------------------------------------------------------------ */
/* Tree building                                                        */
/* ------------------------------------------------------------------ */
struct tree_node    *add_to_node(struct tree_node *, float, float, float);
struct tree_node_nd *add_to_node_nd(struct tree_node_nd *, float *, float *, float *);

void make_tree(struct point *, int, struct tree_node_nd **);

/* ------------------------------------------------------------------ */
/* Tree walking                                                         */
/* ------------------------------------------------------------------ */


/* Flat-buffer neighbour search (no malloc per neighbour) */
void get_multiple_nodes_nd_flat(struct tree_node_nd *, const float *, float, ngb_buf_t *);




/* ------------------------------------------------------------------ */
/* High-level algorithms                                                */
/* ------------------------------------------------------------------ */
void get_points(struct point *, int *, int);
void get_distance_to_nth_nearest_neighbour(struct point *, int, int, struct tree_node_nd *);
void get_kernel_density_estimate(struct point *, int, struct tree_node_nd *);
void get_potential_estimate(struct point *, int, struct tree_node_nd *);

int cmpfunc(const void *, const void *);
float get_rand(int);

/* ------------------------------------------------------------------ */
/* Kernels                                                              */
/* ------------------------------------------------------------------ */
void tabulate_kernel(void);
#ifdef KERNEL_SMOOTHING
void tabulate_kernel_u2(void);
void tabulate_projected_kernel(void);
float cubic_spline_kernel(float);
float cubic_spline_kernel_2d_proj(float);
float get_kernel_value(float);
float get_kernel_value_u2(float);   /* argument is (r/h)^2, no sqrtf needed */
#endif /* KERNEL_SMOOTHING */
float get_projected_kernel_value(float);
void  tabulate_projected_kernel_u2(void);
float get_projected_kernel_value_u2(float);  /* arg is (R/h)^2, no sqrtf */
float cumulative_cubic_spline_interpolant(float);
void tabulate_integral(void);
float get_integral_value(float);

/* ------------------------------------------------------------------ */
/* I/O                                                                  */
/* ------------------------------------------------------------------ */
void check_input_filenames(char *, char *, int, int *);
void read_hdf5_header(char *, struct sim_info *, long long *);
void read_gadget_binary_header(char *, struct sim_info *, long long *);
void read_particles_from_hdf5(char *, float *, float *, float *, int *, int, long long *);
void read_particles_from_gadget_binary(char *, float *, float *, float *, int *, int, long long *);

/*
 * select_particles: keep particles inside the view volume.
 *
 * ptype_mask – bitmask of particle types to keep.
 *              Bit i set => keep type i.
 *              Pass -1 (all bits set) to keep every type.
 *              This replaces the old single-integer ptype_keep, so any
 *              number of types can be combined without changing the
 *              function signature again.
 */
void select_particles(float *, float *, float *, int *, double, long long,
                      float, float, float, float, int ptype_mask, long long *);

void split_across_tasks_as_slabs(float *, float *, float *, long long *,
                                  float, float, float, float, float);
void smooth_to_mesh(long long, float *, float *, float *, float *,
                    float, float, float, float, float, int, int, float *);

/* ------------------------------------------------------------------ */
/* Image output                                                         */
/* ------------------------------------------------------------------ */
void write_to_ppm(char *, int, int, int, float *);

typedef struct { uint8_t red; uint8_t green; uint8_t blue; } pixel_t;
typedef struct { pixel_t *pixels; size_t width; size_t height; }  bitmap_t;

void write_to_png(char *, int, int, float *);
void write_to_png_ex(const char *, int, int, const float *,
                     const render_config_t *);
pixel_t *pixel_at(bitmap_t *, int, int);
int save_png_to_file(bitmap_t *, const char *);


/* ------------------------------------------------------------------ */
/* Direct 2D SPH deposit (deposit_sph_2d.c)                            */
/* ------------------------------------------------------------------ */
#ifdef KERNEL_SMOOTHING
typedef struct { float px; float py; float h; } deposit_particle_t;
void       deposit_sph_2d_init(void);
long long  deposit_sph_2d(long long n_packed,
                           const deposit_particle_t *particles,
                           int full_width, int height,
                           float *data2d);
#endif

/* ------------------------------------------------------------------ */
/* Globals                                                              */
/* ------------------------------------------------------------------ */
extern int ThisTask, NTask;

/* Per-task x-slab boundaries in physical coordinates.
 * slab_x_lo[t] and slab_x_hi[t] are the physical x-range owned by task t.
 * Set by split_across_tasks_as_slabs() and read by smooth_to_mesh().
 * Allocated as NTask+1 floats (slab_x_lo[t]=slab_x_hi[t-1]).
 */
extern float *slab_x_lo;   /* [NTask]: lower x boundary per task */
extern float *slab_x_hi;   /* [NTask]: upper x boundary per task */
extern int SnapFormat;
extern int NThisTask;

#endif /* RENDER_IMAGE_HEADER_H */
