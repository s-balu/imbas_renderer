/*
 * compute_smoothing.c — adaptive smoothing length computation.
 *
 * Extracted from render_image.c.  Wraps the KERNEL_SMOOTHING ifdef
 * so the call site in main() is a single unconditional line.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "compute_smoothing.h"
#include "args.h"

#ifndef IMAGE_DIMENSIONX
#define IMAGE_DIMENSIONX 768
#endif

/* Forward declarations — defined in find_neighbours.c */
void find_neighbours_fast(int num_part, float *smoothing_length, int num_ngb,
                          float *posx, float *posy, float *posz,
                          float xmin, float ymin, float zmin,
                          float sph_eta, const char *cache_file,
                          int render_width);
void find_neighbours_cached(int num_part, float *smoothing_length, int num_ngb,
                             float *posx, float *posy, float *posz,
                             float xmin, float ymin, float zmin,
                             float sph_eta, const char *cache_file);

/* ------------------------------------------------------------------ */

float *compute_smoothing_lengths(const cli_args_t *cfg,
                                 float *x, float *y, float *z,
                                 long long NumPart,
                                 double xc, double yc, double zc)
{
#ifdef KERNEL_SMOOTHING
    /*
     * Adaptive smoothing lengths for SPH kernel deposit.
     * Not needed for CIC mode — smooth_to_mesh ignores smoothing_length
     * when KERNEL_SMOOTHING is not defined.
     */
    float *smoothing_length = (float *)malloc(sizeof(float) * NumPart);
    if (!smoothing_length) {
        fprintf(stderr, "compute_smoothing_lengths: malloc failed\n");
        return NULL;
    }

    /* Pass physical voxel size as xmin so find_neighbours can apply
     * a resolution-based h cap.  The y/z min args are unused. */
    float vox_size = (float)cfg->lbox / (float)IMAGE_DIMENSIONX;

    clock_t tstart = clock();

    if (cfg->fast_smooth) {
        find_neighbours_fast((int)NumPart, smoothing_length, cfg->num_ngb,
                             x, y, z,
                             vox_size, 0.0f, 0.0f, cfg->sph_eta,
                             cfg->sph_cache[0] ? cfg->sph_cache : NULL,
                             IMAGE_DIMENSIONX);
    } else {
        find_neighbours_cached((int)NumPart, smoothing_length, cfg->num_ngb,
                               x, y, z,
                               vox_size, 0.0f, 0.0f, cfg->sph_eta,
                               cfg->sph_cache[0] ? cfg->sph_cache : NULL);
    }

    clock_t tfinish = clock();
    fprintf(stdout, "find_neighbours - time: %.2f s\n",
            (double)(tfinish - tstart) / CLOCKS_PER_SEC);
    fflush(stdout);

    return smoothing_length;

#else
    (void)cfg; (void)x; (void)y; (void)z;
    (void)NumPart; (void)xc; (void)yc; (void)zc;
    fprintf(stdout, "CIC deposit mode (no kernel smoothing)\n");
    fflush(stdout);
    return NULL;   /* smooth_to_mesh ignores NULL in CIC mode */
#endif
}
