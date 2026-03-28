#ifndef COMPUTE_SMOOTHING_H
#define COMPUTE_SMOOTHING_H

/*
 * compute_smoothing.h — adaptive smoothing length computation.
 *
 * Wraps the KERNEL_SMOOTHING ifdef so the caller doesn't need to
 * repeat it.  In CIC mode (KERNEL_SMOOTHING not defined) the function
 * returns NULL and logs a message; the caller passes NULL straight
 * through to smooth_to_mesh(), which ignores it in that mode.
 */

#include "args.h"   /* cli_args_t */

/*
 * compute_smoothing_lengths()
 *
 * Allocates and fills a smoothing-length array for every particle.
 * Uses the fast O(N) grid estimator if cfg->fast_smooth is set,
 * otherwise uses the cached exact KNN search.
 *
 * Returns a malloc'd float array of length NumPart that the caller
 * must free, or NULL if compiled without -DKERNEL_SMOOTHING.
 */
float *compute_smoothing_lengths(const cli_args_t *cfg,
                                 float *x, float *y, float *z,
                                 long long NumPart,
                                 double xc, double yc, double zc);

#endif /* COMPUTE_SMOOTHING_H */
