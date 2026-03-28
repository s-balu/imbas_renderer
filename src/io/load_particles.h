#ifndef LOAD_PARTICLES_H
#define LOAD_PARTICLES_H

/*
 * load_particles.h — snapshot reading and velocity interpolation.
 *
 * load_snapshot_header() reads the file header and resolves boxunits scaling.
 * load_particles()       allocates position arrays, reads particles, applies
 *                        velocity interpolation if requested, then frees the
 *                        temporary velocity arrays.
 */

#include "args.h"    /* cli_args_t  */
#include "io.h"      /* sim_info    */

/*
 * load_snapshot_header()
 *
 * Resolves the input filename, reads the snapshot header, logs particle
 * counts, and applies box-unit scaling to cfg->xcen/ycen/zcen/lbox if
 * cfg->boxunits == 1.
 *
 * filename_out must be at least 256 bytes; receives the resolved path.
 * isDistributed_out is set to 1 if the snapshot is split across files.
 *
 * Returns 0 on success, -1 on error.
 */
int load_snapshot_header(cli_args_t       *cfg,
                         struct sim_info  *header,
                         long long        *NumPart,
                         char             *filename_out,
                         int              *isDistributed_out);

/*
 * load_particles()
 *
 * Allocates *x, *y, *z, *ptype (caller must free).
 * Reads positions (and velocities if interp_frac != 0 and isHDF5).
 * If interpolation is requested, shifts positions by frac*dt*v then
 * frees the velocity arrays internally.
 *
 * Returns the number of particles actually read in NumPartRead_out.
 * Returns 0 on success, -1 on allocation failure.
 */
int load_particles(const cli_args_t     *cfg,
                   const struct sim_info *header,
                   const char            *filename,
                   float               **x,
                   float               **y,
                   float               **z,
                   int                 **ptype,
                   long long            *NumPartRead_out);

#endif /* LOAD_PARTICLES_H */
