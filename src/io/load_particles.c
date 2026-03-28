/*
 * load_particles.c — snapshot reading and velocity interpolation.
 *
 * Extracted from render_image.c.  Contains:
 *   load_snapshot_header()  — filename resolution, header read, boxunits scaling
 *   load_particles()        — allocation, particle read, velocity interpolation
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "load_particles.h"
#include "args.h"
#include "io.h"
#include "globals.h"   /* ThisTask */

/* ------------------------------------------------------------------ */

int load_snapshot_header(cli_args_t      *cfg,
                         struct sim_info *header,
                         long long       *NumPart,
                         char            *filename_out,
                         int             *isDistributed_out)
{
    check_input_filenames(filename_out, cfg->file_root,
                          cfg->isHDF5, isDistributed_out);

    if (ThisTask == 0) {
        if (cfg->xcen > 0 || cfg->lbox > 0) {
            fprintf(stdout, "Centre: (%g|%g|%g)\n",
                    cfg->xcen, cfg->ycen, cfg->zcen);
            fprintf(stdout, "Box Length: %g\n", cfg->lbox);
        }
        fprintf(stdout, "Render configuration:\n");
        render_config_print(&cfg->rcfg);
        fprintf(stdout, "Reading header...\n");
        fflush(stdout);
    }

    if (cfg->isHDF5)
        read_hdf5_header(filename_out, header, NumPart);
    else
        read_gadget_binary_header(filename_out, header, NumPart);

    if (ThisTask == 0) {
        fprintf(stdout, "Number of particle types in file: %d\n",
                header->num_types);
        for (int i = 0; i < (int)(sizeof(int) * 8); i++) {
            if (cfg->ptype_mask != -1 && !(cfg->ptype_mask & (1 << i)))
                continue;
            if (i >= header->num_types || header->nall[i] == 0)
                fprintf(stdout,
                        "Warning: type %d requested but not present\n", i);
        }
        fflush(stdout);
    }

    /* Resolve BoxSize and apply box-unit scaling */
    if (header->BoxSize == 0)
        header->BoxSize = 1.e6;
    else if (cfg->boxunits == 1) {
        cfg->xcen *= header->BoxSize;
        cfg->ycen *= header->BoxSize;
        cfg->zcen *= header->BoxSize;
        cfg->lbox *= header->BoxSize;
    }

    if (ThisTask == 0) {
        fprintf(stdout, "Number of files: %d\n",       header->NumFiles);
        fprintf(stdout, "Number of particles: %lld\n", *NumPart);
        fprintf(stdout, "Time/Expansion Factor: %g\n", header->time);
        fflush(stdout);
    }

    return 0;
}

/* ------------------------------------------------------------------ */

int load_particles(const cli_args_t      *cfg,
                   const struct sim_info  *header,
                   const char             *filename,
                   float                **x,
                   float                **y,
                   float                **z,
                   int                  **ptype,
                   long long             *NumPartRead_out)
{
    long long NumPart = 0;

    /* Count total particles from header so we can allocate */
    for (int t = 0; t < header->num_types; t++)
        NumPart += header->nall[t];

    *x     = (float *)malloc(sizeof(float) * NumPart);
    *y     = (float *)malloc(sizeof(float) * NumPart);
    *z     = (float *)malloc(sizeof(float) * NumPart);
    *ptype = (int   *)malloc(sizeof(int)   * NumPart);

    if (!*x || !*y || !*z || !*ptype) {
        fprintf(stderr, "load_particles: malloc failed for position arrays\n");
        return -1;
    }

    int do_interp = (cfg->interp_frac != 0.0f && cfg->isHDF5);

    float *vx = NULL, *vy = NULL, *vz = NULL;
    if (do_interp) {
        vx = (float *)malloc(sizeof(float) * NumPart);
        vy = (float *)malloc(sizeof(float) * NumPart);
        vz = (float *)malloc(sizeof(float) * NumPart);
        if (!vx || !vy || !vz) {
            fprintf(stderr, "load_particles: malloc failed for velocity arrays\n");
            free(vx); free(vy); free(vz);
            return -1;
        }
    }

    if (ThisTask == 0) {
        fprintf(stdout, "Reading particles...\n");
        fflush(stdout);
    }

    if (cfg->isHDF5)
        read_particles_from_hdf5(cfg->file_root, *x, *y, *z,
                                  vx, vy, vz, *ptype,
                                  header->NumFiles, NumPartRead_out, do_interp);
    else
        read_particles_from_gadget_binary(cfg->file_root, *x, *y, *z, *ptype,
                                           header->NumFiles, NumPartRead_out);

    fprintf(stdout, "NumPart: %llu\tNumPartRead: %llu\n",
            NumPart, *NumPartRead_out);
    fflush(stdout);

    if (do_interp) {
        /*
         * Velocity interpolation: shift particle positions by frac * v * dt
         * to generate sub-snapshot frames between two outputs.
         *
         * This produces smooth transitions in a time sequence without needing
         * to store or read a second snapshot.  The approximation is first-order
         * (straight-line motion), which is accurate for small fractions of the
         * snapshot interval and breaks down near strong interactions.
         *
         * If -snap_dt is omitted, velocities are used as raw offsets (useful
         * when you just want to explore the velocity field visually).
         */
        float scale = (cfg->snap_dt > 0.0f)
                    ? cfg->interp_frac * cfg->snap_dt
                    : cfg->interp_frac;

        fprintf(stdout,
                "Interpolating positions: frac=%.3f dt=%g scale=%g\n",
                cfg->interp_frac, cfg->snap_dt, scale);
        fflush(stdout);

        for (long long k = 0; k < *NumPartRead_out; k++) {
            (*x)[k] += scale * vx[k];
            (*y)[k] += scale * vy[k];
            (*z)[k] += scale * vz[k];
        }

        free(vx); free(vy); free(vz);
    }

    return 0;
}
