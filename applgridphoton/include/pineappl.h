#ifndef PINEAPPL_H
#define PINEAPPL_H

/*
 * PineAPPL - PineAPPL Is Not an Extension of APPLgrid
 * Copyright (C) 2019-2020  Christopher Schwan
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifdef __cplusplus
extern "C"
{
#endif

/// @addtogroup lumi
/// @{

/// @struct pineappl_lumi
/// Struture that captures the definition of a luminosity function.
struct pineappl_lumi;

/// Creates a new luminosity function and returns a pointer to it. If no longer needed, the object
/// should be deleted using @ref pineappl_lumi_delete.
pineappl_lumi* pineappl_lumi_new();

/// Delete luminosity function previously created with @ref pineappl_lumi_new.
void pineappl_lumi_delete(pineappl_lumi* lumi);

/// Adds a linear combination of initial states to the luminosity function `lumi`.
void pineappl_lumi_add(
    pineappl_lumi* lumi,
    unsigned combinations,
    int* pdg_id_pairs,
    double* factors
);

/// @}

/// @addtogroup grid
/// @{

/// Enumeration that determines the meaning of the entries of the array `grid_parameters` in the
/// function @ref pineappl_grid_new. This is important only for @ref pineappl_grid_convolute, which
/// reconstructs the hadronic cross section.
typedef enum
{
    /// For each grid, this requires `grid_parameters` be a tuple of four integers \f$ (a, b, c, d)
    /// \f$, denoting the grid's 1) strong coupling power, \f$ a \f$, 2) the electromagnetic
    /// coupling power, \f$ b \f$, 3) the power of the renormalisation scale logarithm, \f$ c \f$,
    /// and finally, 4) the power of the factorisation scale logarithm, \f$ d \f$. For each bin,
    /// @ref pineappl_grid_convolute reconstructs the cross section as follows:
    /// \f[
    ///     \sigma = \left( \frac{\alpha_\mathrm{s}}{2 \pi} \right)^a
    ///              \left( \alpha \right)^b
    ///              \log^c \left( \xi_\mathrm{R} \right)
    ///              \log^d \left( \xi_\mathrm{F} \right) \sigma_{a,b,c,d}
    /// \f]
    as_a_logmur_logmuf = 1
}
pineappl_grid_format;

/// @struct pineappl_grid
/// Structure representing a PineAPPL grid.
struct pineappl_grid;

/// Create a new @ref pineappl_grid.
pineappl_grid* pineappl_grid_new(
    int n_bins,
    double const* bin_limits,
    pineappl_lumi* lumi,
    pineappl_grid_format format,
    int grids,
    int* grid_parameters,
    unsigned nq2,
    double q2_min,
    double q2_max,
    unsigned q2_order,
    unsigned nx,
    double x_min,
    double x_max,
    unsigned x_order,
    char const* map
);

/// Delete a grid previously created with @ref pineappl_grid_new.
void pineappl_grid_delete(pineappl_grid* grid);

/// Fill `grid` at the position specified with `x1`, `x2`, and `q2`. The array `weight` must be as
/// long the corresponding luminosity function the grid was created with and contain the
/// corresponding weights at each index. The value `grid_index` selects one of the subgrids whose
/// meaning was specified during creation with `grid_parameters` in @ref pineappl_grid_new.
void pineappl_grid_fill(
    pineappl_grid* grid,
    double x1,
    double x2,
    double q2,
    double observable,
    double const* weight,
    unsigned grid_index
);

/// Scale all grids in `grid` by `factor`.
void pineappl_grid_scale(pineappl_grid* grid, double factor);

/// Write `grid` to a file with name `filename`.
void pineappl_grid_write(pineappl_grid* grid, char const* filename);

// TODO: remove the `const&`

///
typedef void (*pineappl_func_xfx)(double const& x1, double const& q2, double* pdfs);

///
typedef double (*pineappl_func_alphas)(double const& q2);

/// Perform a convolution of the APPLgrid given as `grid` with the PDFs `pdf1` and `pdf2`. The
/// value `grid_mask` must either be `NULL` or an array as long as the luminosity function, which
/// can be used to selectively include (set the corresponding value in `grid_mask` to anything
/// non-zero) or exclude (zero) partonic initial states. The factors `scale_ren`, `scale_fac`, and
/// `scale_energy` can be used to scale the renormalisation scale, the factorisation scale, and the
/// energy from the their default (the value the grid was generated with) value. Results must be an
/// array as large as there are bins, and will contain the cross sections for each bin.
void pineappl_grid_convolute(
    pineappl_grid* grid,
    pineappl_func_xfx pdf1,
    pineappl_func_xfx pdf2,
    pineappl_func_alphas alphas,
    int* grid_mask,
    double scale_ren,
    double scale_fac,
    double scale_energy,
    double* results
);

/// @}

#ifdef __cplusplus
}
#endif

#endif
