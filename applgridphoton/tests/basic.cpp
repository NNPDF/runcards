#include <pineappl.h>

#include <catch2/catch.hpp>

#include <array>

void simple_pdf(double const& x, double const&, double* pdf)
{
    std::fill_n(pdf, 14, x);
}

double alphas(double const&)
{
    return 1.0;
}

TEST_CASE("", "")
{
    // create a new and empty luminosity function
    auto* lumi = pineappl_lumi_new();

    // add a new entry to the luminosity function; the combination `1.0 * (up up)`
    std::array<int, 2> pdg_id_pairs = { 2, 2 };
    std::array<double, 1> factors = { 1.0 };
    pineappl_lumi_add(lumi, pdg_id_pairs.size() / 2, pdg_id_pairs.data(), factors.data());

    // we'd like to have a single grid of order alpha^2 - for example LO Drell-Yan
    std::array<int, 4> grid_parameters = { 0, 2, 0, 0 };

    // global grid parameters, not really that important here
    int const    nq2      = 30;
    double const q2_min   = 100;
    double const q2_max   = 1000000;
    int const    q2_order = 1;//3;
    int const    nx       = 50;
    double const x_min    = 2e-7;
    double const x_max    = 1;
    int const    x_order  = 1;//3;

    // a distribution with only one bin
    std::array<double, 2> bin_limits = { 0.0, 1.0 };

    // create a new file
    auto* grid = pineappl_grid_new(
        lumi,
        pineappl_subgrid_format::as_a_logxir_logxif,
        grid_parameters.size() / 4,
        grid_parameters.data(),
        bin_limits.size() - 1,
        bin_limits.data(),
        nx,
        x_min,
        x_max,
        x_order,
        nq2,
        q2_min,
        q2_max,
        q2_order,
        "f2"
    );

    // delete luminosity function
    pineappl_lumi_delete(lumi);

    std::array<double, 1> const weights = { 1.0 };

    pineappl_grid_fill(grid, 0.25, 0.25, 10000.0, 0.25, weights.data(), 0);

    double result;

    pineappl_grid_convolute(
        grid,
        simple_pdf,
        simple_pdf,
        alphas,
        nullptr,
        1.0,
        1.0,
        1.0,
        &result
    );

    CHECK( result == 1.0 );

    pineappl_grid_delete(grid);
}
