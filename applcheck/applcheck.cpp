#include <algorithm>
#include <array>
#include <cstddef>
#include <iostream>
#include <memory>
#include <vector>

#include <LHAPDF/LHAPDF.h>
#include <appl_grid/appl_grid.h>

std::unique_ptr<LHAPDF::PDF> pdf;

enum flavour_map_index : std::size_t
{
    anti_top,      // -6: anti-top
    anti_bottom,   // -5: anti-bottom
    anti_charm,    // -4: anti-charm
    anti_strange,  // -3: anti-strange
    anti_up,       // -2: anti-up
    anti_down,     // -1: anti-down
    gluon,         // 21: gluon
    down,          //  1: down
    up,            //  2: up
    strange,       //  3: strange
    charm,         //  4: charm
    bottom,        //  5: bottom
    top,           //  6: top
    photon,        // 22: photon
};

std::array<bool, 14> flavour_map = {
    true,  // -6: anti-top
    true,  // -5: anti-bottom
    true,  // -4: anti-charm
    true,  // -3: anti-strange
    true,  // -2: anti-up
    true,  // -1: anti-down
    true,  // 21: gluon
    true,  //  1: down
    true,  //  2: up
    true,  //  3: strange
    true,  //  4: charm
    true,  //  5: bottom
    true,  //  6: top
    true,  // 22: photon
};

constexpr int index_to_pdg_id(std::size_t index)
{
    return (index == gluon) ? 21 : ((index == photon) ? 22 : (static_cast <int> (index) - 6));
}

extern "C" void evolvepdf(double const& x, double const& q, double* xfx)
{
    for (std::size_t i = 0; i != flavour_map.size(); ++i)
    {
        if (flavour_map.at(i))
        {
            xfx[i] = pdf->xfxQ(index_to_pdg_id(i), x, q);
        }
        else
        {
            xfx[i] = 0.0;
        }
    }
}

extern "C" double alphaspdf(double const& q)
{
    return pdf->alphasQ(q);
}

int main(int argc, char* argv[])
{
    // set LHAPDF low verbosity
    LHAPDF::setVerbosity(0);

    if (argc != 3)
    {
        std::cerr << "Usage: applcheck PDF-set-name applgrid-file\n";
        return 1;
    }

    appl::grid g(argv[2]);
    const int nloops = g.nloops();

    // initialise the PDF set via LHAPDF6
    pdf.reset(LHAPDF::mkPDF(argv[1], 0));

    // check if PDF set has a photon; disable it this isn't the case
    flavour_map[photon] = pdf->hasFlavor(index_to_pdg_id(photon));

    std::vector<std::vector<double>> const& xsec_appl_orders
        = g.vconvolute_orders(evolvepdf, evolvepdf, alphaspdf);

    std::cout << "\n>>> all bins, all orders:\n\n";

    for (std::size_t i = 0; i != g.order_ids().size(); ++i)
    {
        auto const& xsecs = xsec_appl_orders.at(i);
        auto const& order = g.order_ids().at(i);

        for (std::size_t j = 0; j != g.Nobs_internal(); ++j)
        {
            if ((order.lmur2() != 0) || (order.lmuf2() != 0))
            {
                continue;
            }

            std::cout << " bin #" << std::setw(2) << j << ", O(as^" << order.alphs() << " a^"
                << order.alpha() << "): "<< std::scientific << std::setw(13) << xsecs.at(j)
                << " [pb(/GeV)]\n";
        }
    }

    std::vector<double> const& xsecs = g.vconvolute(evolvepdf, alphaspdf);

    std::cout << "\n>>> all bins:\n\n";

    double sum = 0.0;

    for (std::size_t i = 0; i != g.Nobs_internal(); ++i)
    {
        std::cout << " bin #" << std::setw(2) << i << ": " << std::scientific << std::setw(13)
            << xsecs.at(i) << " [pb(/GeV)]\n";

        sum += xsecs.at(i) * g.deltaobs_internal(i);
    }

    std::cout << "\n>>> sum:\n\n " << std::scientific << sum << " [pb]\n";

    pdf.release();
}
