#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


#include <geogram/basic/common.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/stopwatch.h>
#include <geogram/basic/file_system.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_reorder.h>
#include <geogram/delaunay/delaunay.h>
//#include <geogram/delaunay/delaunay_3d.h>
#include <geogram/delaunay/parallel_delaunay_3d.h>


#include <algorithm>

#include <iostream>

using namespace GEO;


#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)


// // type de retour: std::tuple<std::vector<double>, std::vector<int>>
// // (points, tetrahedra) = delauynet_fast(points)
// // std::tuple
// Delaunay_var get_delaunay_neighbors(std::vector<double> points) {
//     GEO::initialize(GEO::GEOGRAM_INSTALL_ALL);
//     const index_t dimension = 3;   

//     //index_t nb_points = points.size() / dimension;
//     index_t nb_points = points.size();

//     Delaunay_var delaunay = Delaunay::create(dimension,"PDEL");
//     delaunay->set_stores_neighbors(true);
//     delaunay->set_vertices(nb_points, points.data());
//     delaunay->update_neighbors();


//     return delaunay;
// }


// std::vector<std::vector<index_t>> get_delaunay_neighbors(std::vector<double> points) {
//     GEO::initialize(GEO::GEOGRAM_INSTALL_ALL);
//     const index_t dimension = 3;
//     index_t nb_points = points.size(); // / dimension;

//     Delaunay_var delaunay = Delaunay::create(dimension, "PDEL");
//     delaunay->set_stores_neighbors(true);
//     delaunay->set_vertices(nb_points, points.data());

//     GEO::vector<GEO::vector<GEO::index_t>> all_neighbors(nb_points);

//     for(GEO::index_t i = 0; i < nb_points; ++i) {
//         delaunay->get_neighbors(i, all_neighbors[i]); // OK now!
//     }

//     return all_neighbors;
// }


// std::vector<std::vector<GEO::index_t>> get_delaunay_neighbors(std::vector<double> points) {
//     GEO::initialize(GEO::GEOGRAM_INSTALL_ALL);
//     const GEO::index_t dimension = 3;
//     GEO::index_t nb_points = points.size() / dimension;

//     Delaunay_var delaunay = Delaunay::create(dimension, "PDEL");
//     //delaunay->set_stores_neighbors(true);
//     delaunay->set_vertices(nb_points, points.data());

//     GEO::vector<GEO::vector<GEO::index_t>> geo_neighbors(nb_points);

//     // for(GEO::index_t i = 0; i < nb_points; ++i) {
//     //     delaunay->get_neighbors(i, geo_neighbors[i]);
//     // }


//     // // Convert to std::vector for return
//     // std::vector<std::vector<GEO::index_t>> std_neighbors;
//     // std_neighbors.reserve(nb_points);

//     // for(const auto& gvec : geo_neighbors) {
//     //     std_neighbors.emplace_back(gvec.begin(), gvec.end());
//     // }

    
//     simplices = delaunay->cell_to_v();


//     return simplices;
// }

std::vector<std::vector<GEO::index_t>> get_delaunay_simplices(
    const std::vector<double>& points
) {
    GEO::initialize(GEO::GEOGRAM_INSTALL_ALL);
    GEO::index_t dimension = 3; // Assuming 3D points
    GEO::index_t nb_points = points.size() / dimension;

    Delaunay_var delaunay = Delaunay::create(dimension, "PDEL");
    //std::cout << "Delaunay implementation: " << typeid(*delaunay.get()).name() << std::endl;

    
    delaunay->set_vertices(nb_points, points.data());


    // Get the reorder map
    // auto* delaunay3d = dynamic_cast<GEO::ParallelDelaunay3d*>(delaunay.get());
    
    // geo_assert(delaunay3d != nullptr);
    // const GEO::vector<GEO::index_t>& reorder = delaunay3d->get_reorder();

    // // print reorder
    // std::cout << "Reorder: ";
    // for (GEO::index_t i = 0; i < reorder.size(); ++i) {
    //     std::cout << reorder[i] << " ";
    // }
    // std::cout << std::endl;


    // Extract simplices
    GEO::index_t num_cells = delaunay->nb_cells();
    GEO::index_t verts_per_cell = delaunay->cell_size();
    const GEO::signed_index_t* simplices_ptr = delaunay->cell_to_v();

    std::vector<std::vector<GEO::index_t>> simplices(num_cells);
    for (GEO::index_t i = 0; i < num_cells; ++i) {
        simplices[i].reserve(verts_per_cell);
        for (GEO::index_t j = 0; j < verts_per_cell; ++j) {
            GEO::index_t v_reordered = static_cast<GEO::index_t>(simplices_ptr[i * verts_per_cell + j]);
            simplices[i].push_back(v_reordered);
        }
    }

    return simplices;
}

namespace py = pybind11;

PYBIND11_MODULE(diffvoronoi, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: diffvoronoi

        .. autosummary::
           :toctree: _generate

           get_delaunay_simplices
    )pbdoc";

    m.def("get_delaunay_simplices", &get_delaunay_simplices, R"pbdoc(
        get_delaunay_simplices 
    )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
