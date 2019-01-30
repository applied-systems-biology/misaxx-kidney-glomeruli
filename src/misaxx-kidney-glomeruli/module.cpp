//
// Created by rgerst on 14.12.18.
//

#include <misaxx-kidney-glomeruli/module.h>
#include <misaxx-tissue/module_interface.h>
#include <misaxx-tissue/module.h>
#include "algorithms/segmentation2d/segmentation2d_klingberg.h"
#include "algorithms/segmentation2d/segmentation2d_local_otsu.h"
#include "algorithms/segmentation3d/segmentation3d_klingberg.h"
#include "algorithms/quantification/quantification_klingberg.h"
#include "algorithms/quantification/quantification_constrained_klingberg.h"

using namespace misaxx;
using namespace misaxx_kidney_glomeruli;

void module::create_blueprints(misa_dispatcher::blueprint_list &t_blueprints,
                                                       misa_dispatcher::parameter_list &t_parameters) {

    m_segmentation2d_algorithm = t_parameters.create_algorithm_parameter<std::string>("segmentation2d");
    m_segmentation3d_algorithm = t_parameters.create_algorithm_parameter<std::string>("segmentation3d");
    m_quantification_algorithm = t_parameters.create_algorithm_parameter<std::string>("quantification");

    t_blueprints.add(create_submodule_blueprint<misaxx_tissue::module>("tissue", get_module_as<module_interface>()->m_tissue));
    t_blueprints.add(create_blueprint_enum_parameter(m_segmentation2d_algorithm, {
        create_blueprint<segmentation2d_klingberg>("segmentation2d_klingberg"),
        create_blueprint<segmentation2d_local_otsu>("segmentation2d_local_otsu")
    }, "segmentation2d_klingberg"));
    t_blueprints.add(create_blueprint_enum_parameter(m_segmentation3d_algorithm, {
            create_blueprint<segmentation3d_klingberg>("segmentation3d_klingberg")
    }, "segmentation3d_klingberg"));
    t_blueprints.add(create_blueprint_enum_parameter(m_quantification_algorithm, {
            create_blueprint<quantification_klingberg>("quantification_klingberg"),
            create_blueprint<quantification_constrained_klingberg>("quantification_constrained_klingberg")
    }, "quantification_klingberg"));

}

void module::build(const misa_dispatcher::blueprint_builder &t_builder) {

    m_voxel_size = misaxx::ome::misa_ome_voxel_size(*m_input_autofluorescence.get_ome_metadata(),
                                                   0, misaxx::ome::misa_ome_voxel_size::ome_unit_type::MICROMETER);

    group preprocessing;
    preprocessing << t_builder.build<misaxx_tissue::module>("tissue");

    group segmentation2d({{preprocessing}});

    for (size_t plane = 0; plane < m_input_autofluorescence.size(); ++plane) {
        auto &worker = t_builder.build<segmentation2d_base>(m_segmentation2d_algorithm.query());
        worker.m_input_autofluoresence = m_input_autofluorescence.at(plane);
        worker.m_input_tissue = m_tissue->m_output_segmented3d.at(plane);
        worker.m_output_segmented2d = m_output_segmented2d.at(plane);
        segmentation2d << worker;
    }

    chain work3d({segmentation2d});

    {
        auto &worker = t_builder.build<segmentation3d_base>(m_segmentation3d_algorithm.query());
        worker.m_input_segmented2d = m_output_segmented2d;
        worker.m_output_segmented3d = m_output_segmented3d;
        work3d >> worker;
    }
    {
        auto &worker = t_builder.build<quantification_base>(m_quantification_algorithm.query());
        worker.m_input_segmented3d = m_output_segmented3d;
        work3d >> worker;
    }
}
