//
// Created by rgerst on 14.12.18.
//

#include <misaxx_kidney_glomeruli/kidney_glomeruli_detection_module.h>

using namespace misaxx;
using namespace coixx;
using namespace misaxx_kidney_glomeruli;

void kidney_glomeruli_detection_module::misa_init() {

    // Load voxel size
    m_voxel_size = misaxx_ome::misa_ome_voxel_size(*m_input_autofluorescence.get_ome_metadata(),
            0, misaxx_ome::misa_ome_voxel_size::ome_unit_type::MICROMETER);

    group preprocessing;
    preprocessing << misa_dispatch(dispatch_tissue_detection);

    group segmentation2d({{preprocessing}});

    for (size_t plane = 0; plane < m_input_autofluorescence.size(); ++plane) {
        auto &worker = misa_dispatch(dispatch_segmentation2d);
        worker.m_input_autofluoresence = m_input_autofluorescence.at(plane);
//        worker.m_input_tissue = m_tissue_detection. // TODO
        worker.m_output_segmented2d = m_output_segmented2d.at(plane);
        segmentation2d << worker;
    }

    chain work3d({segmentation2d});

    {
        auto &worker = misa_dispatch(dispatch_segmentation3d);
        worker.m_input_segmented2d = m_output_segmented2d;
        worker.m_output_segmented3d = m_output_segmented3d;
        work3d >> worker;
    }
    {
        auto &worker = misa_dispatch(dispatch_quantification);
        worker.m_input_segmented3d = m_output_segmented3d;
        work3d >> worker;
    }
}
