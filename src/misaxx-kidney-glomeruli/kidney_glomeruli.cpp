//
// Created by rgerst on 14.12.18.
//

#include <misaxx-kidney-glomeruli/kidney_glomeruli.h>

using namespace misaxx;
using namespace coixx;
using namespace misaxx_kidney_glomeruli;

glomeruli kidney_glomeruli::get_glomeruli() {
    return m_output_quantification.get_attachment<glomeruli>();
}

void kidney_glomeruli::init_data() {

    m_input_autofluorescence.suggest_import_location(filesystem, "/");
    m_output_segmented2d.suggest_export_location(filesystem, "glomeruli2d", m_input_autofluorescence.derive().of_coixx<images::mask>());
    m_output_segmented3d.suggest_export_location(filesystem, "glomeruli3d", m_output_segmented2d.derive().of_coixx<images::labels>());
    m_output_quantification.suggest_export_location(filesystem, "quantified/quantified.json");

    // Init the submodule
    m_tissue_detection.definition().m_input_autofluorescence = m_input_autofluorescence;
    init_submodule(m_tissue_detection, "tissue");
}
