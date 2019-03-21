//
// Created by rgerst on 14.12.18.
//

#include <misaxx-kidney-glomeruli/module_interface.h>

using namespace misaxx;
using namespace misaxx_kidney_glomeruli;

glomeruli module_interface::get_glomeruli() {
    return m_output_quantification.get_attachment<glomeruli>();
}

void module_interface::setup() {
    m_input_autofluorescence.suggest_import_location(filesystem, "/");
    m_input_autofluorescence.suggest_document_title("Raw images");
    m_input_autofluorescence.suggest_document_description("Grayscale OME TIFF file containg the raw images");

    m_output_segmented2d.suggest_export_location(filesystem, "glomeruli2d", m_input_autofluorescence.derive().of_opencv(CV_8U));
    m_output_segmented2d.suggest_document_title("2D segmented glomeruli");
    m_output_segmented2d.suggest_document_description("Unfiltered glomeruli detected via 2D segmentation");

    m_output_segmented3d.suggest_export_location(filesystem, "glomeruli3d", m_output_segmented2d.derive().of_opencv(CV_32S));
    m_output_segmented3d.suggest_document_title("3D segmented glomeruli");
    m_output_segmented3d.suggest_document_description("2D detected glomeruli filtered by the 3D segmentation algorithm");

    m_output_quantification.suggest_export_location(filesystem, "quantified/quantified.json");
    m_output_quantification.suggest_document_title("Quantification results");

    // Init the submodule
    m_tissue = std::make_shared<misaxx_tissue::module_interface>();
    m_tissue->m_input_autofluorescence = m_input_autofluorescence;
}
