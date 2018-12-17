//
// Created by rgerst on 17.12.18.
//

#include <misaxx_kidney_glomeruli/algorithms/segmentation2d/segmentation2d_klingberg.h>
#include <coixx/toolbox/toolbox_statistics.h>
#include <coixx/toolbox/toolbox_blur.h>
#include <coixx/toolbox/toolbox_normalize.h>
#include <coixx/toolbox/toolbox_morph.h>
#include <coixx/structuring_element.h>

using namespace misaxx;
using namespace misaxx_kidney_glomeruli;
using namespace coixx;

void segmentation2d_klingberg::misa_work() {
    using namespace coixx::toolbox;

    images::mask img_non_tissue_mask = m_input_tissue.clone();

    if(statistics::is_black(img_non_tissue_mask)) { //INFO: Not inverted yet
        // Instead save a black image
        m_output_segmented2d.write(images::mask(img_non_tissue_mask.get_size(), colors::mask::background()));
        return;
    }

    img_non_tissue_mask << values::invert(); // Now we mark all non-tissue

    images::grayscale_float img = m_input_autofluoresence.clone();
    auto img_original = img.clone();

    // Initial median filtering + normalization
    img << blur::median(m_median_filter_size) << normalize::by_max();

    // Generated parameters
    const double voxel_xy_area = module()->m_voxel_size.get_xy_area().get_value();
    int glomeruli_max_morph_disk_radius = static_cast<int>(m_glomeruli_max_rad / voxel_xy_area);
    int glomeruli_min_morph_disk_radius = static_cast<int>((m_glomeruli_min_rad / 2.0) / voxel_xy_area);

    // Morphological operation (opening)
    // Corresponds to only allowing objects > disk_size to be included
    // Also subtract the morph result from the initial to remove uneven background + normalize
    img << morph::tophat(structuring_element::ellipse(glomeruli_max_morph_disk_radius * 2)) << normalize::by_max();

    // We are first extracting tissue data
    auto img_tissue = img.clone();

    // Get rid of non-tissue
    img_tissue << values::set_where(colors::grayscale_float::black(), img_non_tissue_mask);

    // Only select glomeruli if the threshold is higher than 75-percentile of kidney tissue
    double percentile_tissue = statistics::get_percentile_as<double>(img, m_threshold_percentile);

    //////////////
    // Now working in uint8
    //////////////

    // Threshold the main image
    uchar otsu_threshold = 0;
    images::grayscale8u img_as8u = semantic_convert<images::mask>(img) << binarize::otsu(otsu_threshold);
    if((otsu_threshold / 255.0) > percentile_tissue * m_threshold_factor ) {

        // Get rid of non-tissue
        img_as8u << values::set_where(colors::grayscale8u::black(), img_non_tissue_mask);

        // Morphological operation (object should have min. radius)
        img_as8u << morph::open(structuring_element::ellipse(glomeruli_min_morph_disk_radius * 2));
    }
    else {
        img_as8u << values::set(colors::grayscale8u(0));
    }

    // Save the mask
    m_output_segmented2d.write(std::move(img_as8u));
}
