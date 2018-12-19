//
// Created by rgerst on 17.08.18.
//


#pragma once

#include <misaxx/misa_module.h>
#include <misaxx-kidney-glomeruli/algorithms/segmentation2d/segmentation2d_base.h>
#include <misaxx-kidney-glomeruli/algorithms/segmentation3d/segmentation3d_base.h>
#include <misaxx-kidney-glomeruli/algorithms/quantification/quantification_base.h>
#include <misaxx-kidney-glomeruli/algorithms/segmentation2d/segmentation2d_klingberg.h>
#include <misaxx-kidney-glomeruli/algorithms/segmentation2d/segmentation2d_local_otsu.h>
#include <misaxx-kidney-glomeruli/algorithms/quantification/quantification_klingberg.h>
#include <misaxx-kidney-glomeruli/algorithms/quantification/quantification_constrained_klingberg.h>
#include <misaxx-kidney-glomeruli/algorithms/segmentation3d/segmentation3d_klingberg.h>
#include "kidney_glomeruli.h"

namespace misaxx_kidney_glomeruli {

    struct misaxx_kidney_glomeruli_module : public misaxx::misa_module<kidney_glomeruli> {
        using misaxx::misa_module<kidney_glomeruli>::misa_module;

        dispatched <misaxx_tissue::misaxx_tissue_module> dispatch_tissue_detection = future_dispatch(m_tissue_detection);

        dispatched <segmentation2d_base> dispatch_segmentation2d =
                select_from_algorithm_json_or<segmentation2d_base>("segmentation2d",
                                                                   "segmentation2d-klingberg", {
                                                                           future_dispatch<segmentation2d_klingberg>("segmentation2d-klingberg"),
                                                                           future_dispatch<segmentation2d_local_otsu>("segmentation2d-local_otsu")
                                                                   }
                );

        dispatched <segmentation3d_base> dispatch_segmentation3d =
                select_from_algorithm_json_or<segmentation3d_base>("segmentation3d",
                                                                   "segmentation3d-klingberg",
                                                                   {
                                                                           future_dispatch<segmentation3d_klingberg>("segmentation3d-klingberg")
                                                                   });
        dispatched <quantification_base> dispatch_quantification =
                select_from_algorithm_json_or<quantification_base>("quantification",
                                                                   "quantification-klingberg", {
                                                                           future_dispatch<quantification_klingberg>("quantification-klingberg"),
                                                                           future_dispatch<quantification_constrained_klingberg>("quantification-constrained_klingberg")
                                                                   });

        void misa_init() override;

    };
}
