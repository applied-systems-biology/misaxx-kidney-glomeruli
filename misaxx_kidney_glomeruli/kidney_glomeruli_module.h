//
// Created by rgerst on 17.08.18.
//


#pragma once

#include <misaxx/misa_module.h>
#include "quantification/constrained_klingberg.h"
#include "quantification/klingberg.h"
#include "segmentation2d/klingberg.h"
#include "segmentation2d/local_otsu.h"
#include "segmentation2d/segmentation2d_base.h"
#include "segmentation3d/klingberg.h"
#include "kidney_glomeruli.h"

namespace misaxx::module::kidney_glomeruli_detection {

    struct kidney_glomeruli_module : public misa_module<kidney_glomeruli> {
        using misa_module::misa_module;

        dispatched <tissue_detection::tissue_detection_module> dispatch_tissue_detection = future_dispatch(m_tissue_detection);

        dispatched <segmentation2d::segmentation2d_base> dispatch_segmentation2d =
                future_dispatch_any_from_algorithm_json_or<segmentation2d::segmentation2d_base>("segmentation2d",
                                                                                                "segmentation2d-klingberg",
                                                                                                option<segmentation2d::klingberg>("segmentation2d-klingberg"),
                                                                                                option<segmentation2d::local_otsu>("segmentation2d-local_otsu"));

        dispatched <segmentation3d::segmentation3d_base> dispatch_segmentation3d =
                future_dispatch_any_from_algorithm_json_or<segmentation3d::segmentation3d_base>("segmentation3d",
                                                                                                "segmentation3d-klingberg",
                                                                                                option<segmentation3d::klingberg>("segmentation3d-klingberg"));
        dispatched <quantification::quantification_base> dispatch_quantification =
                future_dispatch_any_from_algorithm_json_or<quantification::quantification_base>("quantification",
                                                                                                "quantification-klingberg",
                                                                                                option<quantification::klingberg>("quantification-klingberg"),
                                                                                                option<quantification::constrained_klingberg>("quantification-constrained_klingberg"));

        void misa_init() override {

            if (module().m_input_autofluorescence->files.empty())
                return;

            group preprocessing;
            preprocessing << misa_dispatch(dispatch_tissue_detection);

            group segmentation2d({{preprocessing}});

            for (auto &kv : *module().m_input_autofluorescence) {
                auto &worker = misa_dispatch(dispatch_segmentation2d);
                worker.m_input_autofluoresence = kv.second;
                worker.m_input_tissue = module().m_tissue_detection.instance().m_output_segmented3d->at(kv.first);
                worker.m_output_segmented2d = module().m_output_segmented2d->at(kv.first);
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

            work3d >> run_function([this](pattxx::functional_task &task) {
                this->module().m_output_quantification->save();
            });
        }

    };
}
