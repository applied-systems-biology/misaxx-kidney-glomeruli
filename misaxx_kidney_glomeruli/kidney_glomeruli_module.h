//
// Created by rgerst on 17.08.18.
//


#pragma once

#include <misaxx/misa_module.h>
#include "quantification/klingberg.h"
#include "segmentation2d/klingberg.h"
#include "segmentation2d/local_otsu.h"
#include "segmentation2d/segmentation2d_base.h"
#include "segmentation3d/klingberg.h"
#include "kidney_glomeruli.h"

namespace misaxx::module::kidney_glomeruli_detection {

    struct segmentation2d_algorithms : public cxxh::types::options<std::string> {
        static inline const values available_values = { "klingberg", "local_otsu" };
        static inline const value default_value = "klingberg";
    };

    struct segmentation3d_algorithms : public cxxh::types::options<std::string> {
        static inline const values available_values = { "klingberg" };
        static inline const value default_value = "klingberg";
    };

    struct quantification_algorithms : public cxxh::types::options<std::string> {
        static inline const values available_values = { "klingberg" };
        static inline const value default_value = "klingberg";
    };

    using segmentation2d_algorithm = cxxh::types::option<segmentation2d_algorithms>;
    using segmentation3d_algorithm = cxxh::types::option<segmentation3d_algorithms>;
    using quantification_algorithm = cxxh::types::option<quantification_algorithms>;

    struct kidney_glomeruli_module : public misa_module<kidney_glomeruli> {
        using misa_module::misa_module;

        /**
         * Algorithm for 2D segmentation
         */
        segmentation2d_algorithm m_segmentation2d_algorithm = from_algorithm_json_or<segmentation2d_algorithm>("segmentation2d-algorithm", segmentation2d_algorithm());
        segmentation3d_algorithm m_segmentation3d_algorithm = from_algorithm_json_or<segmentation3d_algorithm>("segmentation3d-algorithm", segmentation3d_algorithm());
        quantification_algorithm m_quantification_algorithm = from_algorithm_json_or<quantification_algorithm>("quantification-algorithm", quantification_algorithm());

        void init() {

            if(module().m_input_autofluorescence->files.empty())
                return;

            group preprocessing;
            preprocessing << misa_dispatch(m_tissue_detection);

            group segmentation2d({{ preprocessing }});

            for(auto &kv : *module().m_input_autofluorescence) {
                auto &worker = dispatch_segmentation2d();
                worker.m_input_autofluoresence = kv.second;
                worker.m_input_tissue = module().m_tissue_detection.instance().m_output_segmented3d->at(kv.first);
                worker.m_output_segmented2d = module().m_output_segmented2d->at(kv.first);
                segmentation2d << worker;
            }

            chain work3d({segmentation2d});

            {
                auto &worker = dispatch_segmentation3d();
                worker.m_input_segmented2d = m_output_segmented2d;
                worker.m_output_segmented3d = m_output_segmented3d;
                work3d >> worker;
            }
            {
                auto &worker = dispatch_quantification();
                worker.m_input_segmented3d = m_output_segmented3d;
                work3d >> worker;
            }

            work3d >> run_function([this](pattxx::functional_task &task) {
                this->module().m_output_quantification->save();
            });
        }

        segmentation2d::segmentation2d_base &dispatch_segmentation2d() {
            if(m_segmentation2d_algorithm.get() == "klingberg") {
                return misa_dispatch<segmentation2d::klingberg>("segmentation2d-klingberg");
            }
            else  if(m_segmentation2d_algorithm.get() == "local_otsu") {
                return misa_dispatch<segmentation2d::local_otsu>("segmentation2d-local_otsu");
            }
            else {
                throw std::runtime_error("Unknown 2D segmentation algorithm");
            }
        }

        segmentation3d::segmentation3d_base &dispatch_segmentation3d() {
            if(m_segmentation3d_algorithm.get() == "klingberg") {
                return misa_dispatch<segmentation3d::klingberg>("segmentation3d-klingberg");
            }
            else {
                throw std::runtime_error("Unknown 3D segmentation algorithm");
            }
        }

        quantification::quantification_base &dispatch_quantification() {
            if(m_quantification_algorithm.get() == "klingberg") {
                return misa_dispatch<quantification::klingberg>("quantification-klingberg");
            }
            else {
                throw std::runtime_error("Unknown quantification algorithm");
            }
        }
    };
}
