//
// Created by rgerst on 17.08.18.
//


#pragma once

#include <misaxx/misa_module.h>
#include <misaxx_kidney_glomeruli/segmentation2d/klingberg.h>
#include "segmentation2d/segmentation2d_base.h"
#include "kidney_glomeruli.h"

namespace misaxx::module::kidney_glomeruli_detection {

    /**
     * Avaialable 2d segmentation algorithms
     */
//    struct segmentation2d_algorithms : public cxxh::types::options<std::string> {
//        static inline const values available_values = {{ "klingberg", "local_otsu" }};
//        static inline const value default_value = "klingberg";
//    };

//    using segmentation2d_algorithm = cxxh::types::option<segmentation2d_algorithms>;


    struct kidney_glomeruli_module : public misa_module<kidney_glomeruli> {
        using misa_module::misa_module;

        /**
         * Algorithm for 2D segmentation
         */
//        segmentation2d_algorithm m_segmentation2d_algorithm = from_algorithm_json_or<segmentation2d_algorithm>("segmentation2d-algorithm", segmentation2d_algorithms::default_value);

        void init() {

            // Autofluorescence contained directly in folder
            init_data(m_input_autofluorescence) << filesystem.imported;
            init_data(m_output_quantification) << filesystem.exported / filesystem::as<filesystem::file >("quantified.json");

            // Init submodule
            init_data(m_tissue_detection); // TODO: How to init?

            // Init dependent data
            init_data(m_output_segmented2d) << m_input_autofluorescence;
            init_data(m_output_segmented3d) << m_output_segmented2d;

            if(module().m_input_autofluorescence.files.empty())
                return;

            misa_dispatch(m_tissue_detection);

//            group segmentation2d;
//            for(auto &kv : module().m_input_autofluorescence) {
//                std::cout << kv.first << std::endl;
//                auto &worker = dispatch_segmentation2d();
//                worker.m_input_autofluoresence = kv.second;
//                worker.m_output_segmented2d = module().m_output_segmented2d.files.at(kv.first);
//                segmentation2d << worker;
//            }

//            auto &segmentation3d = misa_dispatch<segmentation3d::no_segmentation3d>("segmentation3d");
//            segmentation3d.m_input_segmented2d = m_output_segmented2d;
//            segmentation3d.m_output_segmented3d = m_output_segmented3d;
//
//            auto &quant = misa_dispatch<quantification>("quantification");
//            chain workflow;
//            workflow.add_dependency(segmentation2d);
//            workflow >> segmentation3d >> quant;
        }

        segmentation2d::segmentation2d_base &dispatch_segmentation2d() {
//            if(m_segmentation2d_algorithm.get() == "klingberg") {
                return misa_dispatch<segmentation2d::klingberg>("segmentation2d-klingberg");
//            }
//            else {
//                throw std::runtime_error("Unknown 2D segmentation algorithm");
//            }
        }
    };
}
