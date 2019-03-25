#include <misaxx/core/misa_mutable_module_info.h>
#include <misaxx/core/module_info.h>
#include <misaxx-kidney-glomeruli/module_info.h>
#include <misaxx-tissue/module_info.h>

misaxx::misa_module_info misaxx_kidney_glomeruli::module_info() {
    misaxx::misa_mutable_module_info info;
    info.set_id("misaxx-kidney-glomeruli");
    info.set_version("1.0.0");
    info.set_name("MISA++ Kidney Glomeruli Segmentation");
    info.set_description("Segments glomeruli within kidney LSFM images");
    info.add_author("Ruman Gerst");
    info.set_license("BSD-2-Clause");
    info.set_organization("Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Institute (HKI), Jena, Germany");
    info.set_url("https://asb-git.hki-jena.de/RGerst/misaxx-kidney-glomeruli/");

    info.add_dependency(misaxx::module_info());
    info.add_dependency(misaxx_tissue::module_info());
    return info;
}