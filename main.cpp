//
// Created by rgerst on 30.07.18.
//

#include <misaxx/misa_cli.h>
#include <misaxx_kidney_glomeruli/kidney_glomeruli_module.h>

using namespace misaxx;
using namespace misaxx::module::kidney_glomeruli_detection;

int main(int argc, const char** argv) {
    misa_cli<misa_multiobject_root<kidney_glomeruli_module>> cli("kidney-glomeruli-detection");
    return cli.prepare_and_run(argc, argv);
}