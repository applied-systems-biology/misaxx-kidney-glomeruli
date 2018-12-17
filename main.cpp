//
// Created by rgerst on 30.07.18.
//

#include <misaxx/runtime/misa_cli.h>
#include <misaxx_kidney_glomeruli/kidney_glomeruli_detection_module.h>

using namespace misaxx;
using namespace misaxx_kidney_glomeruli;

int main(int argc, const char** argv) {
    misa_cli<misa_multiobject_root<kidney_glomeruli_detection_module>> cli("kidney-glomeruli");
    return cli.prepare_and_run(argc, argv);
}