#include <misaxx-kidney-glomeruli/misaxx_kidney_glomeruli_module.h>
#include <misaxx/runtime/misa_cli.h>

using namespace misaxx;
using namespace misaxx_kidney_glomeruli;

int main(int argc, const char** argv) {
    misa_cli<misa_multiobject_root<misaxx_kidney_glomeruli_module>> cli("misaxx-kidney-glomeruli");
    return cli.prepare_and_run(argc, argv);
}