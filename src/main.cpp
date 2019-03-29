#include <misaxx-kidney-glomeruli/module.h>
#include <misaxx/core/runtime/misa_cli.h>
#include <misaxx-kidney-glomeruli/module_info.h>

using namespace misaxx;
using namespace misaxx_kidney_glomeruli;

int main(int argc, const char** argv) {
    misa_cli cli {};
    cli.set_module_info(misaxx_kidney_glomeruli::module_info());
    cli.set_root_module<misaxx_kidney_glomeruli::module>("misaxx-kidney-glomeruli");
    return cli.prepare_and_run(argc, argv);
}