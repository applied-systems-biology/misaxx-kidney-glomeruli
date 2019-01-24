#include <misaxx-kidney-glomeruli/module.h>
#include <misaxx/core/runtime/misa_cli.h>
#include <misaxx-kidney-glomeruli/module_info.h>

using namespace misaxx;
using namespace misaxx_kidney_glomeruli;

int main(int argc, const char** argv) {
    misa_cli<module> cli(misaxx_kidney_glomeruli::module_info());
    return cli.prepare_and_run(argc, argv);
}