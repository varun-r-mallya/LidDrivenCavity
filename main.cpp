#include <omp.h>
#include "LidDrivenCavity.h"

int main() {
    omp_set_num_threads(omp_get_max_threads());
    LidDrivenCavity cavity(500.0, 2.0, 100, 0.01, 1000);
    cavity.solve();

    return 0;
}
