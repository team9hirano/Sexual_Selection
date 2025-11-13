#include <stdio.h>
#include <omp.h>

int main() {
    omp_set_num_threads(10);
    printf("omp_get_max_threads() = %d\n", omp_get_max_threads());

    #pragma omp parallel
    {
        printf("Hello from thread %d / %d\n", omp_get_thread_num(), omp_get_num_threads());
    }
}
