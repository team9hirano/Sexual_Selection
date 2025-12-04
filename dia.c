/* nearest neighbor interaction */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "MT.h"

#define LH 1000 // 10000
#define LV 1000 // 10000
#define K 0.005
#define u 0.3
#define a1 3.0
#define tend 20000 // 2000

int main(void)
{
    int x;
    double y, z;
    FILE *gp, *data1, *data2, *data3;
    char *data_file1, *data_file2, *data_file3;
    data_file1 = malloc(100);
    sprintf(data_file1, "Two_env_2dime_T2P2_K_%f.dat", K);

    data_file2 = malloc(100);
    sprintf(data_file2, "Two_env_2dime_T2P2_flow_K_%f.dat", K);

    data_file3 = malloc(100);
    sprintf(data_file3, "Two_env_3dime_final_T2P2_K_%f.dat", K);

    // T1P1-T2P2図
    gp = popen("gnuplot -persist", "w");
    fprintf(gp, "set terminal png\n");
    fprintf(gp, "set term pngcairo size 1000,700\n");
    fprintf(gp, "set output 'Genotype_Twoalleles/K_%f_a1_%f_T1P1_T2P2.png'\n", K, a1);
    fprintf(gp, "set xrange [0:%f]\n", 1.0);
    fprintf(gp, "set xlabel 'T1P1'\n");
    fprintf(gp, "set yrange [0:%f]\n", 1.0);
    fprintf(gp, "set ylabel 'T2P2'\n");
    // fprintf(gp, "set zrange [0:%f]\n", 1.0);
    // fprintf(gp, "set zlabel 'T2P1'\n");

    fprintf(gp, "plot \'%s\' using 2:5 with points pointtype 7 lc rgb 'blue' title \"survivalrateK=%f\",\'%s\' using 1:3:($4-$1):($6-$3) with vectors head filled lc rgb 'blue',\'%s\' using 2:5 with points pointtype 7 lc rgb 'red' title \"finalarrival\"\n", data_file1, K, data_file2, data_file3);
    // fprintf(gp, "splot \'%s\' using 2:5:4 with points pointtype 7 lc rgb 'blue' title \"survivalrateK=%f\",\'%s\' using 1:3:2:($4-$1):($6-$3):($5-$2) with vectors head filled lc rgb 'blue',\'%s\' using 2:5:4 with points pointtype 7 lc rgb 'red' title \"finalarrival\"\n", data_file1, K, data_file2, data_file3);

    pclose(gp);

    gp = popen("gnuplot -persist", "w");
    fprintf(gp, "set terminal png\n");
    fprintf(gp, "set term pngcairo size 1000,700\n");
    fprintf(gp, "set output 'Genotype_Twoalleles/K_%f_a1_%f_T1P1_T2P1.png'\n", K, a1);
    fprintf(gp, "set xrange [0:%f]\n", 1.0);
    fprintf(gp, "set xlabel 'T1P1'\n");
    fprintf(gp, "set yrange [0:%f]\n", 1.0);
    fprintf(gp, "set ylabel 'T2P1'\n");
    // fprintf(gp, "set zrange [0:%f]\n", 1.0);
    // fprintf(gp, "set zlabel 'T2P1'\n");

    fprintf(gp, "plot \'%s\' using 2:4 with points pointtype 7 lc rgb 'blue' title \"survivalrateK=%f\",\'%s\' using 1:2:($4-$1):($5-$2) with vectors head filled lc rgb 'blue',\'%s\' using 2:4 with points pointtype 7 lc rgb 'red' title \"finalarrival\"\n", data_file1, K, data_file2, data_file3);
    // fprintf(gp, "splot \'%s\' using 2:5:4 with points pointtype 7 lc rgb 'blue' title \"survivalrateK=%f\",\'%s\' using 1:3:2:($4-$1):($6-$3):($5-$2) with vectors head filled lc rgb 'blue',\'%s\' using 2:5:4 with points pointtype 7 lc rgb 'red' title \"finalarrival\"\n", data_file1, K, data_file2, data_file3);

    pclose(gp);
}
