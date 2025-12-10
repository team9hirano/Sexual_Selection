/* nearest neighbor interaction */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "MT.h"

#define LH 1000 // 10000
#define LV 1000 // 10000

#define u 0.3

#define tend 70000 // 2000

int main(void)
{
    double initT2, initP2, initT3, initP3, initT2P1;
    int k, k2, i, j, i2, j2, t, ok, x1, x2, a, b, n;
    int maleI, maleJ, femaleI, femaleJ;
    int numMT1, numMT2, numMT3, numMP1, numMP2, numMP3;
    int numFT1, numFT2, numFT3, numFP1, numFP2, numFP3;
    double sum[4], gsum[4], sum_pair[6];
    double rnd, rnd2, rnd3, rnd4, rnd5, sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9, gen1, gen2, gen3, gen4, geno1, geno2, geno3, geno4, init, init1;
    double gsum1, gsum2, gsum3, gsum4, gsum5, gsum6, gsum7, gsum8, gsum9;
    int iK, iV, il, ia2, ia1;
    double K, V, l, a2, a1;
    double x11, x21, x22, x31, x32, x33;
    double nx11, nx21, nx22, nx31, nx32, nx33;
    FILE *gp, *data1, *data2, *data3, *data4, *data5, *data6, *data7;
    FILE *snapshot1, *snapshot2, *snapshot3, *snapshot4, *snapshot5, *snapshot6;
    char *data_file1, *data_file2, *data_file3, *data_file4, *data_file5, *data_file6, *data_file7;
    char *snapshot_file1, *snapshot_file2, *snapshot_file3, *snapshot_file4, *snapshot_file5, *snapshot_file6;
    char *Figaxis[9] = {"x_11", "x_21", "x_31", "x_22", "x_32", "x_33", "T1P1", "T2P1", "T2P2"};

    K = 0.0;
    V = K / 2;
    a2 = 3.0;
    for (iK = 0; iK <= 0; iK++)
    { // iK=1;iK<=3;iK++
        K = 0.03 + (double)0.01 * (double)iK;
        // K=0.11+(double)iK*0.001;

        // V=0.20+(double)iV*0.01;

        for (ia1 = 3; ia1 <= 3; ia1++)
        {
            // if(ia1==3)continue;
            // else a1=(double)ia1;
            a1 = (double)ia1;

            char data_file4[256];
            snprintf(data_file4, sizeof(data_file4),
                     "Two_env_2dime_pair_K_%f.dat",
                     K);

            gp = fopen(data_file4, "r");
            if (!gp)
            {
                printf("ファイル読めん: %s\n", data_file4);
                continue;
            }
            char data_file6[256];
            snprintf(data_file6, sizeof(data_file6),
                     "Two_env_2dime_pair_final_K_%f.dat",
                     K);

            gp = fopen(data_file6, "r");
            if (!gp)
            {
                printf("ファイル読めん: %s\n", data_file6);
                continue;
            }
            data_file5 = malloc(100);
            sprintf(data_file5, "Two_env_2dime_pair_flow_K_%f.dat", K);

            // data_file1= f("Three_env_T3P3_K_%f_V_%f_l_%f_a1_%f_a2_%f.dat", K,V,l,a1,a2);

            data4 = fopen(data_file4, "r");
            data5 = fopen(data_file5, "w");
            if (fscanf(data4, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &x1, &x11, &x21, &x31, &x22, &x32, &x33, &sum1, &sum3, &sum4, &init) != 11)
                return 1;
            while (fscanf(data4, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &x2, &nx11, &nx21, &nx31, &nx22, &nx32, &nx33, &gsum1, &gsum3, &gsum4, &init1) == 11)
            {
                if (fabs(init - init1) < 1e-9)
                {
                    fprintf(data5, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
                            x11, x21, x31, x22, x32, x33, sum1, sum3, sum4, nx11, nx21, nx31, nx22, nx32, nx33, gsum1, gsum3, gsum4);
                }
                x1 = x2;
                x11 = nx11;
                x21 = nx21;
                x31 = nx31;
                x22 = nx22;
                x32 = nx32;
                x33 = nx33;
                sum1 = gsum1;

                sum3 = gsum3;
                sum4 = gsum4;
                init = init1;
            }
            fclose(data4);
            fclose(data5);

            for (i = 0; i < 9; i++)
            {
                for (j = i + 1; j < 9; j++)
                {
                    gp = popen("gnuplot -persist", "w");
                    fprintf(gp, "set terminal png\n");
                    fprintf(gp, "set term pngcairo size 1000,700\n");
                    fprintf(gp, "set output 'Genotype_Twoalleles_pair/K_%f_a1_%f_%s_%s.png'\n", K, a1, Figaxis[i], Figaxis[j]);
                    fprintf(gp, "set xrange [0:%f]\n", 1.0);
                    fprintf(gp, "set xlabel \'%s\'\n", Figaxis[i]);
                    fprintf(gp, "set yrange [0:%f]\n", 1.0);
                    fprintf(gp, "set ylabel \'%s\'\n", Figaxis[j]);
                    fprintf(gp, "plot \'%s\' using %d:%d with points pointtype 7 lc rgb 'blue' title \
                    \"survivalrateK=%f\",\
                    \'%s\' using %d:%d:($%d-$%d):($%d-$%d) with vectors head filled lc rgb 'blue',\
                    \'%s\' using %d:%d with points pointtype 7 lc rgb 'red' title \"finalarrival\"\n",
                            data_file4, i + 2, j + 2,
                            K,
                            data_file5, i + 1, j + 1, i + 10, i + 1, j + 10, j + 1,
                            data_file6, i + 2, j + 2);
                    pclose(gp);
                }
            }
        }
    }
}
