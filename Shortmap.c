/* nearest neighbor interaction */
#define _POSIX_C_SOURCE 199309L
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <time.h>
#include "MT.h"

#define LH 3 // 1000
#define LV 3 // 1000
// #define K 0.075 //P3メスのコスト
// #define V 0.074 //P2メスのコスト(0<V<K)
// #define u 0.3 // T3オスのコスト0.3
// #define l 0.15  //T2オスのコスト(0<l<u)
#define a1 3.0 // P2メスがT2オスを選好する倍率3.0
// #define a2 6.0    // P3メスがT3オスを選好する倍率
#define tend 300 // 4000 80000 10000 70000 30000
#define mapinitP 0.25
#define initialP 3
#define initialT 1
#define SAVE_INTERVAL 1 // 100世代ごとに書き出し
#define MAX_SAVE ((tend / SAVE_INTERVAL) + 2)
#define N_SEED 5
// #define MAX_THREADS 12

void Map(const char *sex, const char *filename, double K, double initP, int t)
{
    FILE *gp;
    gp = popen("gnuplot -persist", "w");
    fprintf(gp, "set term pngcairo size 1000,1000\n");
    //  fprintf(gp,"set terminal png\n");
    if (strcmp(sex, "male") == 0)
    {
        if (strstr(filename, "stmap"))
            fprintf(gp, "set output 'Genotype_twoalleles_map/MapSt_%s_two_%g_initP_%g.png'\n", sex, K, initP);
        else if (strstr(filename, "intmap"))
            fprintf(gp, "set output 'Genotype_twoalleles_map/Shortmap_K_0.05/MapInt_200_%s_two_%g_initP_%g_t_%06d.png'\n", sex, K, initP, t);
        else if (strstr(filename, "finmap"))
            fprintf(gp, "set output 'Genotype_twoalleles_map/MapFin_%s_two_%g_initP_%g.png'\n", sex, K, initP);
        fprintf(gp, "unset key\n");
        fprintf(gp, "set size ratio -1\n");
        fprintf(gp, "set xrange [0:%d]\n", 200); // LH - 1
        // fprintf(gp,"set xlabel 'T2'\n");
        fprintf(gp, "set yrange [0:%d]\n", 200); // LV - 1
        fprintf(gp, "set palette defined(1 'blue',2 'green',3 'orange',4 'red')\n");
        fprintf(gp, "set cbtics ('T1P1' 1, 'T1P2' 2, 'T2P1' 3, 'T2P2' 4)\n");
        fprintf(gp, "unset autoscale cb\n");
        fprintf(gp, "set cbrange [1:4]\n");
        // fprintf(gp,"set multiplot layout 1,2 title 'Genotype map (T×P: 0=T1P1, 1=T1P2, 2=T2P1, 3=T2P2)'\n");
        fprintf(gp, "set title 'Male map'\n");
        fprintf(gp, "unset xtics;unset ytics\n");
        // fprintf(gp,"unset yticks\n");
        fprintf(gp, "plot \'%s\' using 1:2:3 with image\n", filename);
    }
    else if (strcmp(sex, "female") == 0)
    {
        if (strstr(filename, "stmap"))
            fprintf(gp, "set output 'Genotype_twoalleles_map/MapSt_%s_two_%g_initP_%g.png'\n", sex, K, initP);
        else if (strstr(filename, "intmap"))
            fprintf(gp, "set output 'Genotype_twoalleles_map/MapInt_%s_two_%g_initP_%g_t_%06d.png'\n", sex, K, initP, t);
        else if (strstr(filename, "finmap"))
            fprintf(gp, "set output 'Genotype_twoalleles_map/MapFin_%s_two_%g_initP_%g.png'\n", sex, K, initP);
        fprintf(gp, "unset key\n");
        fprintf(gp, "set size ratio -1\n");
        fprintf(gp, "set xrange [0:%d]\n", LH - 1);
        // fprintf(gp,"set xlabel 'T2'\n");
        fprintf(gp, "set yrange [0:%d]\n", LV - 1);
        fprintf(gp, "set palette defined(1 'blue',2 'green',3 'orange',4 'red')\n");
        fprintf(gp, "set cbtics ('T1P1' 1, 'T1P2' 2, 'T2P1' 3, 'T2P2' 4)\n");
        fprintf(gp, "unset autoscale cb\n");
        fprintf(gp, "set cbrange [1:4]\n");
        // fprintf(gp,"set multiplot layout 1,2 title 'Genotype map (T×P: 0=T1P1, 1=T1P2, 2=T2P1, 3=T2P2)'\n");
        fprintf(gp, "set title 'Female map'\n");
        fprintf(gp, "unset xtics;unset ytics\n");
        // fprintf(gp,"unset yticks\n");
        fprintf(gp, "plot \'%s\' using 1:2:4 with image\n", filename);
    }
    else
        printf("それはだめよ");

    pclose(gp);
}

static inline void calc_male_sum(double *sum, int i, int j, int **restrict maleT, int **restrict maleP, int female)
{
    int di[5] = {-1, 0, 1, 0, 0};
    int dj[5] = {0, 1, 0, -1, 0};
    int i2, j2, n, t, p, id;
    double w;
    for (n = 0; n < 4; n++)
        sum[n] = 0.0;
    // printf("OK1");
    for (n = 0; n < 5; n++)
    {
        i2 = (i + di[n] + LH) % LH;
        j2 = (j + dj[n] + LV) % LV;
        // printf("OK2");
        t = maleT[i2][j2];
        p = maleP[i2][j2];
        // printf("OK3");
        w = 1.0;
        if (female == 2 && t == 2)
            w = a1;
        // printf("T=0");
        // if(t==0)printf("T=0");
        id = 2 * (t - 1) + p - 1;
        sum[id] += w;
    }
}

static inline void calc_female_sum(double *sum, int i, int j, int **restrict femaleT, int **restrict femaleP)
{
    int di[5] = {-1, 0, 1, 0, 0};
    int dj[5] = {0, 1, 0, -1, 0};
    int i2, j2, n, t, p, id;
    for (n = 0; n < 4; n++)
        sum[n] = 0.0;
    for (n = 0; n < 5; n++)
    {
        i2 = (i + di[n] + LH) % LH;
        j2 = (j + dj[n] + LV) % LV;
        t = femaleT[i2][j2];
        p = femaleP[i2][j2];
        // if(t==0)printf("T=0");
        id = 2 * (t - 1) + p - 1;
        sum[id] += 1.0;
    }
}

static inline void genotype(double *sum, int *sexT0, int *sexP0, mt_state *rng_states)
{
    int n;
    // int tid = omp_get_thread_num();
    double rnd = genrand_real2_mt(rng_states);
    double total, acc;
    total = 0.0;
    acc = 0.0;
    for (n = 0; n < 4; n++)
        total += sum[n];
    if (total == 0.0)
    {
        *sexT0 = 1;
        *sexP0 = 1;
        // printf("fault");
        return;
    }
    if (rnd < sum[0] / total)
    {
        *sexT0 = 1;
        *sexP0 = 1;
    }
    else if (rnd < (sum[0] + sum[1]) / total)
    {
        *sexT0 = 1;
        *sexP0 = 2;
    }
    else if (rnd < (sum[0] + sum[1] + sum[2]) / total)
    {
        *sexT0 = 2;
        *sexP0 = 1;
    }
    else
    {
        *sexT0 = 2;
        *sexP0 = 2;
    }
}

static inline void calc_pair(double *sum_pair, int i, int j, int **restrict maleT, int **restrict maleP)
{
    int di[4] = {-1, 0, 1, 0};
    int dj[4] = {0, 1, 0, -1};
    int i2, j2, n, t, p, id;
    double w;

    for (n = 0; n < 4; n++)
    {
        i2 = (i + di[n] + LH) % LH;
        j2 = (j + dj[n] + LV) % LV;
        // printf("OK2");
        t = maleT[i2][j2];
        p = maleP[i2][j2];
        // printf("OK3");
        if (maleT[i][j] == 1 && maleP[i][j] == 1)
        {
            if (t == 1 && p == 1)
                sum_pair[0] += 1.0;
            else if (t == 2 && p == 1)
                sum_pair[1] += 1.0;
            else if (t == 2 && p == 2)
                sum_pair[2] += 1.0;
        }
        else if (maleT[i][j] == 2 && maleP[i][j] == 1)
        {
            if (t == 2 && p == 1)
                sum_pair[3] += 1.0;
            else if (t == 2 && p == 2)
                sum_pair[4] += 1.0;
        }
        else if (maleT[i][j] == 2 && maleP[i][j] == 2)
        {
            if (t == 2 && p == 2)
                sum_pair[5] += 1.0;
        }
    }
}

static inline void Initialize(int initT2P2, int **restrict maleT, int **restrict maleP, int **restrict femaleT, int **restrict femaleP)
{
    int di[5] = {-1, 0, 1, 0, 0};
    int dj[5] = {0, 1, 0, -1, 0};
    int i2, j2, n, t, p, id, i, j;
    double w;

    // printf("OK2");
    // t = maleT[i2][j2];
    // p = maleP[i2][j2];
    // printf("OK3");
    if (initT2P2 == 0)
    {
        for (i = 0; i < LH; i++)
        {
            for (j = 0; j < LV; j++)
            {
                if (i == LH / 2 && j == LV / 2)
                {
                    maleT[i][j] = 2;
                    maleP[i][j] = 2;
                    femaleT[i][j] = 1;
                    femaleP[i][j] = 1;
                }
                else
                {
                    maleT[i][j] = 1;
                    maleP[i][j] = 1;
                    femaleT[i][j] = 1;
                    femaleP[i][j] = 1;
                }
            }
        }
    }
    else if (initT2P2 == 1)
    {
        for (i = 0; i < LH; i++)
        {
            for (j = 0; j < LV; j++)
            {
                if (i == LH / 2 && j == LV / 2)
                {
                    maleT[i][j] = 2;
                    maleP[i][j] = 2;
                    femaleT[i][j] = 1;
                    femaleP[i][j] = 1;
                }
                else if (i == (LH / 2 + di[0]) % LH && j == (LV / 2 + dj[0]) % LV)
                {
                    maleT[i][j] = 2;
                    maleP[i][j] = 2;
                    femaleT[i][j] = 2;
                    femaleP[i][j] = 2;
                }
                else
                {
                    maleT[i][j] = 1;
                    maleP[i][j] = 1;
                    femaleT[i][j] = 1;
                    femaleP[i][j] = 1;
                }
            }
        }
    }
    else if (initT2P2 == 2)
    {
        for (i = 0; i < LH; i++)
        {
            for (j = 0; j < LV; j++)
            {
                if (i == LH / 2 && j == LV / 2)
                {
                    maleT[i][j] = 2;
                    maleP[i][j] = 2;
                    femaleT[i][j] = 1;
                    femaleP[i][j] = 1;
                }
                else if ((i == (LH / 2 + di[0]) % LH && j == (LV / 2 + dj[0]) % LV) || i == (LH / 2 + di[1]) % LH && j == (LV / 2 + dj[1]) % LV)
                {
                    maleT[i][j] = 2;
                    maleP[i][j] = 2;
                    femaleT[i][j] = 2;
                    femaleP[i][j] = 2;
                }
                else
                {
                    maleT[i][j] = 1;
                    maleP[i][j] = 1;
                    femaleT[i][j] = 1;
                    femaleP[i][j] = 1;
                }
            }
        }
    }
    else if (initT2P2 == 3)
    {
        for (i = 0; i < LH; i++)
        {
            for (j = 0; j < LV; j++)
            {
                if (i == LH / 2 && j == LV / 2)
                {
                    maleT[i][j] = 2;
                    maleP[i][j] = 2;
                    femaleT[i][j] = 1;
                    femaleP[i][j] = 1;
                }
                else if ((i == (LH / 2 + di[0]) % LH && j == (LV / 2 + dj[0]) % LV) || i == (LH / 2 + di[1]) % LH && j == (LV / 2 + dj[1]) % LV || (i == (LH / 2 + di[2]) % LH && j == (LV / 2 + dj[2]) % LV))
                {
                    maleT[i][j] = 2;
                    maleP[i][j] = 2;
                    femaleT[i][j] = 2;
                    femaleP[i][j] = 2;
                }
                else
                {
                    maleT[i][j] = 1;
                    maleP[i][j] = 1;
                    femaleT[i][j] = 1;
                    femaleP[i][j] = 1;
                }
            }
        }
    }
    else if (initT2P2 == 4)
    {
        for (i = 0; i < LH; i++)
        {
            for (j = 0; j < LV; j++)
            {
                if (i == LH / 2 && j == LV / 2)
                {
                    maleT[i][j] = 2;
                    maleP[i][j] = 2;
                    femaleT[i][j] = 1;
                    femaleP[i][j] = 1;
                }
                else if ((i == (LH / 2 + di[0]) % LH && j == (LV / 2 + dj[0]) % LV) || i == (LH / 2 + di[1]) % LH && j == (LV / 2 + dj[1]) % LV || (i == (LH / 2 + di[2]) % LH && j == (LV / 2 + dj[2]) % LV) || (i == (LH / 2 + di[3]) % LH && j == (LV / 2 + dj[3]) % LV))
                {
                    maleT[i][j] = 2;
                    maleP[i][j] = 2;
                    femaleT[i][j] = 2;
                    femaleP[i][j] = 2;
                }
                else
                {
                    maleT[i][j] = 1;
                    maleP[i][j] = 1;
                    femaleT[i][j] = 1;
                    femaleP[i][j] = 1;
                }
            }
        }
    }
    else if (initT2P2 == 5)
    {
        for (i = 0; i < LH; i++)
        {
            for (j = 0; j < LV; j++)
            {
                if (i == LH / 2 && j == LV / 2)
                {
                    maleT[i][j] = 2;
                    maleP[i][j] = 2;
                    femaleT[i][j] = 2;
                    femaleP[i][j] = 2;
                }
                else if ((i == (LH / 2 + di[0]) % LH && j == (LV / 2 + dj[0]) % LV) || i == (LH / 2 + di[1]) % LH && j == (LV / 2 + dj[1]) % LV || (i == (LH / 2 + di[2]) % LH && j == (LV / 2 + dj[2]) % LV) || (i == (LH / 2 + di[3]) % LH && j == (LV / 2 + dj[3]) % LV))
                {
                    maleT[i][j] = 2;
                    maleP[i][j] = 2;
                    femaleT[i][j] = 2;
                    femaleP[i][j] = 2;
                }
                else
                {
                    maleT[i][j] = 1;
                    maleP[i][j] = 1;
                    femaleT[i][j] = 1;
                    femaleP[i][j] = 1;
                }
            }
        }
    }
}

int main(void)
{
    // 新しく挿入
    //  int tend;
    int **restrict maleT, *base_maleT;
    int **restrict maleP, *base_maleP;
    int **restrict femaleT, *base_femaleT;
    int **restrict femaleP, *base_femaleP;
    int **restrict maleTdummy, *base_maleTdummy;
    int **restrict malePdummy, *base_malePdummy;
    int **restrict femaleTdummy, *base_femaleTdummy;
    int **restrict femalePdummy, *base_femalePdummy;
    double initT2, initP2, initT3, initP3, initT2P1;
    int k, k2, i, j, i2, j2, t, ok, x1, x2, a, b, n, initT2P2;
    int maleI, maleJ, femaleI, femaleJ;
    int numMT1, numMT2, numMT3, numMP1, numMP2, numMP3;
    int numFT1, numFT2, numFT3, numFP1, numFP2, numFP3;
    double sum[4], gsum[4], sum_pair[6];
    double rnd, rnd2, rnd3, rnd4, rnd5, sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9, gen1, gen2, gen3, gen4, geno1, geno2, geno3, geno4, init, init1;
    double gsum1, gsum2, gsum3, gsum4, gsum5, gsum6, gsum7, gsum8, gsum9;
    double s1, s2, s3, curs1, curs2, curs3;
    int iK, iV, il, ia2, ia1, iu, situ;
    double K, V, l, a2, u, k_max;
    double x11, x21, x22, x31, x32, x33;
    double nx11, nx21, nx22, nx31, nx32, nx33;
    int maleT0, maleP0, femaleT0, femaleP0, mgenotype, fgenotype, count;
    FILE *gp, *data1, *data2, *data3, *data4, *data5, *data6, *data7, *data8;
    FILE *snapshot1, *snapshot2, *snapshot3, *snapshot4, *snapshot5, *snapshot6;
    char *data_file1, *data_file2, *data_file3, *data_file4, *data_file5, *data_file6, *data_file7, *data_file8, *data_file9;
    char *snapshot_file1, *snapshot_file2, *snapshot_file3, *snapshot_file4, *snapshot_file5, *snapshot_file6;
    char *Figaxis[9] = {"x_11", "x_21", "x_31", "x_22", "x_32", "x_33", "T1P1", "T2P1", "T2P2"};
    struct timespec start, end;

    data_file9 = malloc(100);
    sprintf(data_file9, "Twoalleles_threshold_small.csv");

    for (iK = 0; iK <= 0; iK++)
    {
        // K = (double)(iK * 2 - 1) * 0.10;

        if (iK == 0)
            K = 0.05;
        else if (iK == 1)
            K = 0.05;
        else if (iK == 2)
            K = 0.1;
        else if (iK == 3)
            K = 0.13;
        else if (iK == 4)
            K = 0.14;
        else if (iK == 5)
            K = 0.15;
        else if (iK == 6)
            K = 0.16;
        else if (iK == 7)
            K = 0.2;
        else if (iK == 8)
            K = 0.3;
        else if (iK == 9)
            K = 0.5;
        for (i = 0; i < 50000; i++)
        {
            char data_file4[256];
            snprintf(data_file4, sizeof(data_file4),
                     "Two_intmap_t_%d_cost_%f_initT2P1_0.3.dat",
                     i, K);

            gp = fopen(data_file4, "r");
            if (!gp)
            {
                printf("ファイル読めん: %s\n", data_file4);
                continue;
            }

            Map("male", data_file4, K, mapinitP, i);
        }
        return 0;
    }
}
