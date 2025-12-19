/* nearest neighbor interaction */
#define _POSIX_C_SOURCE 199309L
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MT.h"

#define LH 1000 // 10000
#define LV 1000 // 10000
// #define K 0.075 //P3メスのコスト
// #define V 0.074 //P2メスのコスト(0<V<K)
#define u 0.3 // T3オスのコスト0.3
// #define l 0.15  //T2オスのコスト(0<l<u)
// #define a1 3.0    // P2メスがT2オスを選好する倍率
// #define a2 6.0    // P3メスがT3オスを選好する倍率
#define tend 150000 // 4000 80000 10000 80000
#define mapinitP 0.5
#define initialP 3
#define initialT 1

void Map(const char *sex, const char *filename, double K, double V, double initP, int t)
{
    FILE *gp;
    gp = popen("gnuplot -persist", "w");
    fprintf(gp, "set term pngcairo size 1000,1000\n");
    //  fprintf(gp,"set terminal png\n");
    if (strcmp(sex, "male") == 0)
    {
        if (strstr(filename, "stmap"))
            fprintf(gp, "set output 'Genotype_Threealleles_map/ThreeMapSt_%s_env_K_%g_V_%g_initP_%g.png'\n", sex, K, V, initP);
        else if (strstr(filename, "intmap"))
            fprintf(gp, "set output 'Genotype_Threealleles_map/ThreeMapInt_%s_env_K_%g_V_%g_initP_%g_t_%06d.png'\n", sex, K, V, initP, t);
        else if (strstr(filename, "finmap"))
            fprintf(gp, "set output 'Genotype_Threealleles_map/ThreeMapFin_%s_env_K_%g_V_%g_initP_%g.png'\n", sex, K, V, initP);
        fprintf(gp, "unset key\n");
        fprintf(gp, "set size ratio -1\n");
        fprintf(gp, "set xrange [0:%d]\n", LH - 1);
        // fprintf(gp,"set xlabel 'T2'\n");
        fprintf(gp, "set yrange [0:%d]\n", LV - 1);
        fprintf(gp, "set palette defined(1 \"#FFFFFF\", 2 \"#02befcff\", 3 \"#0000FF\", 4 \"#00FF00\", \
             5 \"#FFFF00\", 6 \"#FFA500\", 7 \"#FF0000\", 8 \"#FF69B4\", 9 \"#000000\")\n");

        fprintf(gp, "set cbtics ('T1P1' 1, 'T1P2' 2, 'T1P3' 3, 'T2P1' 4, \
            'T2P2' 5, 'T2P3' 6, 'T3P1' 7, 'T3P2' 8, 'T3P3' 9)\n");
        fprintf(gp, "unset autoscale cb\n");
        fprintf(gp, "set cbrange [1:9]\n");
        // fprintf(gp,"set multiplot layout 1,2 title 'Genotype map (T×P: 0=T1P1, 1=T1P2, 2=T2P1, 3=T2P2)'\n");
        fprintf(gp, "set title 'Male map'\n");
        fprintf(gp, "unset xtics;unset ytics\n");
        // fprintf(gp,"unset yticks\n");
        fprintf(gp, "plot \'%s\' using 1:2:3 with image\n", filename);
    }
    else if (strcmp(sex, "female") == 0)
    {

        if (strstr(filename, "stmap"))
            fprintf(gp, "set output 'Genotype_Threealleles_map/ThreeMapSt_%s_env_K_%g_V_%g_initP_%g.png'\n", sex, K, V, initP);
        else if (strstr(filename, "intmap"))
            fprintf(gp, "set output 'Genotype_Threealleles_map/ThreeMapInt_%s_env_K_%g_V_%g_initP_%g_t_%06d.png'\n", sex, K, V, initP, t);
        else if (strstr(filename, "finmap"))
            fprintf(gp, "set output 'Genotype_Threealleles_map/ThreeMapFin_%s_env_K_%g_V_%g_initP_%g.png'\n", sex, K, V, initP);
        fprintf(gp, "unset key\n");
        fprintf(gp, "set size ratio -1\n");
        fprintf(gp, "set xrange [0:%d]\n", LH - 1);
        // fprintf(gp,"set xlabel 'T2'\n");
        fprintf(gp, "set yrange [0:%d]\n", LV - 1);
        fprintf(gp, "set palette defined(1 \"#FFFFFF\", 2 \"#02befcff\", 3 \"#0000FF\", 4 \"#00FF00\", \
             5 \"#FFFF00\", 6 \"#FFA500\", 7 \"#FF0000\", 8 \"#FF69B4\", 9 \"#000000\")\n");

        fprintf(gp, "set cbtics ('T1P1' 1, 'T1P2' 2, 'T1P3' 3, 'T2P1' 4, \
            'T2P2' 5, 'T2P3' 6, 'T3P1' 7, 'T3P2' 8, 'T3P3' 9)\n");
        fprintf(gp, "unset autoscale cb\n");
        fprintf(gp, "set cbrange [1:9]\n");
        // fprintf(gp,"set multiplot layout 1,2 title 'Genotype map (T×P: 0=T1P1, 1=T1P2, 2=T2P1, 3=T2P2)'\n");
        fprintf(gp, "set title 'Male map'\n");
        fprintf(gp, "unset xtics;unset ytics\n");
        // fprintf(gp,"unset yticks\n");
        fprintf(gp, "plot \'%s\' using 1:2:3 with image\n", filename);
    }
    else
        printf("それはだめよ");

    pclose(gp);
}

void calc_male_sum(double *sum, int i, int j, int **maleT, int **maleP, int female, double a1, double a2)
{
    int di[5] = {-1, 0, 1, 0, 0};
    int dj[5] = {0, 1, 0, -1, 0};
    int i2, j2, n, t, p, id;
    double w;
    for (n = 0; n < 9; n++)
        sum[n] = 0.0;
    for (n = 0; n < 5; n++)
    {
        i2 = (i + di[n] + LH) % LH;
        j2 = (j + dj[n] + LV) % LV;
        t = maleT[i2][j2];
        p = maleP[i2][j2];
        w = 1.0;
        if (female == 2 && t == 2)
            w = a1;
        else if (female == 3 && t == 3)
            w = a2;
        if (t == 0)
            printf("T=0");
        id = 3 * (t - 1) + p - 1;
        sum[id] += w;
    }
}

void calc_female_sum(double *sum, int i, int j, int **femaleT, int **femaleP)
{
    int di[5] = {-1, 0, 1, 0, 0};
    int dj[5] = {0, 1, 0, -1, 0};
    int i2, j2, n, t, p, id;
    for (n = 0; n < 9; n++)
        sum[n] = 0.0;
    for (n = 0; n < 5; n++)
    {
        i2 = (i + di[n] + LH) % LH;
        j2 = (j + dj[n] + LV) % LV;
        t = femaleT[i2][j2];
        p = femaleP[i2][j2];
        if (t == 0)
            printf("T=0");
        id = 3 * (t - 1) + p - 1;
        sum[id] += 1.0;
    }
}

void genotype(double *sum, int *sexT0, int *sexP0, mt_state *rng_states)
{
    int n;
    double rnd = genrand_real2_mt(rng_states);
    double total, acc;
    total = 0.0;
    acc = 0.0;
    for (n = 0; n < 9; n++)
        total += sum[n];
    if (total == 0.0)
    {
        *sexT0 = 1;
        *sexP0 = 1;
        printf("fault");
        return;
    }
    //
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
        *sexT0 = 1;
        *sexP0 = 3;
    }
    else if (rnd < (sum[0] + sum[1] + sum[2] + sum[3]) / total)
    {
        *sexT0 = 2;
        *sexP0 = 1;
    }
    else if (rnd < (sum[0] + sum[1] + sum[2] + sum[3] + sum[4]) / total)
    {
        *sexT0 = 2;
        *sexP0 = 2;
    }
    else if (rnd < (sum[0] + sum[1] + sum[2] + sum[3] + sum[4] + sum[5]) / total)
    {
        *sexT0 = 2;
        *sexP0 = 3;
    }
    else if (rnd < (sum[0] + sum[1] + sum[2] + sum[3] + sum[4] + sum[5] + sum[6]) / total)
    {
        *sexT0 = 3;
        *sexP0 = 1;
    }
    else if (rnd < (sum[0] + sum[1] + sum[2] + sum[3] + sum[4] + sum[5] + sum[6] + sum[7]) / total)
    {
        *sexT0 = 3;
        *sexP0 = 2;
    }
    else
    {
        *sexT0 = 3;
        *sexP0 = 3;
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
            else if (t == 3 && p == 3)
                sum_pair[3] += 1.0;
        }
        else if (maleT[i][j] == 2 && maleP[i][j] == 1)
        {
            if (t == 2 && p == 1)
                sum_pair[4] += 1.0;
            else if (t == 2 && p == 2)
                sum_pair[5] += 1.0;
            else if (t == 3 && p == 3)
                sum_pair[6] += 1.0;
        }
        else if (maleT[i][j] == 2 && maleP[i][j] == 2)
        {
            if (t == 2 && p == 2)
                sum_pair[7] += 1.0;
            else if (t == 3 && p == 3)
                sum_pair[8] += 1.0;
        }
        else if (maleT[i][j] == 3 && maleP[i][j] == 3)
        {
            if (t == 3 && p == 3)
                sum_pair[9] += 1.0;
        }
    }
}

int main(void)
{
    // 新しく挿入
    //  int tend;
    int **maleT, *base_maleT;
    int **maleP, *base_maleP;
    int **femaleT, *base_femaleT;
    int **femaleP, *base_femaleP;
    int **maleTdummy, *base_maleTdummy;
    int **malePdummy, *base_malePdummy;
    int **femaleTdummy, *base_femaleTdummy;
    int **femalePdummy, *base_femalePdummy;
    double initT2, initP2, initT3, initP3, initT2P2;
    int k, k2, i, j, i2, j2, t, ok, x1, x2, a, b, n;
    int maleI, maleJ, femaleI, femaleJ;
    int numMT1, numMT2, numMT3, numMP1, numMP2, numMP3;
    int numFT1, numFT2, numFT3, numFP1, numFP2, numFP3;
    double sum[9], gsum[9], sum_pair[10];
    double rnd, rnd2, rnd3, rnd4, rnd5, sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9, y1, z1, y2, z2, init, init1;
    double gsum1, gsum2, gsum3, gsum4, gsum5, gsum6, gsum7, gsum8, gsum9;
    int iK, iV, il, ia2, ia1;
    double K, V, l, a1, a2;
    double x11, x21, x22, x31, x32, x33, x41, x51, x91, x44, x54, x94, x55, x95, x99;
    double nx11, nx21, nx22, nx31, nx32, nx33;
    int maleT0, maleP0, femaleT0, femaleP0, mgenotype, fgenotype, count;
    FILE *gp, *data1, *data2, *data3, *data4, *data5, *data6, *data7, *data8;
    FILE *snapshot1, *snapshot2, *snapshot3, *snapshot4, *snapshot5, *snapshot6;
    char *data_file1, *data_file2, *data_file3, *data_file4, *data_file5, *data_file6, *data_file7, *data_file8;
    char *snapshot_file1, *snapshot_file2, *snapshot_file3, *snapshot_file4, *snapshot_file5, *snapshot_file6;
    char *Figaxis[14] = {"x_11", "x_41", "x_51", "x_81", "x_44", "x_54", "x_84", "x_55", "x_85", "x_88", "T1P1", "T2P1", "T2P2", "T3P3"};
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);
    // data1(T2P2頻度図の書き込み用)
    typedef struct
    {
        int t;
        double sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9, initT2P2;
    } RecordT2P2;

    // RecordT2P2 *buffer = malloc(sizeof(RecordT2P2) * 9*(tend + 1));
    // RecordT2P2 *buft3p3 = malloc(sizeof(RecordT2P2) * 9*(tend + 1));
    int buf_count = 0;
    // data7(遺伝子頻度の書き込み用)
    typedef struct
    {
        int t;
        double sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9;
    } Recordgenoport;

    // Recordgenoport *genorepo = malloc(sizeof(Recordgenoport) * 9*(tend + 1));
    int geno_count = 0;

    // snapshot2(空間構造の書き込み用)
    typedef struct
    {
        int i, j, mgenotype, fgenotype;

    } Recordmap;

    Recordmap *recomap = malloc(sizeof(Recordmap) * LH * LV);

    // 局所的ペア頻度
    typedef struct
    {
        int t;
        double x11, x41, x51, x91, x44, x54, x94, x55, x95, x99, initT2P2;
    } Localpair;
    Localpair *Pair = malloc(sizeof(Localpair) * 4 * (tend + 1));
    int pair_count = 0;
    Localpair *Pair_freq = malloc(sizeof(Localpair) * 4 * (tend + 1));
    int freq_count = 0;

    maleT = malloc(sizeof(int *) * LH);
    maleP = malloc(sizeof(int *) * LH);
    femaleT = malloc(sizeof(int *) * LH);
    femaleP = malloc(sizeof(int *) * LH);
    maleTdummy = malloc(sizeof(int *) * LH);
    malePdummy = malloc(sizeof(int *) * LH);
    femaleTdummy = malloc(sizeof(int *) * LH);
    femalePdummy = malloc(sizeof(int *) * LH);
    base_maleT = malloc(sizeof(int) * LH * LV);
    base_maleP = malloc(sizeof(int) * LH * LV);
    base_femaleT = malloc(sizeof(int) * LH * LV);
    base_femaleP = malloc(sizeof(int) * LH * LV);
    base_maleTdummy = malloc(sizeof(int) * LH * LV);
    base_malePdummy = malloc(sizeof(int) * LH * LV);
    base_femaleTdummy = malloc(sizeof(int) * LH * LV);
    base_femalePdummy = malloc(sizeof(int) * LH * LV);
    for (i = 0; i < LH; i++)
    {
        maleT[i] = base_maleT + i * LV;
        maleP[i] = base_maleP + i * LV;
        femaleT[i] = base_femaleT + i * LV;
        femaleP[i] = base_femaleP + i * LV;
        maleTdummy[i] = base_maleTdummy + i * LV;
        malePdummy[i] = base_malePdummy + i * LV;
        femaleTdummy[i] = base_femaleTdummy + i * LV;
        femalePdummy[i] = base_femalePdummy + i * LV;
    }

    int num_threads = omp_get_num_procs(); // 最大利用可能スレッド数（論理コア数）
    mt_state *rng_states = malloc(sizeof(mt_state) * num_threads);
    if (!rng_states)
    {
        perror("malloc rng_states");
        return 1;
    }

    for (int tid = 0; tid < num_threads; tid++)
    {
        init_genrand_mt(&rng_states[tid], 5489UL + (unsigned long)tid * 12345UL);
    }

    omp_set_num_threads(num_threads);
    printf("Using %d threads\n", num_threads);
    fflush(stdout);
    RecordT2P2 *buffer = malloc(sizeof(RecordT2P2) * 9 * (tend + 1));
    RecordT2P2 *buft3p3 = malloc(sizeof(RecordT2P2) * 9 * (tend + 1));
    Recordgenoport *genorepo = malloc(sizeof(Recordgenoport) * 9 * (tend + 1));

    for (iK = 0; iK <= 11; iK++)
    { // iK=1;iK<=3;iK++
        if (iK == 0)
            K = 0.05;
        else if (iK == 1)
            K = 0.1;
        else if (iK == 2)
            K = 0.11;
        else if (iK == 3)
            K = 0.115;
        else if (iK == 4)
            K = 0.12;
        else if (iK == 5)
            K = 0.125;
        else if (iK == 6)
            K = 0.13;
        else if (iK == 7)
            K = 0.14;
        else if (iK == 8)
            K = 0.15;
        else if (iK == 9)
            K = 0.2;
        else if (iK == 10)
            K = 0.25;
        else if (iK == 11)
            K = 0.3;
        // K = 0.05 + (double)0.01 * (double)iK;
        // K=0.11+(double)iK*0.001;
        for (iV = 2; iV <= 2; iV++)
        {
            V = (double)(iV * 2 - 1) * K / 6.0;
            // V=0.20+(double)iV*0.01;
            for (il = 2; il <= 2; il++)
            {
                l = (double)(il * 2 - 1) * u / 6.0;

                for (ia1 = 3; ia1 <= 3; ia1++)
                {
                    // if(ia1==3)continue;
                    // else a1=(double)ia1;
                    a1 = (double)ia1;
                    for (ia2 = 3; ia2 <= 3; ia2++)
                    {
                        a2 = (double)ia2;
                        printf("K V l a2:%f %f %f %f\n", K, V, l, a2);

                        init_genrand(0);

                        // snapshot_file1 = malloc(100);
                        // sprintf(snapshot_file1, "Three_env_stmap_%f_initP_%g.dat", K, mapinitP);

                        // snapshot_file2 = malloc(100);
                        // // sprintf(snapshot2, "efcs_intmap_%f.dat", K);

                        // snapshot_file3 = malloc(100);
                        // sprintf(snapshot_file3, "Three_env_finmap_%f_initP_%g.dat", K, mapinitP);

                        data_file1 = malloc(100);
                        sprintf(data_file1, "Three_env_T3P3_K_%f_V_%f_l_%f_a1_%f_a2_%f.dat", K, V, l, a1, a2);

                        data_file2 = malloc(100);
                        sprintf(data_file2, "Three_env_T3P3_flow_K_%f_V_%f_l_%f_a1_%f_a2_%f.dat", K, V, l, a1, a2);

                        data_file4 = malloc(100);
                        sprintf(data_file4, "Three_env_pair_K_%f_V_%f_l_%f_a1_%f_a2_%f.dat", K, V, l, a1, a2);

                        data_file5 = malloc(100);
                        sprintf(data_file5, "Three_env_pair_flow_K_%f_V_%f_l_%f_a1_%f_a2_%f.dat", K, V, l, a1, a2);

                        data_file6 = malloc(100);
                        sprintf(data_file6, "Three_env_pair_final_K_%f_V_%f_l_%f_a1_%f_a2_%f.dat", K, V, l, a1, a2);

                        data_file7 = malloc(100);
                        sprintf(data_file7, "Three_env_genoport_K_%f_V_%f_l_%f_a1_%f_a2_%f.dat", K, V, l, a1, a2);

                        data_file8 = malloc(100);
                        sprintf(data_file8, "Three_env_pair_full_K_%f_V_%f_l_%f_a1_%f_a2_%f.dat", K, V, l, a1, a2);

                        snapshot_file2 = malloc(100);

                        buf_count = 0;
                        geno_count = 0;
                        pair_count = 0;
                        freq_count = 0;
                        for (k = 1; k <= 1; k++)
                        {
                            // initT2 = 0.1 * initialT;
                            // initT3 = (1.0 - initT2) / 2;
                            for (k2 = 1; k2 <= 4; k2++) // k2=1; k2<=9; k2++
                            {
                                // data1 = fopen(data_file1, "a");
                                // data4 = fopen(data_file4, "a");
                                // data7 = fopen(data_file7, "a");

                                // initP2 = 0.1 * k2;
                                // initP3 = (1.0 - initP2) / 2;
                                initT2P2 = 0.1 * 2.0 * (double)k2;
                                for (i = 0; i < LH; i++)
                                {
                                    for (j = 0; j < LV; j++)
                                    {
                                        // T2P2で分ける。他はランダム
                                        rnd = genrand_real2();
                                        if (rnd < initT2P2)
                                        {
                                            maleT[i][j] = 2;
                                            maleP[i][j] = 2;
                                        }
                                        else
                                        {
                                            rnd2 = genrand_real2();
                                            if (rnd2 < 0.125)
                                            {
                                                maleT[i][j] = 1;
                                                maleP[i][j] = 1;
                                            }
                                            else if (rnd2 < 0.25)
                                            {
                                                maleT[i][j] = 1;
                                                maleP[i][j] = 2;
                                            }
                                            else if (rnd2 < 0.375)
                                            {
                                                maleT[i][j] = 1;
                                                maleP[i][j] = 3;
                                            }
                                            else if (rnd2 < 0.5)
                                            {
                                                maleT[i][j] = 2;
                                                maleP[i][j] = 1;
                                            }
                                            else if (rnd2 < 0.625)
                                            {
                                                maleT[i][j] = 2;
                                                maleP[i][j] = 3;
                                            }
                                            else if (rnd2 < 0.75)
                                            {
                                                maleT[i][j] = 3;
                                                maleP[i][j] = 1;
                                            }
                                            else if (rnd2 < 0.875)
                                            {
                                                maleT[i][j] = 3;
                                                maleP[i][j] = 2;
                                            }
                                            else
                                            {
                                                maleT[i][j] = 3;
                                                maleP[i][j] = 3;
                                            }
                                        }
                                        rnd = genrand_real2();
                                        if (rnd < initT2P2)
                                        {
                                            femaleT[i][j] = 2;
                                            femaleP[i][j] = 2;
                                        }
                                        else
                                        {
                                            rnd2 = genrand_real2();
                                            if (rnd2 < 0.125)
                                            {
                                                femaleT[i][j] = 1;
                                                femaleP[i][j] = 1;
                                            }
                                            else if (rnd2 < 0.25)
                                            {
                                                femaleT[i][j] = 1;
                                                femaleP[i][j] = 2;
                                            }
                                            else if (rnd2 < 0.375)
                                            {
                                                femaleT[i][j] = 1;
                                                femaleP[i][j] = 3;
                                            }
                                            else if (rnd2 < 0.5)
                                            {
                                                femaleT[i][j] = 2;
                                                femaleP[i][j] = 1;
                                            }
                                            else if (rnd2 < 0.625)
                                            {
                                                femaleT[i][j] = 2;
                                                femaleP[i][j] = 3;
                                            }
                                            else if (rnd2 < 0.75)
                                            {
                                                femaleT[i][j] = 3;
                                                femaleP[i][j] = 1;
                                            }
                                            else if (rnd2 < 0.875)
                                            {
                                                femaleT[i][j] = 3;
                                                femaleP[i][j] = 2;
                                            }
                                            else
                                            {
                                                femaleT[i][j] = 3;
                                                femaleP[i][j] = 3;
                                            }
                                        }

                                        // T1P1vsT2P1

                                        // T2P1vsT2P2

                                        // T1P1vsT2P2
                                    }
                                }

                                // numMT1 = numMT2 = numMT3 = numFT1 = numFT2 = numFT3 = numMP1 = numMP2 = numMP3 = numFP1 = numFP2 = numFP3 = 0;
                                // for (i = 0; i < LH; i++)
                                // {
                                //     for (j = 0; j < LV; j++)
                                //     {
                                //         if (maleT[i][j] == 1)
                                //             numMT1++;
                                //         else if (maleT[i][j] == 2)
                                //             numMT2++;
                                //         else if (maleT[i][j] == 3)
                                //             numMT3++;
                                //         if (femaleT[i][j] == 1)
                                //             numFT1++;
                                //         else if (femaleT[i][j] == 2)
                                //             numFT2++;
                                //         else if (femaleT[i][j] == 3)
                                //             numFT3++;
                                //         if (maleP[i][j] == 1)
                                //             numMP1++;
                                //         else if (maleP[i][j] == 2)
                                //             numMP2++;
                                //         else if (maleP[i][j] == 3)
                                //             numMP3++;
                                //         if (femaleP[i][j] == 1)
                                //             numFP1++;
                                //         else if (femaleP[i][j] == 2)
                                //             numFP2++;
                                //         else if (femaleP[i][j] == 3)
                                //             numFP3++;
                                //     }
                                // }
                                // printf("0 T2 P2:%f %f \n", (double)numMT2 / (double)(LH * LV), (double)numMP2 / (double)(LH * LV));

                                // 最初の割合を出力
                                sum1 = sum2 = sum3 = sum4 = sum5 = sum6 = sum7 = sum8 = sum9 = 0.0;
                                // data7=fopen(data_file7,"a");
                                for (i = 0; i < LH; i++)
                                {
                                    for (j = 0; j < LV; j++)
                                    {
                                        if ((maleT[i][j] == 1 && maleP[i][j] == 1))
                                            sum1 += 1;
                                        else if ((maleT[i][j] == 1 && maleP[i][j] == 2))
                                            sum2 += 1;
                                        else if ((maleT[i][j] == 1 && maleP[i][j] == 3))
                                            sum3 += 1;
                                        else if ((maleT[i][j] == 2 && maleP[i][j] == 1))
                                            sum4 += 1;
                                        else if ((maleT[i][j] == 2 && maleP[i][j] == 2))
                                            sum5 += 1;
                                        else if ((maleT[i][j] == 2 && maleP[i][j] == 3))
                                            sum6 += 1;
                                        else if ((maleT[i][j] == 3 && maleP[i][j] == 1))
                                            sum7 += 1;
                                        else if ((maleT[i][j] == 3 && maleP[i][j] == 2))
                                            sum8 += 1;
                                        else if ((maleT[i][j] == 3 && maleP[i][j] == 3))
                                            sum9 += 1;
                                        calc_pair(sum_pair, i, j, maleT, maleP);
                                    }
                                }
                                // fprintf(data7, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 0,(double)sum1/(double)(LH*LV),(double)sum2/(double)(LH*LV),\
                                // (double)sum3/(double)(LH*LV),(double)sum4/(double)(LH*LV),(double)sum5/(double)(LH*LV),(double)sum6/(double)(LH*LV),\
                                // (double)sum7/(double)(LH*LV),(double)sum8/(double)(LH*LV),(double)sum9/(double)(LH*LV));
                                // fclose(data7);
                                genorepo[geno_count].t = 0;
                                genorepo[geno_count].sum1 = (double)sum1 / (double)(LH * LV);
                                genorepo[geno_count].sum2 = (double)sum2 / (double)(LH * LV);
                                genorepo[geno_count].sum3 = (double)sum3 / (double)(LH * LV);
                                genorepo[geno_count].sum4 = (double)sum4 / (double)(LH * LV);
                                genorepo[geno_count].sum5 = (double)sum5 / (double)(LH * LV);
                                genorepo[geno_count].sum6 = (double)sum6 / (double)(LH * LV);
                                genorepo[geno_count].sum7 = (double)sum7 / (double)(LH * LV);
                                genorepo[geno_count].sum8 = (double)sum8 / (double)(LH * LV);
                                genorepo[geno_count].sum9 = (double)sum9 / (double)(LH * LV);
                                geno_count++;

                                buffer[buf_count].t = 0;
                                buffer[buf_count].sum1 = (double)sum1 / (double)(LH * LV);
                                buffer[buf_count].sum2 = (double)sum2 / (double)(LH * LV);
                                buffer[buf_count].sum3 = (double)sum3 / (double)(LH * LV);
                                buffer[buf_count].sum4 = (double)sum4 / (double)(LH * LV);
                                buffer[buf_count].sum5 = (double)sum5 / (double)(LH * LV);
                                buffer[buf_count].sum6 = (double)sum6 / (double)(LH * LV);
                                buffer[buf_count].sum7 = (double)sum7 / (double)(LH * LV);
                                buffer[buf_count].sum8 = (double)sum8 / (double)(LH * LV);
                                buffer[buf_count].sum9 = (double)sum9 / (double)(LH * LV);
                                buffer[buf_count].initT2P2 = initT2P2;
                                buf_count++;

                                Pair[pair_count].t = 0;
                                Pair[pair_count].x11 = sum_pair[0] / ((double)4.0 * sum1);
                                Pair[pair_count].x41 = sum_pair[1] / ((double)4.0 * sum1);
                                Pair[pair_count].x51 = sum_pair[2] / ((double)4.0 * sum1);
                                Pair[pair_count].x91 = sum_pair[3] / ((double)4.0 * sum1);
                                Pair[pair_count].x44 = sum_pair[4] / ((double)4.0 * sum4);
                                Pair[pair_count].x54 = sum_pair[5] / ((double)4.0 * sum4);
                                Pair[pair_count].x94 = sum_pair[6] / ((double)4.0 * sum4);
                                Pair[pair_count].x55 = sum_pair[7] / ((double)4.0 * sum5);
                                Pair[pair_count].x95 = sum_pair[8] / ((double)4.0 * sum5);
                                Pair[pair_count].x99 = sum_pair[9] / ((double)4.0 * sum9);
                                Pair[pair_count].initT2P2 = initT2P2;
                                pair_count++;

                                Pair_freq[freq_count].t = 0;
                                Pair_freq[freq_count].x11 = sum_pair[0] / ((double)4.0 * sum1);
                                Pair_freq[freq_count].x41 = sum_pair[1] / ((double)4.0 * sum1);
                                Pair_freq[freq_count].x51 = sum_pair[2] / ((double)4.0 * sum1);
                                Pair_freq[freq_count].x91 = sum_pair[3] / ((double)4.0 * sum1);
                                Pair_freq[freq_count].x44 = sum_pair[4] / ((double)4.0 * sum4);
                                Pair_freq[freq_count].x54 = sum_pair[5] / ((double)4.0 * sum4);
                                Pair_freq[freq_count].x94 = sum_pair[6] / ((double)4.0 * sum4);
                                Pair_freq[freq_count].x55 = sum_pair[7] / ((double)4.0 * sum5);
                                Pair_freq[freq_count].x95 = sum_pair[8] / ((double)4.0 * sum5);
                                Pair_freq[freq_count].x99 = sum_pair[9] / ((double)4.0 * sum9);
                                Pair_freq[freq_count].initT2P2 = initT2P2;
                                freq_count++;

                                // --- 追加: dummy 配列を初期化（未書き込み領域を防ぐ） ---
#pragma omp parallel for collapse(2) schedule(static) default(none) shared(maleT, maleP, femaleT, femaleP,                         \
                                                                               maleTdummy, malePdummy, femaleTdummy, femalePdummy) \
    private(i, j)
                                for (i = 0; i < LH; i++)
                                {
                                    for (j = 0; j < LV; j++)
                                    {
                                        maleTdummy[i][j] = maleT[i][j];
                                        malePdummy[i][j] = maleP[i][j];
                                        femaleTdummy[i][j] = femaleT[i][j];
                                        femalePdummy[i][j] = femaleP[i][j];
                                    }
                                }

                                t = 0;
                                for (t = 1; t <= tend; t++)
                                {

#pragma omp parallel default(none) shared(maleT, maleP, femaleT, femaleP, maleTdummy, malePdummy, femaleTdummy, femalePdummy, rng_states, K, V, l, a1, a2) \
    private(i, j, sum, gsum, rnd2, rnd4, rnd5, maleT0, maleP0, femaleT0, femaleP0)
                                    {
                                        int tid = omp_get_thread_num();
                                        mt_state *rng = &rng_states[tid];
#pragma omp for collapse(2) schedule(static)
                                        for (i = 0; i < LH; i++)
                                            for (j = 0; j < LV; j++)
                                            {
                                                int tid = omp_get_thread_num();
                                                calc_male_sum(sum, i, j, maleT, maleP, femaleP[i][j], a1, a2);

                                                // メス遺伝のための重みづけ
                                                //  メスのコスト

                                                calc_female_sum(gsum, i, j, femaleT, femaleP);
                                                // printf("calc_female_sum OK");

                                                rnd4 = genrand_real2_mt(rng);
                                                // オス遺伝
                                                if (rnd4 < 0.5)
                                                {
                                                    rnd5 = genrand_real2_mt(rng);
                                                    if (rnd5 < 0.5)
                                                    {
                                                        do
                                                        {
                                                            rnd2 = genrand_real2_mt(rng);
                                                            genotype(sum, &maleT0, &maleP0, rng);
                                                            // printf("genotype OK");

                                                            if (maleT0 == 1)
                                                                break;
                                                            if (maleT0 == 2 && rnd2 < l)
                                                                continue;
                                                            if (maleT0 == 3 && rnd2 < u)
                                                                continue;
                                                            break;
                                                        } while (1); // 次世代がオスの場合

                                                        maleTdummy[i][j] = maleT0;
                                                        malePdummy[i][j] = maleP0;
                                                    }
                                                    else
                                                    {

                                                        do
                                                        {
                                                            rnd2 = genrand_real2_mt(rng);
                                                            genotype(sum, &maleT0, &maleP0, rng);
                                                            // printf("genotype OK");

                                                            if (maleP0 == 1)
                                                                break;
                                                            if (maleP0 == 2 && rnd2 < V)
                                                                continue;
                                                            if (maleP0 == 3 && rnd2 < K)
                                                                continue;
                                                            break;

                                                        } while (1);

                                                        femaleTdummy[i][j] = maleT0;
                                                        femalePdummy[i][j] = maleP0;
                                                    }
                                                }
                                                else
                                                { // メス遺伝
                                                    rnd5 = genrand_real2_mt(rng);
                                                    if (rnd5 < 0.5)
                                                    { // 次世代がオス
                                                        femaleT0 = femaleT[i][j];
                                                        femaleP0 = femaleP[i][j];
                                                        if (femaleT0 != 1)
                                                        {
                                                            rnd2 = genrand_real2_mt(rng);
                                                            if (femaleT0 == 2 && rnd2 < l || femaleT0 == 3 && rnd2 < u)
                                                            {
                                                                do
                                                                {
                                                                    rnd2 = genrand_real2_mt(rng);
                                                                    genotype(gsum, &femaleT0, &femaleP0, rng);
                                                                    // printf("genotype OK");

                                                                    if (femaleT0 == 1)
                                                                        break;
                                                                    if (femaleT0 == 2 && rnd2 < l)
                                                                        continue;
                                                                    if (femaleT0 == 3 && rnd2 < u)
                                                                        continue;
                                                                    break;
                                                                } while (1);
                                                            }
                                                        }
                                                        maleTdummy[i][j] = femaleT0;
                                                        malePdummy[i][j] = femaleP0;
                                                    }
                                                    else
                                                    { // 次世代がメス
                                                        femaleT0 = femaleT[i][j];
                                                        femaleP0 = femaleP[i][j];
                                                        if (femaleP0 != 1)
                                                        {
                                                            rnd2 = genrand_real2_mt(rng);
                                                            if (femaleP0 == 2 && rnd2 < V || femaleP0 == 3 && rnd2 < K)
                                                            {
                                                                do
                                                                {
                                                                    rnd2 = genrand_real2_mt(rng);
                                                                    genotype(gsum, &femaleT0, &femaleP0, rng);
                                                                    // printf("genotype OK");

                                                                    if (femaleP0 == 1)
                                                                        break;
                                                                    if (femaleP0 == 2 && rnd2 < V)
                                                                        continue;
                                                                    if (femaleP0 == 3 && rnd2 < K)
                                                                        continue;
                                                                    break;

                                                                } while (1);
                                                            }
                                                        }
                                                        femaleTdummy[i][j] = femaleT0;
                                                        femalePdummy[i][j] = femaleP0;
                                                    }
                                                }
                                            }
                                    }
                                    int **tmp;

                                    tmp = maleT;
                                    maleT = maleTdummy;
                                    maleTdummy = tmp;

                                    tmp = maleP;
                                    maleP = malePdummy;
                                    malePdummy = tmp;

                                    tmp = femaleT;
                                    femaleT = femaleTdummy;
                                    femaleTdummy = tmp;

                                    tmp = femaleP;
                                    femaleP = femalePdummy;
                                    femalePdummy = tmp;

                                    // 途中の図
                                    // 途中の図
                                    // if (t % 100 == 0 && fabs(initT2P2 - mapinitP) < 1e-12 && t <= 50000 && fabs(V - K / 2) < 1e-12)
                                    // {
                                    //     sprintf(snapshot_file2, "Three_intmap_K_%f_V_%f_initP_%g_t_%d.dat", K, V, initT2P2, t);

                                    //     mgenotype = fgenotype = 0;

                                    //     for (i = 0; i < LH; i++)
                                    //     {
                                    //         for (j = 0; j < LV; j++)
                                    //         {
                                    //             if (maleT[i][j] == 1 && maleP[i][j] == 1)
                                    //                 mgenotype = 1;
                                    //             else if (maleT[i][j] == 1 && maleP[i][j] == 2)
                                    //                 mgenotype = 2;
                                    //             else if (maleT[i][j] == 1 && maleP[i][j] == 3)
                                    //                 mgenotype = 3;
                                    //             else if (maleT[i][j] == 2 && maleP[i][j] == 1)
                                    //                 mgenotype = 4;
                                    //             else if (maleT[i][j] == 2 && maleP[i][j] == 2)
                                    //                 mgenotype = 5;
                                    //             else if (maleT[i][j] == 2 && maleP[i][j] == 3)
                                    //                 mgenotype = 6;
                                    //             else if (maleT[i][j] == 3 && maleP[i][j] == 1)
                                    //                 mgenotype = 7;
                                    //             else if (maleT[i][j] == 3 && maleP[i][j] == 2)
                                    //                 mgenotype = 8;
                                    //             else if (maleT[i][j] == 3 && maleP[i][j] == 3)
                                    //                 mgenotype = 9;
                                    //             if (femaleT[i][j] == 1 && femaleP[i][j] == 1)
                                    //                 fgenotype = 1;
                                    //             else if (femaleT[i][j] == 1 && femaleP[i][j] == 2)
                                    //                 fgenotype = 2;
                                    //             else if (femaleT[i][j] == 1 && femaleP[i][j] == 3)
                                    //                 fgenotype = 3;
                                    //             else if (femaleT[i][j] == 2 && femaleP[i][j] == 1)
                                    //                 fgenotype = 4;
                                    //             else if (femaleT[i][j] == 2 && femaleP[i][j] == 2)
                                    //                 fgenotype = 5;
                                    //             else if (femaleT[i][j] == 2 && femaleP[i][j] == 3)
                                    //                 fgenotype = 6;
                                    //             else if (femaleT[i][j] == 3 && femaleP[i][j] == 1)
                                    //                 fgenotype = 7;
                                    //             else if (femaleT[i][j] == 3 && femaleP[i][j] == 2)
                                    //                 fgenotype = 8;
                                    //             else if (femaleT[i][j] == 3 && femaleP[i][j] == 3)
                                    //                 fgenotype = 9;

                                    //             recomap[i * LH + j].i = i;
                                    //             recomap[i * LH + j].j = j;
                                    //             recomap[i * LH + j].mgenotype = mgenotype;
                                    //             recomap[i * LH + j].fgenotype = fgenotype;
                                    //         }
                                    //     }
                                    //     snapshot2 = fopen(snapshot_file2, "w");
                                    //     for (n = 0; n < LH * LV; n++)
                                    //     {
                                    //         fprintf(snapshot2, "%d\t%d\t%d\t%d\n", recomap[n].i, recomap[n].j, recomap[n].mgenotype, recomap[n].fgenotype);
                                    //     }
                                    //     fclose(snapshot2);
                                    //     Map("male", snapshot_file2, K, V, mapinitP, t);
                                    //     Map("female", snapshot_file2, K, V, mapinitP, t);
                                    // }
                                    // 遺伝子型の割合出力
                                    sum1 = sum2 = sum3 = sum4 = sum5 = sum6 = sum7 = sum8 = sum9 = 0.0;
                                    // data7=fopen(data_file7,"a");
                                    for (i = 0; i < LH; i++)
                                    {
                                        for (j = 0; j < LV; j++)
                                        {
                                            if ((maleT[i][j] == 1 && maleP[i][j] == 1))
                                                sum1 += 1;
                                            else if ((maleT[i][j] == 1 && maleP[i][j] == 2))
                                                sum2 += 1;
                                            else if ((maleT[i][j] == 1 && maleP[i][j] == 3))
                                                sum3 += 1;
                                            else if ((maleT[i][j] == 2 && maleP[i][j] == 1))
                                                sum4 += 1;
                                            else if ((maleT[i][j] == 2 && maleP[i][j] == 2))
                                                sum5 += 1;
                                            else if ((maleT[i][j] == 2 && maleP[i][j] == 3))
                                                sum6 += 1;
                                            else if ((maleT[i][j] == 3 && maleP[i][j] == 1))
                                                sum7 += 1;
                                            else if ((maleT[i][j] == 3 && maleP[i][j] == 2))
                                                sum8 += 1;
                                            else if ((maleT[i][j] == 3 && maleP[i][j] == 3))
                                                sum9 += 1;
                                            calc_pair(sum_pair, i, j, maleT, maleP);
                                        }
                                    }
                                    genorepo[geno_count].t = t;
                                    genorepo[geno_count].sum1 = (double)sum1 / (double)(LH * LV);
                                    genorepo[geno_count].sum2 = (double)sum2 / (double)(LH * LV);
                                    genorepo[geno_count].sum3 = (double)sum3 / (double)(LH * LV);
                                    genorepo[geno_count].sum4 = (double)sum4 / (double)(LH * LV);
                                    genorepo[geno_count].sum5 = (double)sum5 / (double)(LH * LV);
                                    genorepo[geno_count].sum6 = (double)sum6 / (double)(LH * LV);
                                    genorepo[geno_count].sum7 = (double)sum7 / (double)(LH * LV);
                                    genorepo[geno_count].sum8 = (double)sum8 / (double)(LH * LV);
                                    genorepo[geno_count].sum9 = (double)sum9 / (double)(LH * LV);
                                    geno_count++;
                                    // fclose(data7);
                                    if (t % 100 == 0)
                                    {
                                        Pair[pair_count].t = t;
                                        if (sum1 == 0 && sum4 == 0 && sum5 == 0 && sum9 == 0)
                                        {
                                            Pair[pair_count].x11 = 0.0;
                                            Pair[pair_count].x41 = 0.0;
                                            Pair[pair_count].x51 = 0.0;
                                            Pair[pair_count].x91 = 0.0;
                                            Pair[pair_count].x44 = 0.0;
                                            Pair[pair_count].x54 = 0.0;
                                            Pair[pair_count].x94 = 0.0;
                                            Pair[pair_count].x55 = 0.0;
                                            Pair[pair_count].x95 = 0.0;
                                            Pair[pair_count].x99 = 0.0;
                                        }
                                        else if (sum1 == 0 && sum4 == 0 && sum5 == 0)
                                        {
                                            Pair[pair_count].x11 = 0.0;
                                            Pair[pair_count].x41 = 0.0;
                                            Pair[pair_count].x51 = 0.0;
                                            Pair[pair_count].x91 = 0.0;
                                            Pair[pair_count].x44 = 0.0;
                                            Pair[pair_count].x54 = 0.0;
                                            Pair[pair_count].x94 = 0.0;
                                            Pair[pair_count].x55 = 0.0;
                                            Pair[pair_count].x95 = 0.0;
                                            Pair[pair_count].x99 = sum_pair[9] / ((double)4.0 * sum9);
                                        }
                                        else if (sum1 == 0 && sum4 == 0 && sum9 == 0)
                                        {
                                            Pair[pair_count].x11 = 0.0;
                                            Pair[pair_count].x41 = 0.0;
                                            Pair[pair_count].x51 = 0.0;
                                            Pair[pair_count].x91 = 0.0;
                                            Pair[pair_count].x44 = 0.0;
                                            Pair[pair_count].x54 = 0.0;
                                            Pair[pair_count].x94 = 0.0;
                                            Pair[pair_count].x55 = sum_pair[7] / ((double)4.0 * sum5);
                                            Pair[pair_count].x95 = sum_pair[8] / ((double)4.0 * sum5);
                                            Pair[pair_count].x99 = 0.0;
                                        }
                                        else if (sum1 == 0 && sum5 == 0 && sum9 == 0)
                                        {
                                            Pair[pair_count].x11 = 0.0;
                                            Pair[pair_count].x41 = 0.0;
                                            Pair[pair_count].x51 = 0.0;
                                            Pair[pair_count].x91 = 0.0;
                                            Pair[pair_count].x44 = sum_pair[4] / ((double)4.0 * sum4);
                                            Pair[pair_count].x54 = sum_pair[5] / ((double)4.0 * sum4);
                                            Pair[pair_count].x94 = sum_pair[6] / ((double)4.0 * sum4);
                                            Pair[pair_count].x55 = 0.0;
                                            Pair[pair_count].x95 = 0.0;
                                            Pair[pair_count].x99 = 0.0;
                                        }
                                        else if (sum4 == 0 && sum5 == 0 && sum9 == 0)
                                        {
                                            Pair[pair_count].x11 = sum_pair[0] / ((double)4.0 * sum1);
                                            Pair[pair_count].x41 = sum_pair[1] / ((double)4.0 * sum1);
                                            Pair[pair_count].x51 = sum_pair[2] / ((double)4.0 * sum1);
                                            Pair[pair_count].x91 = sum_pair[3] / ((double)4.0 * sum1);
                                            Pair[pair_count].x44 = 0.0;
                                            Pair[pair_count].x54 = 0.0;
                                            Pair[pair_count].x94 = 0.0;
                                            Pair[pair_count].x55 = 0.0;
                                            Pair[pair_count].x95 = 0.0;
                                            Pair[pair_count].x99 = 0.0;
                                        }
                                        else if (sum1 == 0 && sum4 == 0)
                                        {
                                            Pair[pair_count].x11 = 0.0;
                                            Pair[pair_count].x41 = 0.0;
                                            Pair[pair_count].x51 = 0.0;
                                            Pair[pair_count].x91 = 0.0;
                                            Pair[pair_count].x44 = 0.0;
                                            Pair[pair_count].x54 = 0.0;
                                            Pair[pair_count].x94 = 0.0;
                                            Pair[pair_count].x55 = sum_pair[7] / ((double)4.0 * sum5);
                                            Pair[pair_count].x95 = sum_pair[8] / ((double)4.0 * sum5);
                                            Pair[pair_count].x99 = sum_pair[9] / ((double)4.0 * sum9);
                                        }
                                        else if (sum1 == 0 && sum5 == 0)
                                        {
                                            Pair[pair_count].x11 = 0.0;
                                            Pair[pair_count].x41 = 0.0;
                                            Pair[pair_count].x51 = 0.0;
                                            Pair[pair_count].x91 = 0.0;
                                            Pair[pair_count].x44 = sum_pair[4] / ((double)4.0 * sum4);
                                            Pair[pair_count].x54 = sum_pair[5] / ((double)4.0 * sum4);
                                            Pair[pair_count].x94 = sum_pair[6] / ((double)4.0 * sum4);
                                            Pair[pair_count].x55 = 0.0;
                                            Pair[pair_count].x95 = 0.0;
                                            Pair[pair_count].x99 = sum_pair[9] / ((double)4.0 * sum9);
                                        }
                                        else if (sum1 == 0 && sum9 == 0)
                                        {
                                            Pair[pair_count].x11 = 0.0;
                                            Pair[pair_count].x41 = 0.0;
                                            Pair[pair_count].x51 = 0.0;
                                            Pair[pair_count].x91 = 0.0;
                                            Pair[pair_count].x44 = sum_pair[4] / ((double)4.0 * sum4);
                                            Pair[pair_count].x54 = sum_pair[5] / ((double)4.0 * sum4);
                                            Pair[pair_count].x94 = sum_pair[6] / ((double)4.0 * sum4);
                                            Pair[pair_count].x55 = sum_pair[7] / ((double)4.0 * sum5);
                                            Pair[pair_count].x95 = sum_pair[8] / ((double)4.0 * sum5);
                                            Pair[pair_count].x99 = 0.0;
                                        }
                                        else if (sum4 == 0 && sum5 == 0)
                                        {
                                            Pair[pair_count].x11 = sum_pair[0] / ((double)4.0 * sum1);
                                            Pair[pair_count].x41 = sum_pair[1] / ((double)4.0 * sum1);
                                            Pair[pair_count].x51 = sum_pair[2] / ((double)4.0 * sum1);
                                            Pair[pair_count].x91 = sum_pair[3] / ((double)4.0 * sum1);
                                            Pair[pair_count].x44 = 0.0;
                                            Pair[pair_count].x54 = 0.0;
                                            Pair[pair_count].x94 = 0.0;
                                            Pair[pair_count].x55 = 0.0;
                                            Pair[pair_count].x95 = 0.0;
                                            Pair[pair_count].x99 = sum_pair[9] / ((double)4.0 * sum9);
                                        }
                                        else if (sum4 == 0 && sum9 == 0)
                                        {
                                            Pair[pair_count].x11 = sum_pair[0] / ((double)4.0 * sum1);
                                            Pair[pair_count].x41 = sum_pair[1] / ((double)4.0 * sum1);
                                            Pair[pair_count].x51 = sum_pair[2] / ((double)4.0 * sum1);
                                            Pair[pair_count].x91 = sum_pair[3] / ((double)4.0 * sum1);
                                            Pair[pair_count].x44 = 0.0;
                                            Pair[pair_count].x54 = 0.0;
                                            Pair[pair_count].x94 = 0.0;
                                            Pair[pair_count].x55 = sum_pair[7] / ((double)4.0 * sum5);
                                            Pair[pair_count].x95 = sum_pair[8] / ((double)4.0 * sum5);
                                            Pair[pair_count].x99 = 0.0;
                                        }
                                        else if (sum5 == 0 && sum9 == 0)
                                        {
                                            Pair[pair_count].x11 = sum_pair[0] / ((double)4.0 * sum1);
                                            Pair[pair_count].x41 = sum_pair[1] / ((double)4.0 * sum1);
                                            Pair[pair_count].x51 = sum_pair[2] / ((double)4.0 * sum1);
                                            Pair[pair_count].x91 = sum_pair[3] / ((double)4.0 * sum1);
                                            Pair[pair_count].x44 = sum_pair[4] / ((double)4.0 * sum4);
                                            Pair[pair_count].x54 = sum_pair[5] / ((double)4.0 * sum4);
                                            Pair[pair_count].x94 = sum_pair[6] / ((double)4.0 * sum4);
                                            Pair[pair_count].x55 = 0.0;
                                            Pair[pair_count].x95 = 0.0;
                                            Pair[pair_count].x99 = 0.0;
                                        }
                                        else if (sum1 == 0)
                                        {
                                            Pair[pair_count].x11 = 0.0;
                                            Pair[pair_count].x41 = 0.0;
                                            Pair[pair_count].x51 = 0.0;
                                            Pair[pair_count].x91 = 0.0;
                                            Pair[pair_count].x44 = sum_pair[4] / ((double)4.0 * sum4);
                                            Pair[pair_count].x54 = sum_pair[5] / ((double)4.0 * sum4);
                                            Pair[pair_count].x94 = sum_pair[6] / ((double)4.0 * sum4);
                                            Pair[pair_count].x55 = sum_pair[7] / ((double)4.0 * sum5);
                                            Pair[pair_count].x95 = sum_pair[8] / ((double)4.0 * sum5);
                                            Pair[pair_count].x99 = sum_pair[9] / ((double)4.0 * sum9);
                                        }
                                        else if (sum4 == 0)
                                        {
                                            Pair[pair_count].x11 = sum_pair[0] / ((double)4.0 * sum1);
                                            Pair[pair_count].x41 = sum_pair[1] / ((double)4.0 * sum1);
                                            Pair[pair_count].x51 = sum_pair[2] / ((double)4.0 * sum1);
                                            Pair[pair_count].x91 = sum_pair[3] / ((double)4.0 * sum1);
                                            Pair[pair_count].x44 = 0.0;
                                            Pair[pair_count].x54 = 0.0;
                                            Pair[pair_count].x94 = 0.0;
                                            Pair[pair_count].x55 = sum_pair[7] / ((double)4.0 * sum5);
                                            Pair[pair_count].x95 = sum_pair[8] / ((double)4.0 * sum5);
                                            Pair[pair_count].x99 = sum_pair[9] / ((double)4.0 * sum9);
                                        }
                                        else if (sum5 == 0)
                                        {
                                            Pair[pair_count].x11 = sum_pair[0] / ((double)4.0 * sum1);
                                            Pair[pair_count].x41 = sum_pair[1] / ((double)4.0 * sum1);
                                            Pair[pair_count].x51 = sum_pair[2] / ((double)4.0 * sum1);
                                            Pair[pair_count].x91 = sum_pair[3] / ((double)4.0 * sum1);
                                            Pair[pair_count].x44 = sum_pair[4] / ((double)4.0 * sum4);
                                            Pair[pair_count].x54 = sum_pair[5] / ((double)4.0 * sum4);
                                            Pair[pair_count].x94 = sum_pair[6] / ((double)4.0 * sum4);
                                            Pair[pair_count].x55 = 0.0;
                                            Pair[pair_count].x95 = 0.0;
                                            Pair[pair_count].x99 = sum_pair[9] / ((double)4.0 * sum9);
                                        }
                                        else if (sum9 == 0)
                                        {
                                            Pair[pair_count].x11 = sum_pair[0] / ((double)4.0 * sum1);
                                            Pair[pair_count].x41 = sum_pair[1] / ((double)4.0 * sum1);
                                            Pair[pair_count].x51 = sum_pair[2] / ((double)4.0 * sum1);
                                            Pair[pair_count].x91 = sum_pair[3] / ((double)4.0 * sum1);
                                            Pair[pair_count].x44 = sum_pair[4] / ((double)4.0 * sum4);
                                            Pair[pair_count].x54 = sum_pair[5] / ((double)4.0 * sum4);
                                            Pair[pair_count].x94 = sum_pair[6] / ((double)4.0 * sum4);
                                            Pair[pair_count].x55 = sum_pair[7] / ((double)4.0 * sum5);
                                            Pair[pair_count].x95 = sum_pair[8] / ((double)4.0 * sum5);
                                            Pair[pair_count].x99 = 0.0;
                                        }
                                        else
                                        {
                                            Pair[pair_count].x11 = sum_pair[0] / ((double)4.0 * sum1);
                                            Pair[pair_count].x41 = sum_pair[1] / ((double)4.0 * sum1);
                                            Pair[pair_count].x51 = sum_pair[2] / ((double)4.0 * sum1);
                                            Pair[pair_count].x91 = sum_pair[3] / ((double)4.0 * sum1);
                                            Pair[pair_count].x44 = sum_pair[4] / ((double)4.0 * sum4);
                                            Pair[pair_count].x54 = sum_pair[5] / ((double)4.0 * sum4);
                                            Pair[pair_count].x94 = sum_pair[6] / ((double)4.0 * sum4);
                                            Pair[pair_count].x55 = sum_pair[7] / ((double)4.0 * sum5);
                                            Pair[pair_count].x95 = sum_pair[8] / ((double)4.0 * sum5);
                                            Pair[pair_count].x99 = sum_pair[9] / ((double)4.0 * sum9);
                                        }

                                        Pair[pair_count].initT2P2 = initT2P2;
                                        pair_count++;

                                        buffer[buf_count].t = t;
                                        buffer[buf_count].sum1 = (double)sum1 / (double)(LH * LV);
                                        buffer[buf_count].sum2 = (double)sum2 / (double)(LH * LV);
                                        buffer[buf_count].sum3 = (double)sum3 / (double)(LH * LV);
                                        buffer[buf_count].sum4 = (double)sum4 / (double)(LH * LV);
                                        buffer[buf_count].sum5 = (double)sum5 / (double)(LH * LV);
                                        buffer[buf_count].sum6 = (double)sum6 / (double)(LH * LV);
                                        buffer[buf_count].sum7 = (double)sum7 / (double)(LH * LV);
                                        buffer[buf_count].sum8 = (double)sum8 / (double)(LH * LV);
                                        buffer[buf_count].sum9 = (double)sum9 / (double)(LH * LV);
                                        buffer[buf_count].initT2P2 = initT2P2;
                                        buf_count++;
                                    }
                                    Pair_freq[freq_count].t = t;
                                    if (sum1 == 0 && sum4 == 0 && sum5 == 0 && sum9 == 0)
                                    {
                                        Pair_freq[freq_count].x11 = 0.0;
                                        Pair_freq[freq_count].x41 = 0.0;
                                        Pair_freq[freq_count].x51 = 0.0;
                                        Pair_freq[freq_count].x91 = 0.0;
                                        Pair_freq[freq_count].x44 = 0.0;
                                        Pair_freq[freq_count].x54 = 0.0;
                                        Pair_freq[freq_count].x94 = 0.0;
                                        Pair_freq[freq_count].x55 = 0.0;
                                        Pair_freq[freq_count].x95 = 0.0;
                                        Pair_freq[freq_count].x99 = 0.0;
                                    }
                                    else if (sum1 == 0 && sum4 == 0 && sum5 == 0)
                                    {
                                        Pair_freq[freq_count].x11 = 0.0;
                                        Pair_freq[freq_count].x41 = 0.0;
                                        Pair_freq[freq_count].x51 = 0.0;
                                        Pair_freq[freq_count].x91 = 0.0;
                                        Pair_freq[freq_count].x44 = 0.0;
                                        Pair_freq[freq_count].x54 = 0.0;
                                        Pair_freq[freq_count].x94 = 0.0;
                                        Pair_freq[freq_count].x55 = 0.0;
                                        Pair_freq[freq_count].x95 = 0.0;
                                        Pair_freq[freq_count].x99 = sum_pair[9] / ((double)4.0 * sum9);
                                    }
                                    else if (sum1 == 0 && sum4 == 0 && sum9 == 0)
                                    {
                                        Pair_freq[freq_count].x11 = 0.0;
                                        Pair_freq[freq_count].x41 = 0.0;
                                        Pair_freq[freq_count].x51 = 0.0;
                                        Pair_freq[freq_count].x91 = 0.0;
                                        Pair_freq[freq_count].x44 = 0.0;
                                        Pair_freq[freq_count].x54 = 0.0;
                                        Pair_freq[freq_count].x94 = 0.0;
                                        Pair_freq[freq_count].x55 = sum_pair[7] / ((double)4.0 * sum5);
                                        Pair_freq[freq_count].x95 = sum_pair[8] / ((double)4.0 * sum5);
                                        Pair_freq[freq_count].x99 = 0.0;
                                    }
                                    else if (sum1 == 0 && sum5 == 0 && sum9 == 0)
                                    {
                                        Pair_freq[freq_count].x11 = 0.0;
                                        Pair_freq[freq_count].x41 = 0.0;
                                        Pair_freq[freq_count].x51 = 0.0;
                                        Pair_freq[freq_count].x91 = 0.0;
                                        Pair_freq[freq_count].x44 = sum_pair[4] / ((double)4.0 * sum4);
                                        Pair_freq[freq_count].x54 = sum_pair[5] / ((double)4.0 * sum4);
                                        Pair_freq[freq_count].x94 = sum_pair[6] / ((double)4.0 * sum4);
                                        Pair_freq[freq_count].x55 = 0.0;
                                        Pair_freq[freq_count].x95 = 0.0;
                                        Pair_freq[freq_count].x99 = 0.0;
                                    }
                                    else if (sum4 == 0 && sum5 == 0 && sum9 == 0)
                                    {
                                        Pair_freq[freq_count].x11 = sum_pair[0] / ((double)4.0 * sum1);
                                        Pair_freq[freq_count].x41 = sum_pair[1] / ((double)4.0 * sum1);
                                        Pair_freq[freq_count].x51 = sum_pair[2] / ((double)4.0 * sum1);
                                        Pair_freq[freq_count].x91 = sum_pair[3] / ((double)4.0 * sum1);
                                        Pair_freq[freq_count].x44 = 0.0;
                                        Pair_freq[freq_count].x54 = 0.0;
                                        Pair_freq[freq_count].x94 = 0.0;
                                        Pair_freq[freq_count].x55 = 0.0;
                                        Pair_freq[freq_count].x95 = 0.0;
                                        Pair_freq[freq_count].x99 = 0.0;
                                    }
                                    else if (sum1 == 0 && sum4 == 0)
                                    {
                                        Pair_freq[freq_count].x11 = 0.0;
                                        Pair_freq[freq_count].x41 = 0.0;
                                        Pair_freq[freq_count].x51 = 0.0;
                                        Pair_freq[freq_count].x91 = 0.0;
                                        Pair_freq[freq_count].x44 = 0.0;
                                        Pair_freq[freq_count].x54 = 0.0;
                                        Pair_freq[freq_count].x94 = 0.0;
                                        Pair_freq[freq_count].x55 = sum_pair[7] / ((double)4.0 * sum5);
                                        Pair_freq[freq_count].x95 = sum_pair[8] / ((double)4.0 * sum5);
                                        Pair_freq[freq_count].x99 = sum_pair[9] / ((double)4.0 * sum9);
                                    }
                                    else if (sum1 == 0 && sum5 == 0)
                                    {
                                        Pair_freq[freq_count].x11 = 0.0;
                                        Pair_freq[freq_count].x41 = 0.0;
                                        Pair_freq[freq_count].x51 = 0.0;
                                        Pair_freq[freq_count].x91 = 0.0;
                                        Pair_freq[freq_count].x44 = sum_pair[4] / ((double)4.0 * sum4);
                                        Pair_freq[freq_count].x54 = sum_pair[5] / ((double)4.0 * sum4);
                                        Pair_freq[freq_count].x94 = sum_pair[6] / ((double)4.0 * sum4);
                                        Pair_freq[freq_count].x55 = 0.0;
                                        Pair_freq[freq_count].x95 = 0.0;
                                        Pair_freq[freq_count].x99 = sum_pair[9] / ((double)4.0 * sum9);
                                    }
                                    else if (sum1 == 0 && sum9 == 0)
                                    {
                                        Pair_freq[freq_count].x11 = 0.0;
                                        Pair_freq[freq_count].x41 = 0.0;
                                        Pair_freq[freq_count].x51 = 0.0;
                                        Pair_freq[freq_count].x91 = 0.0;
                                        Pair_freq[freq_count].x44 = sum_pair[4] / ((double)4.0 * sum4);
                                        Pair_freq[freq_count].x54 = sum_pair[5] / ((double)4.0 * sum4);
                                        Pair_freq[freq_count].x94 = sum_pair[6] / ((double)4.0 * sum4);
                                        Pair_freq[freq_count].x55 = sum_pair[7] / ((double)4.0 * sum5);
                                        Pair_freq[freq_count].x95 = sum_pair[8] / ((double)4.0 * sum5);
                                        Pair_freq[freq_count].x99 = 0.0;
                                    }
                                    else if (sum4 == 0 && sum5 == 0)
                                    {
                                        Pair_freq[freq_count].x11 = sum_pair[0] / ((double)4.0 * sum1);
                                        Pair_freq[freq_count].x41 = sum_pair[1] / ((double)4.0 * sum1);
                                        Pair_freq[freq_count].x51 = sum_pair[2] / ((double)4.0 * sum1);
                                        Pair_freq[freq_count].x91 = sum_pair[3] / ((double)4.0 * sum1);
                                        Pair_freq[freq_count].x44 = 0.0;
                                        Pair_freq[freq_count].x54 = 0.0;
                                        Pair_freq[freq_count].x94 = 0.0;
                                        Pair_freq[freq_count].x55 = 0.0;
                                        Pair_freq[freq_count].x95 = 0.0;
                                        Pair_freq[freq_count].x99 = sum_pair[9] / ((double)4.0 * sum9);
                                    }
                                    else if (sum4 == 0 && sum9 == 0)
                                    {
                                        Pair_freq[freq_count].x11 = sum_pair[0] / ((double)4.0 * sum1);
                                        Pair_freq[freq_count].x41 = sum_pair[1] / ((double)4.0 * sum1);
                                        Pair_freq[freq_count].x51 = sum_pair[2] / ((double)4.0 * sum1);
                                        Pair_freq[freq_count].x91 = sum_pair[3] / ((double)4.0 * sum1);
                                        Pair_freq[freq_count].x44 = 0.0;
                                        Pair_freq[freq_count].x54 = 0.0;
                                        Pair_freq[freq_count].x94 = 0.0;
                                        Pair_freq[freq_count].x55 = sum_pair[7] / ((double)4.0 * sum5);
                                        Pair_freq[freq_count].x95 = sum_pair[8] / ((double)4.0 * sum5);
                                        Pair_freq[freq_count].x99 = 0.0;
                                    }
                                    else if (sum5 == 0 && sum9 == 0)
                                    {
                                        Pair_freq[freq_count].x11 = sum_pair[0] / ((double)4.0 * sum1);
                                        Pair_freq[freq_count].x41 = sum_pair[1] / ((double)4.0 * sum1);
                                        Pair_freq[freq_count].x51 = sum_pair[2] / ((double)4.0 * sum1);
                                        Pair_freq[freq_count].x91 = sum_pair[3] / ((double)4.0 * sum1);
                                        Pair_freq[freq_count].x44 = sum_pair[4] / ((double)4.0 * sum4);
                                        Pair_freq[freq_count].x54 = sum_pair[5] / ((double)4.0 * sum4);
                                        Pair_freq[freq_count].x94 = sum_pair[6] / ((double)4.0 * sum4);
                                        Pair_freq[freq_count].x55 = 0.0;
                                        Pair_freq[freq_count].x95 = 0.0;
                                        Pair_freq[freq_count].x99 = 0.0;
                                    }
                                    else if (sum1 == 0)
                                    {
                                        Pair_freq[freq_count].x11 = 0.0;
                                        Pair_freq[freq_count].x41 = 0.0;
                                        Pair_freq[freq_count].x51 = 0.0;
                                        Pair_freq[freq_count].x91 = 0.0;
                                        Pair_freq[freq_count].x44 = sum_pair[4] / ((double)4.0 * sum4);
                                        Pair_freq[freq_count].x54 = sum_pair[5] / ((double)4.0 * sum4);
                                        Pair_freq[freq_count].x94 = sum_pair[6] / ((double)4.0 * sum4);
                                        Pair_freq[freq_count].x55 = sum_pair[7] / ((double)4.0 * sum5);
                                        Pair_freq[freq_count].x95 = sum_pair[8] / ((double)4.0 * sum5);
                                        Pair_freq[freq_count].x99 = sum_pair[9] / ((double)4.0 * sum9);
                                    }
                                    else if (sum4 == 0)
                                    {
                                        Pair_freq[freq_count].x11 = sum_pair[0] / ((double)4.0 * sum1);
                                        Pair_freq[freq_count].x41 = sum_pair[1] / ((double)4.0 * sum1);
                                        Pair_freq[freq_count].x51 = sum_pair[2] / ((double)4.0 * sum1);
                                        Pair_freq[freq_count].x91 = sum_pair[3] / ((double)4.0 * sum1);
                                        Pair_freq[freq_count].x44 = 0.0;
                                        Pair_freq[freq_count].x54 = 0.0;
                                        Pair_freq[freq_count].x94 = 0.0;
                                        Pair_freq[freq_count].x55 = sum_pair[7] / ((double)4.0 * sum5);
                                        Pair_freq[freq_count].x95 = sum_pair[8] / ((double)4.0 * sum5);
                                        Pair_freq[freq_count].x99 = sum_pair[9] / ((double)4.0 * sum9);
                                    }
                                    else if (sum5 == 0)
                                    {
                                        Pair_freq[freq_count].x11 = sum_pair[0] / ((double)4.0 * sum1);
                                        Pair_freq[freq_count].x41 = sum_pair[1] / ((double)4.0 * sum1);
                                        Pair_freq[freq_count].x51 = sum_pair[2] / ((double)4.0 * sum1);
                                        Pair_freq[freq_count].x91 = sum_pair[3] / ((double)4.0 * sum1);
                                        Pair_freq[freq_count].x44 = sum_pair[4] / ((double)4.0 * sum4);
                                        Pair_freq[freq_count].x54 = sum_pair[5] / ((double)4.0 * sum4);
                                        Pair_freq[freq_count].x94 = sum_pair[6] / ((double)4.0 * sum4);
                                        Pair_freq[freq_count].x55 = 0.0;
                                        Pair_freq[freq_count].x95 = 0.0;
                                        Pair_freq[freq_count].x99 = sum_pair[9] / ((double)4.0 * sum9);
                                    }
                                    else if (sum9 == 0)
                                    {
                                        Pair_freq[freq_count].x11 = sum_pair[0] / ((double)4.0 * sum1);
                                        Pair_freq[freq_count].x41 = sum_pair[1] / ((double)4.0 * sum1);
                                        Pair_freq[freq_count].x51 = sum_pair[2] / ((double)4.0 * sum1);
                                        Pair_freq[freq_count].x91 = sum_pair[3] / ((double)4.0 * sum1);
                                        Pair_freq[freq_count].x44 = sum_pair[4] / ((double)4.0 * sum4);
                                        Pair_freq[freq_count].x54 = sum_pair[5] / ((double)4.0 * sum4);
                                        Pair_freq[freq_count].x94 = sum_pair[6] / ((double)4.0 * sum4);
                                        Pair_freq[freq_count].x55 = sum_pair[7] / ((double)4.0 * sum5);
                                        Pair_freq[freq_count].x95 = sum_pair[8] / ((double)4.0 * sum5);
                                        Pair_freq[freq_count].x99 = 0.0;
                                    }
                                    else
                                    {
                                        Pair_freq[freq_count].x11 = sum_pair[0] / ((double)4.0 * sum1);
                                        Pair_freq[freq_count].x41 = sum_pair[1] / ((double)4.0 * sum1);
                                        Pair_freq[freq_count].x51 = sum_pair[2] / ((double)4.0 * sum1);
                                        Pair_freq[freq_count].x91 = sum_pair[3] / ((double)4.0 * sum1);
                                        Pair_freq[freq_count].x44 = sum_pair[4] / ((double)4.0 * sum4);
                                        Pair_freq[freq_count].x54 = sum_pair[5] / ((double)4.0 * sum4);
                                        Pair_freq[freq_count].x94 = sum_pair[6] / ((double)4.0 * sum4);
                                        Pair_freq[freq_count].x55 = sum_pair[7] / ((double)4.0 * sum5);
                                        Pair_freq[freq_count].x95 = sum_pair[8] / ((double)4.0 * sum5);
                                        Pair_freq[freq_count].x99 = sum_pair[9] / ((double)4.0 * sum9);
                                    }
                                    Pair_freq[freq_count].initT2P2 = initT2P2;
                                    freq_count++;
                                    for (i = 0; i < 10; i++)
                                        sum_pair[i] = 0.0;
                                }
                            }
                        }
                        data1 = fopen(data_file1, "w");
                        for (int n = 0; n < buf_count; n++)
                        {
                            fprintf(data1, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
                                    buffer[n].t, buffer[n].sum1, buffer[n].sum2, buffer[n].sum3, buffer[n].sum4,
                                    buffer[n].sum5, buffer[n].sum6, buffer[n].sum7, buffer[n].sum8, buffer[n].sum9,
                                    buffer[n].initT2P2);
                        }
                        fclose(data1);

                        // data4 = fopen(data_file4, "w");
                        // for (int n =0; n < buf_count; n++) {
                        //     fprintf(data4, "%d\t%f\t%f\t%f\n",
                        //     buffer[n].t, buffer[n].mt2, buffer[n].mp2, buffer[n].initP2);
                        // }
                        // fclose(data4);

                        data7 = fopen(data_file7, "w");
                        for (int n = 0; n < geno_count; n++)
                        {
                            fprintf(data7, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
                                    genorepo[n].t, genorepo[n].sum1, genorepo[n].sum2, genorepo[n].sum3, genorepo[n].sum4, genorepo[n].sum5,
                                    genorepo[n].sum6, genorepo[n].sum7, genorepo[n].sum8, genorepo[n].sum9);
                        }
                        fclose(data7);

                        data4 = fopen(data_file4, "w");
                        for (int n = 0; n < pair_count; n++)
                        {
                            fprintf(data4, "%d\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\n",
                                    Pair[n].t, Pair[n].x11, Pair[n].x41, Pair[n].x51, Pair[n].x91,
                                    Pair[n].x44, Pair[n].x54, Pair[n].x94, Pair[n].x55, Pair[n].x95, Pair[n].x99,
                                    buffer[n].sum1, buffer[n].sum4, buffer[n].sum5, buffer[n].sum9, Pair[n].initT2P2);
                        }
                        fclose(data4);

                        data8 = fopen(data_file8, "w");
                        for (int n = 0; n < freq_count; n++)
                        {
                            fprintf(data8, "%d\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\n",
                                    Pair_freq[n].t, Pair_freq[n].x11, Pair_freq[n].x41, Pair_freq[n].x51, Pair_freq[n].x91,
                                    Pair_freq[n].x44, Pair_freq[n].x54, Pair_freq[n].x94, Pair_freq[n].x55, Pair_freq[n].x95, Pair_freq[n].x99,
                                    buffer[n].sum1, buffer[n].sum4, buffer[n].sum5, buffer[n].sum9, Pair_freq[n].initT2P2);
                        }
                        fclose(data8);

                        printf("Ok\n");
                        // T2P2-T3P3図
                        data_file3 = malloc(100);
                        sprintf(data_file3, "Three_env_final_K_%f_V_%f_l_%f_a1_%f_a2_%f.dat", K, V, l, a1, a2);
                        // data_file3 = "final_env.dat";
                        gp = fopen(data_file1, "r");
                        data3 = fopen(data_file3, "w");
                        while (fscanf(gp, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &x1, &sum1, &sum2, &sum3, &sum4, &sum5,
                                      &sum6, &sum7, &sum8, &sum9, &init) == 11)
                        {
                            if (x1 == (tend))
                            {
                                fprintf(data3, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
                                        x1, sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9);
                            }
                        }
                        fclose(data3);
                        fclose(gp);

                        data2 = fopen(data_file2, "w");
                        data1 = fopen(data_file1, "r");
                        if (fscanf(data1, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &x1, &sum1, &sum2, &sum3, &sum4, &sum5,
                                   &sum6, &sum7, &sum8, &sum9, &init) != 11)
                            return 1;
                        while (fscanf(data1, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &x2, &gsum1, &gsum2, &gsum3, &gsum4, &gsum5,
                                      &gsum6, &gsum7, &gsum8, &gsum9, &init1) == 11)
                        {
                            if (fabs(init - init1) < 1e-12)
                            {
                                fprintf(data2, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
                                        sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9,
                                        gsum1, gsum2, gsum3, gsum4, gsum5, gsum6, gsum7, gsum8, gsum9);
                            }
                            x1 = x2;
                            sum1 = gsum1;
                            sum2 = gsum2;
                            sum3 = gsum3;
                            sum4 = gsum4;
                            sum5 = gsum5;
                            sum6 = gsum6;
                            sum7 = gsum7;
                            sum8 = gsum8;
                            sum9 = gsum9;
                            init = init1;
                        }
                        fclose(data2);
                        fclose(data1);

                        // 図示
                        //  gp = popen("gnuplot -persist", "w");
                        //  fprintf(gp, "set terminal png\n");
                        //  fprintf(gp, "set term pngcairo size 1000,700\n");
                        //  fprintf(gp, "set output 'Genotype_Threealleles/Three_env_T2P2_T3P3_K_%f_V_%f_l_%f_a1_%f_a2_%f.png'\n", K,V,l,a1,a2);
                        //  fprintf(gp, "set xrange [0:%f]\n", 1.0);
                        //  fprintf(gp, "set xlabel 'T2P2'\n");
                        //  fprintf(gp, "set yrange [0:%f]\n", 1.0);
                        //  fprintf(gp, "set ylabel 'T3P3'\n");
                        //  fprintf(gp, "plot \'%s\' using 7:10 with points pointtype 7 lc rgb 'blue' title \"survivalrateK=%f\",\'%s\' using 1:2:($3-$1):($4-$2) with vectors head filled lc rgb 'blue',\'%s\' using 2:3 with points pointtype 7 lc rgb 'red' title \"finalarrival\"\n", data_file1, K, data_file2, data_file3);
                        //  pclose(gp);

                        // T2P1-T3P3図

                        // 図示
                        //  gp = popen("gnuplot -persist", "w");
                        //  fprintf(gp, "set terminal png\n");
                        //  fprintf(gp, "set term pngcairo size 1000,700\n");
                        //  fprintf(gp, "set output 'Genotype_Threealleles/Three_env_T2P1_T3P3_K_%f_V_%f_l_%f_a1_%f_a2_%f.png'\n", K,V,l,a1,a2);
                        //  fprintf(gp, "set xrange [0:%f]\n", 1.0);
                        //  fprintf(gp, "set xlabel 'T2P1'\n");
                        //  fprintf(gp, "set yrange [0:%f]\n", 1.0);
                        //  fprintf(gp, "set ylabel 'T3P3'\n");
                        //  fprintf(gp, "plot \'%s\' using 6:10 with points pointtype 7 lc rgb 'blue' title \"survivalrateK=%f\",\'%s\' using 1:2:($3-$1):($4-$2) with vectors head filled lc rgb 'blue',\'%s\' using 2:3 with points pointtype 7 lc rgb 'red' title \"finalarrival\"\n", data_file1, K, data_file2, data_file3);
                        //  pclose(gp);

                        // T2P1-T2P2図

                        // 図示
                        //  gp = popen("gnuplot -persist", "w");
                        //  fprintf(gp, "set terminal png\n");
                        //  fprintf(gp, "set term pngcairo size 1000,700\n");
                        //  fprintf(gp, "set output 'Genotype_Threealleles/Three_env_T2P1_T2P2_K_%f_V_%f_l_%f_a1_%f_a2_%f.png'\n", K,V,l,a1,a2);
                        //  fprintf(gp, "set xrange [0:%f]\n", 1.0);
                        //  fprintf(gp, "set xlabel 'T2P1'\n");
                        //  fprintf(gp, "set yrange [0:%f]\n", 1.0);
                        //  fprintf(gp, "set ylabel 'T2P2'\n");
                        //  fprintf(gp, "plot \'%s\' using 6:7 with points pointtype 7 lc rgb 'blue' title \"survivalrateK=%f\",\'%s\' using 1:2:($3-$1):($4-$2) with vectors head filled lc rgb 'blue',\'%s\' using 2:3 with points pointtype 7 lc rgb 'red' title \"finalarrival\"\n", data_file1, K, data_file2, data_file3);
                        //  pclose(gp);

                        // 時間変化図

                        for (i = 1; i <= 4; i++)
                        {
                            initT2P2 = (double)0.1 * 2.0 * i;
                            gp = popen("gnuplot -persist", "w");
                            fprintf(gp, "set terminal png\n");
                            fprintf(gp, "set output 'Genotype_Threealleles_genoport/Three_env_genoport_K_%f_V_%f_l_%f_a1_%f_a2_%f_initT2P2_%f.png'\n", K, V, l, a1, a2, initT2P2);
                            fprintf(gp, "set xrange [0:%d]\n", tend);
                            fprintf(gp, "set xlabel 't'\n");
                            fprintf(gp, "set yrange [0:%f]\n", 1.0);
                            fprintf(gp, "set ylabel 'genotype_frequency'\n");
                            fprintf(gp, "titles='T1P1 T1P2 T1P3 T2P1 T2P2 T2P3 T3P1 T3P2 T3P3'\n");
                            fprintf(gp, "set style line 1 lc rgb \"#A9A9A9\" lw 2\n");
                            fprintf(gp, "set style line 2 lc rgb \"#ADD8E6\" lw 2\n");
                            fprintf(gp, "set style line 3 lc rgb \"#0000FF\" lw 2\n");
                            fprintf(gp, "set style line 4 lc rgb \"#00FF00\" lw 2\n");
                            fprintf(gp, "set style line 5 lc rgb \"#FFFF00\" lw 2\n");
                            fprintf(gp, "set style line 6 lc rgb \"#FFA500\" lw 2\n");
                            fprintf(gp, "set style line 7 lc rgb \"#FF0000\" lw 2\n");
                            fprintf(gp, "set style line 8 lc rgb \"#FF69B4\" lw 2\n");
                            fprintf(gp, "set style line 9 lc rgb \"#000000\" lw 2\n");

                            // fprintf(gp, "plot \'%s\' using 2:3 with points pointtype 7 lc rgb 'blue' title \"survivalrateV=%f\",\'%s\' using 1:2:($3-$1):($4-$2) with vectors head filled lc rgb 'blue',\'%s\' using 2:3 with points pointtype 7 lc rgb 'red' title \"finalarrival\"\n", data_file4, V, data_file5, data_file6);
                            fprintf(gp, "plot for [j=2:10] \'%s\' every ::%d::%d using 1:j with lines ls (j-1) title word(titles, j-1)\n", data_file7, (i - 1) * (tend + 1), i * (tend + 1) - 1);
                            pclose(gp);
                        }

                        // pairのグラフ
                        gp = fopen(data_file4, "r");
                        data6 = fopen(data_file6, "w");
                        while (fscanf(gp, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                                      &x1, &x11, &x41, &x51, &x91, &x44, &x54, &x94, &x55, &x95, &x99, &sum1, &sum4, &sum5, &sum9, &init) == 16)
                        {
                            if (x1 == (tend))
                            {

                                fprintf(data6, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", x1, x11, x21, x31, x22, x32, x33, sum1, sum3, sum4);
                            }
                        }
                        fclose(data6);
                        fclose(gp);

                        // data4 = fopen(data_file4, "r");
                        // data5 = fopen(data_file5, "w");
                        // if (fscanf(data4, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &x1, &x11, &x21, &x31, &x22, &x32, &x33, &sum1, &sum3, &sum4, &init) != 11)
                        //     return 1;
                        // while (fscanf(data4, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &x2, &nx11, &nx21, &nx31, &nx22, &nx32, &nx33, &gsum1, &gsum3, &gsum4, &init1) == 11)
                        // {
                        //     if (fabs(init - init1) < 1e-9)
                        //     {
                        //         fprintf(data5, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
                        //                 x11, x21, x31, x22, x32, x33, sum1, sum3, sum4, nx11, nx21, nx31, nx22, nx32, nx33, gsum1, gsum3, gsum4);
                        //     }
                        //     x1 = x2;
                        //     x11 = nx11;
                        //     x21 = nx21;
                        //     x31 = nx31;
                        //     x22 = nx22;
                        //     x32 = nx32;
                        //     x33 = nx33;
                        //     sum1 = gsum1;

                        //     sum3 = gsum3;
                        //     sum4 = gsum4;
                        //     init = init1;
                        // }
                        // fclose(data4);
                        // fclose(data5);

                        // 12/14 2:50　図はmathemathicaでいいかなと思いdatファイルだけ出すわ
                        //  for (i = 0; i < 9; i++)
                        //  {
                        //      for (j = i + 1; j < 9; j++)
                        //      {
                        //          gp = popen("gnuplot -persist", "w");
                        //          fprintf(gp, "set terminal png\n");
                        //          fprintf(gp, "set term pngcairo size 1000,700\n");
                        //          fprintf(gp, "set output 'Genotype_Twoalleles_pair/K_%f_a1_%f_%s_%s.png'\n", K, a1, Figaxis[i], Figaxis[j]);
                        //          fprintf(gp, "set xrange [0:%f]\n", 1.0);
                        //          fprintf(gp, "set xlabel \'%s\'\n", Figaxis[i]);
                        //          fprintf(gp, "set yrange [0:%f]\n", 1.0);
                        //          fprintf(gp, "set ylabel \'%s\'\n", Figaxis[j]);
                        //          fprintf(gp, "plot \'%s\' using %d:%d with points pointtype 7 lc rgb 'blue' title \
                        //              \"survivalrateK=%f\",\
                        //              \'%s\' using %d:%d:($%d-$%d):($%d-$%d) with vectors head filled lc rgb 'blue',\
                        //              \'%s\' using %d:%d with points pointtype 7 lc rgb 'red' title \"finalarrival\"\n",
                        //                  data_file4, i + 2, j + 2,
                        //                  K,
                        //                  data_file5, i + 1, j + 1, i + 10, i + 1, j + 10, j + 1,
                        //                  data_file6, i + 2, j + 2);
                        //          pclose(gp);
                        //      }
                        //  }

                        for (i = 1; i <= 4; i++)
                        {
                            initT2P2 = (double)0.2 * i;
                            gp = popen("gnuplot -persist", "w");
                            fprintf(gp, "set terminal png\n");
                            fprintf(gp, "set output 'Genotype_Threealleles_pair/Three_env_pair_K_%f_a1_%f_initT2P2_%f.png'\n", K, a1, initT2P2);
                            fprintf(gp, "set xrange [0:%d]\n", tend);
                            fprintf(gp, "set xlabel 't'\n");
                            fprintf(gp, "set yrange [0:%f]\n", 1.0);
                            fprintf(gp, "set ylabel 'pair_frequency'\n");
                            fprintf(gp, "titles='x_1/1 x_2/1 x_3/1 x_4/1 x_5/1 x_9/1 x_4/4 x_5/4 x_9/4 x_5/5 x_9/5 x_9/9'\n");
                            // 1. 赤 (Red): 最も目立つ基本色
                            fprintf(gp, "set style line 1 lc rgb \"#FF0000\" lw 2\n");

                            // 2. 青 (Blue): 赤と対比する基本色
                            fprintf(gp, "set style line 2 lc rgb \"#0000FF\" lw 2\n");

                            // 3. 濃い緑 (Dark Green): 黄緑だと見にくいので、濃い緑を採用
                            fprintf(gp, "set style line 3 lc rgb \"#008000\" lw 2\n");

                            // 4. マゼンタ (Magenta): 紫より明るく、赤や青と区別しやすい
                            fprintf(gp, "set style line 4 lc rgb \"#FF00FF\" lw 2\n");

                            // 5. オレンジ (Orange-Red): 黄色の代わり。赤とは区別できる濃さ
                            fprintf(gp, "set style line 5 lc rgb \"#FF4500\" lw 2\n");

                            // 6. 黒 (Black): 最も収縮する色。全体を引き締める
                            fprintf(gp, "set style line 6 lc rgb \"#000000\" lw 2\n");

                            // 7. ティール/青緑 (Teal): 水色は見にくいので、濃い青緑
                            fprintf(gp, "set style line 7 lc rgb \"#008080\" lw 2\n");

                            // 8. 紫 (Purple): マゼンタや青とは違う、深い紫
                            fprintf(gp, "set style line 8 lc rgb \"#800080\" lw 2\n");

                            // 9. 茶色 (Saddle Brown): 黄色系だが暗いので白背景ではっきり見える
                            fprintf(gp, "set style line 9 lc rgb \"#8B4513\" lw 2\n");

                            // 10. 濃いグレー (Dark Gray): 黒とは区別できるが、線として認識できる濃さ
                            fprintf(gp, "set style line 10 lc rgb \"#555555\" lw 2\n");
                            // fprintf(gp, "set style line 1 lc rgb \"#0000FF\" lw 2\n");
                            // fprintf(gp, "set style line 2 lc rgb \"#00CC00\" lw 2\n");
                            // fprintf(gp, "set style line 3 lc rgb \"#FF8800\" lw 2\n");
                            // fprintf(gp, "set style line 4 lc rgb \"#FF0000\" lw 2\n");
                            // fprintf(gp, "set style line 5 lc rgb \"#FF00FF\" lw 2\n");
                            // fprintf(gp, "set style line 6 lc rgb \"#000000\" lw 2\n");
                            // fprintf(gp, "set style line 7 lc rgb \"#00FFFF\" lw 2\n");
                            // fprintf(gp, "set style line 8 lc rgb \"#FFFF00\" lw 2\n");
                            // fprintf(gp, "set style line 9 lc rgb \"#FF00FF\" lw 2\n");
                            // fprintf(gp, "set style line 10 lc rgb \"#FF8800\" lw 2\n");

                            fprintf(gp, "plot for [j=2:10] \'%s\' every ::%d::%d using 1:j with lines ls (j-1) title word(titles, j-1)\n", data_file8, (i - 1) * (tend + 1), i * (tend + 1) - 1);
                            pclose(gp);
                        }

                        // Map("male", snapshot_file1, initP2, 0);
                        // Map("female", snapshot_file1, initP2, 0);
                        // Map("male", snapshot_file3, initP2, 0);
                        // Map("female", snapshot_file3, initP2, 0);

                        free(data_file1);
                        free(data_file2);
                        free(data_file3);
                        free(data_file4);
                        free(data_file5);
                        free(data_file6);
                        free(data_file7);
                        free(snapshot_file2);
                        free(buffer);
                        free(buft3p3);
                        free(recomap);
                        free(genorepo);
                    }
                }
            }
        }
    }
    free(rng_states);
    free(base_femalePdummy);
    free(base_femaleTdummy);
    free(base_malePdummy);
    free(base_maleTdummy);
    free(base_femaleP);
    free(base_femaleT);
    free(base_maleP);
    free(base_maleT);
    free(femalePdummy);
    free(femaleTdummy);
    free(malePdummy);
    free(maleTdummy);
    free(femaleP);
    free(femaleT);
    free(maleP);
    free(maleT);

    clock_gettime(CLOCK_MONOTONIC, &end);
    double elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) * 1e-9;
    printf("simulation elapsed: %.3f s\n", elapsed);
}
