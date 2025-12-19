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

#define LH 1000// 1000
#define LV 1000// 1000
// #define K 0.075 //P3メスのコスト
// #define V 0.074 //P2メスのコスト(0<V<K)
// #define u 0.3 // T3オスのコスト0.3
// #define l 0.15  //T2オスのコスト(0<l<u)
#define a1 3.0 // P2メスがT2オスを選好する倍率3.0
// #define a2 6.0    // P3メスがT3オスを選好する倍率
#define tend 70000 // 4000 80000 10000 70000
#define mapinitP 0.25
#define initialP 3
#define initialT 1
#define SAVE_INTERVAL 100 // 100世代ごとに書き出し
#define MAX_SAVE ((tend / SAVE_INTERVAL) + 2)
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
    int k, k2, i, j, i2, j2, t, ok, x1, x2, a, b, n;
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
    clock_gettime(CLOCK_MONOTONIC, &start);

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

    // data1(T2P2頻度図の書き込み用)
    typedef struct
    {
        int t;
        double geno1, geno2, geno3, geno4, initT2P1;
    } RecordT2P2;

    RecordT2P2 *buffer = malloc(sizeof(RecordT2P2) * 9 * (tend + 1));
    int buf_count = 0;
    // data7(遺伝子頻度の書き込み用)
    typedef struct
    {
        int t;
        double sum1, sum2, sum3, sum4;
    } Recordgenoport;

    Recordgenoport *genorepo = malloc(sizeof(Recordgenoport) * 9 * (tend + 1));
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
        double x11, x21, x22, x31, x32, x33, initT2P1;
    } Localpair;
    Localpair *Pair = malloc(sizeof(Localpair) * 9 * (tend + 1));
    int pair_count = 0;
    Localpair *Pair_freq = malloc(sizeof(Localpair) * 9 * (tend + 1));
    int freq_count = 0;

    // 閾値記述
    typedef struct
    {
        int situ;
        double k,u,gen1,gen2,gen3,gen4;
    } Threshold;
    Threshold *Threshold1 = malloc(sizeof(Threshold) * 170);
    int Threshold_count = 0;
    //"%d\t%d\t%d\t%d\n",i,j,mgenotype,fgenotype
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
    data_file9 = "Twoalleles_threshold_result_Rogistic.csv";

    for (iK = 0; iK <= 7; iK++)
    {
        // K=(double)(iK*2-1)*0.00;
        if (iK == 0)
            K = 0.0;
        else if (iK == 1)
            K = 0.01;
        else if (iK == 2)
            K = 0.05;
        else if (iK == 3)
            K = 0.1;
        else if (iK == 4)
            K = 0.15;
        else if (iK == 5)
            K = 0.2;
        else if (iK == 6)
            K = 0.3;
        else if (iK == 7)
            K = 0.5;
        // K = 0.1 + (double)iK * 0.01; // K=0.05~0.20まで0.01刻み

        s1 = s2 = curs1 = curs2 = curs3 = 2.0;
        k_max = 0.9;
        situ = 0;
        printf("K:%f\n", K);
        count = 0;
        for (iu = 0; iu <= 20; iu++)
        {
            u=0.0+(double)iu*0.025;
            // printf("更新前 iu:%d,u:%f,curs1:%f,curs2:%f,curs3:%f,count:%dsitu:%d\n",
            //        iu, u, curs1, curs2, curs3, count, situ);

            // if (count == 0)
            // {
            //     if (iu == 0)
            //     { // 一番最初
            //         u = 0.0;
            //     }
            //     else if(iu>8){
            //         s1=curs1;
            //         count=1;
            //         situ=0;
            //         continue;
            //     }
            //     else if (iu == 1 && situ == 2)
            //     { // U=0で共存で終了
            //         s1 = 0.0;
            //         curs1 = 0.0;
            //         curs2 = 0.0;
            //         count = 1;
            //         situ = 0;
            //         continue;
            //     }
            //     else if (iu == 1 && situ == 3)
            //     { // U=0で絶滅(T1P1のみ)で終了
            //         s1 = 0.0;
            //         s2 = 0.0;
            //         curs1 = curs2 = 0.0;
            //         curs3 = 0.0;
            //         count = 1;
            //         situ = 0;
            //         continue;
            //     }

            //     else if (fabs(curs1 - k_max) < 0.00001 && fabs(curs2 - 2.0) < 0.00001 && fabs(curs3 - 2.0) < 0.00001)
            //     { // １回も共存していないため共存状態を探す
            //         s1 = k_max;
            //         s2 = k_max;
            //         curs1 = k_max;
            //         curs2 = k_max;
            //         curs3 = k_max;
            //         count = 1;
            //         situ = 0;
            //         continue;
            //     }
            //     else if (situ == 1 && fabs(curs2 - 2.0) < 0.00001 && fabs(curs3 - 2.0) < 0.00001)
            //     { // １回も共存していないため共存状態を探す

            //         curs1 = u;
            //         u = u + 0.5;
            //         if (u > k_max)
            //             u = k_max;
            //     }

            //     else if (fabs(curs2 - curs1) < 0.01)
            //     { // 十分に近づいた場合閾値発見で終了
            //         s1 = curs1;

            //         count = 1;
            //         situ = 0;
            //         continue;
            //     }
            //     else if (fabs(curs3 - curs1) < 0.01&& fabs(curs2 - 2.0) < 0.00001)
            //     { // 十分に近づいた場合閾値発見で終了
            //         s1 = curs1;
            //         s2=curs1;

            //         count = 1;
            //         situ = 0;
            //         continue;
            //     }

            //     else if (fabs(curs2 - 2.0) < 0.00001 && situ == 1)
            //     { // uではT2P1orT2P2のみである場合
            //         curs1 = u;
            //         u = (curs1 + curs3) / 2.0;
            //     }
            //     else if(situ==1){
            //         curs1 = u;
            //         u = (curs1 + curs2) / 2.0;
            //     }
            //     else if (situ == 2)
            //     { // uで共存の場合
            //         curs2 = u;
            //         u = (curs1 + curs2) / 2.0;
            //     }
            //     else if (situ == 3)
            //     { // uで絶滅である場合
            //         curs3 = u;
            //         u = (curs1 + u) / 2.0;
            //     }
            // }
            // else if (count == 1)
            // {

            //     if ((fabs(s1 - 0.0) < 0.00001 && fabs(s2 - 0.0) < 0.00001) || (fabs(s1 - 1.0) < 0.00001 && fabs(s2 - 1.0) < 0.00001))
            //     { // 絶滅閾値発見で終了
            //         count = 2;
            //         break;
            //     }
            //     else if(iu>16){
            //         s2=curs2;
            //         if(fabs(curs2-2.0)<0.00001){
            //             s2=curs3;
            //             break;
            //         }
            //         count=2;
            //         break;
            //     }
            //     else if(fabs(s1 - s2) < 0.00001 ){
            //         count = 2;
            //         break;
            //     }
            //     else if (fabs(curs2 - k_max) < 0.00001 && fabs(curs3 - 2.0) < 0.00001)
            //     { // 共存状態がなかったので終了
            //         s2 = k_max;
            //         curs3 = 1.0;
            //         count = 2;
            //         break;
            //     }

            //     else if (situ == 2 && fabs(curs3 - 2.0) < 0.00001)
            //     { // u=1.0で共存で終了
            //         curs2 = u;
            //         u = u + 0.5;
            //         if (u > k_max)
            //             u = k_max;
            //     }
            //     else if (situ == 0 && fabs(curs3 - 2.0) < 0.00001)
            //     { // u=1.0で共存で終了
            //         u = k_max;
            //     }
            //     else if(situ==0){
            //         if(curs3<curs2){
            //             s2=curs2;
            //             count=2;
            //             break;
            //         }else{
            //             u = (curs2 + curs3) / 2.0;
            //         }
            //     }

            //     else if (fabs(curs3 - curs2) < 0.01)
            //     { // 十分に近づいた場合閾値発見で終了
            //         s2 = curs3;

            //         count = 2;
            //         break;
            //     }

            //     else if (situ == 3)
            //     { // uで絶滅の場合
            //         curs3 = u;
            //         u = (curs2 + curs3) / 2.0;
            //     }
            //     else if (situ == 2)
            //     { // uで共存の場合
            //         curs2 = u;
            //         u = (curs2 + curs3) / 2.0;
            //     }
            //     else if (situ == 1)
            //     {
            //         s2 = curs2;
            //         count = 2;
            //         break;
            //     }
            // }
            // printf("更新後 iu:%d,u:%f,curs1:%f,curs2:%f,curs3:%f,count:%dsitu:%d\n",
            //        iu, u, curs1, curs2, curs3, count, situ);
            //  u=0.3+(double)iu*0.1;

            init_genrand(0);

            // snapshot_file1 = malloc(100);
            // sprintf(snapshot_file1, "Three_env_stmap_%f_initP_%g.dat", K, mapinitP);

            // snapshot_file2 = malloc(100);
            // // sprintf(snapshot2, "efcs_intmap_%f.dat", K);

            // snapshot_file3 = malloc(100);
            // sprintf(snapshot_file3, "Three_env_finmap_%f_initP_%g.dat", K, mapinitP);

            data_file1 = malloc(100);
            sprintf(data_file1, "Two_env_2dime_T2P2_K_%f_u_%f.dat", K, u);

            data_file2 = malloc(100);
            sprintf(data_file2, "Two_env_2dime_T2P2_flow_K_%f_u_%f.dat", K, u);

            data_file4 = malloc(100);
            sprintf(data_file4, "Two_env_2dime_pair_K_%f_u_%f.dat", K, u);

            data_file5 = malloc(100);
            sprintf(data_file5, "Two_env_2dime_pair_flow_K_%f_u_%f.dat", K, u);

            data_file6 = malloc(100);
            sprintf(data_file6, "Two_env_2dime_pair_final_K_%f_u_%f.dat", K, u);

            data_file7 = malloc(100);
            sprintf(data_file7, "Two_env_genoport_K_%f_u_%f.dat", K, u);

            data_file8 = malloc(100);
            sprintf(data_file8, "Two_env_2dime_pair_full_K_%f_u_%f.dat", K, u);

            // snapshot_file2 = malloc(100);

            buf_count = 0;
            geno_count = 0;
            pair_count = 0;
            freq_count = 0;
            // Threshold_count = 0;
            for (k = 1; k <= 1; k++)
            {
                // initT2 = 0.1 * initialT;
                for (k2 = 1; k2 <= 1; k2++) // k2=1; k2<=9; k2++
                {
                    // data1 = fopen(data_file1, "a");
                    // data7 = fopen(data_file7, "a");
                    // initT2 = 0.1 * k2;
                    // initP2 = 0.1 * k2;
                    // initT2P1 = 0.1 * 2.0 * k2;
                    initT2P1 = 0.25;
                    for (i = 0; i < LH; i++)
                    {
                        for (j = 0; j < LV; j++)
                        {
                            // T2P1で分ける。他はランダム
                            rnd = genrand_real2();
                            if (rnd < initT2P1)
                            {
                                maleT[i][j] = 2;
                                maleP[i][j] = 1;
                            }
                            else
                            {
                                rnd2 = genrand_real2();
                                if (rnd2 < 0.333)
                                {
                                    maleT[i][j] = 1;
                                    maleP[i][j] = 1;
                                }
                                else if (rnd2 < 0.666)
                                {
                                    maleT[i][j] = 1;
                                    maleP[i][j] = 2;
                                }
                                else
                                {
                                    maleT[i][j] = 2;
                                    maleP[i][j] = 2;
                                }
                            }
                            rnd = genrand_real2();
                            if (rnd < initT2P1)
                            {
                                femaleT[i][j] = 2;
                                femaleP[i][j] = 1;
                            }
                            else
                            {
                                rnd2 = genrand_real2();
                                if (rnd2 < 0.333)
                                {
                                    femaleT[i][j] = 1;
                                    femaleP[i][j] = 1;
                                }
                                else if (rnd2 < 0.666)
                                {
                                    femaleT[i][j] = 1;
                                    femaleP[i][j] = 2;
                                }
                                else
                                {
                                    femaleT[i][j] = 2;
                                    femaleP[i][j] = 2;
                                }
                            }

                            // 通常
                            // rnd = genrand_real2();
                            // if (rnd < initT2)
                            //     maleT[i][j] = 2;
                            // else
                            //     maleT[i][j] = 1;
                            // rnd = genrand_real2();
                            // if (rnd < initP2)
                            //     maleP[i][j] = 2;
                            // else
                            //     maleP[i][j] = 1;
                            // rnd = genrand_real2();
                            // if (rnd < initT2)
                            //     femaleT[i][j] = 2;
                            // else
                            //     femaleT[i][j] = 1;
                            // rnd = genrand_real2();
                            // if (rnd < initP2)
                            //     femaleP[i][j] = 2;
                            // else
                            //     femaleP[i][j] = 1;

                            // T1P1vsT2P1

                            // T2P1vsT2P2

                            // T1P1vsT2P2
                        }
                    }

                    // 最初の割合を出力
                    sum1 = sum2 = sum3 = sum4 = 0.0;
                    for (i = 0; i < 6; i++)
                        sum_pair[i] = 0.0;
                    // data7=fopen(data_file7,"a");
                    for (i = 0; i < LH; i++)
                    {
                        for (j = 0; j < LV; j++)
                        {
                            //
                            if (maleT[i][j] == 1 && maleP[i][j] == 1)
                                sum1 += 1.0;
                            else if (maleT[i][j] == 1 && maleP[i][j] == 2)
                                sum2 += 1.0;
                            else if (maleT[i][j] == 2 && maleP[i][j] == 1)
                                sum3 += 1.0;
                            else if (maleT[i][j] == 2 && maleP[i][j] == 2)
                                sum4 += 1.0;
                            calc_pair(sum_pair, i, j, maleT, maleP);
                        }
                    }
                    // fprintf(data7, "%d\t%lf\t%lf\t%lf\t%lf\n", 0,(double)sum1/(double)(LH*LV),(double)sum2/(double)(LH*LV),\
                    // (double)sum3/(double)(LH*LV),(double)sum4/(double)(LH*LV));
                    genorepo[geno_count].t = 0;
                    genorepo[geno_count].sum1 = (double)sum1 / (double)(LH * LV);
                    genorepo[geno_count].sum2 = (double)sum2 / (double)(LH * LV);
                    genorepo[geno_count].sum3 = (double)sum3 / (double)(LH * LV);
                    genorepo[geno_count].sum4 = (double)sum4 / (double)(LH * LV);
                    geno_count++;

                    buffer[buf_count].t = 0;
                    buffer[buf_count].geno1 = (double)sum1 / (double)(LH * LV);
                    buffer[buf_count].geno2 = (double)sum2 / (double)(LH * LV);
                    buffer[buf_count].geno3 = (double)sum3 / (double)(LH * LV);
                    buffer[buf_count].geno4 = (double)sum4 / (double)(LH * LV);
                    buffer[buf_count].initT2P1 = initT2P1;
                    buf_count++;

                    Pair[pair_count].t = 0;
                    Pair[pair_count].x11 = sum_pair[0] / ((double)4.0 * sum1);
                    Pair[pair_count].x21 = sum_pair[1] / ((double)4.0 * sum1);
                    Pair[pair_count].x31 = sum_pair[2] / ((double)4.0 * sum1);
                    Pair[pair_count].x22 = sum_pair[3] / ((double)4.0 * sum3);
                    Pair[pair_count].x32 = sum_pair[4] / ((double)4.0 * sum3);
                    Pair[pair_count].x33 = sum_pair[5] / ((double)4.0 * sum4);
                    Pair[pair_count].initT2P1 = initT2P1;
                    pair_count++;

                    Pair_freq[freq_count].t = 0;
                    Pair_freq[freq_count].x11 = sum_pair[0] / ((double)4.0 * sum1);
                    Pair_freq[freq_count].x21 = sum_pair[1] / ((double)4.0 * sum1);
                    Pair_freq[freq_count].x31 = sum_pair[2] / ((double)4.0 * sum1);
                    Pair_freq[freq_count].x22 = sum_pair[3] / ((double)4.0 * sum3);
                    Pair_freq[freq_count].x32 = sum_pair[4] / ((double)4.0 * sum3);
                    Pair_freq[freq_count].x33 = sum_pair[5] / ((double)4.0 * sum4);
                    Pair_freq[freq_count].initT2P1 = initT2P1;
                    freq_count++;
                    for (i = 0; i < 6; i++)
                        sum_pair[i] = 0.0;

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

                    for (t = 1; t <= tend; t++) // while(t<tend)
                    {

#pragma omp parallel default(none) shared(maleT, maleP, femaleT, femaleP, maleTdummy, malePdummy, femaleTdummy, femalePdummy, rng_states, K, u) \
    private(i, j, sum, gsum, rnd2, rnd4, rnd5, maleT0, maleP0, femaleT0, femaleP0)
                        {
                            int tid = omp_get_thread_num();
                            mt_state *rng = &rng_states[tid];
#pragma omp for collapse(2) schedule(static)
                            for (i = 0; i < LH; i++)
                                for (j = 0; j < LV; j++)
                                {

                                    calc_male_sum(sum, i, j, maleT, maleP, femaleP[i][j]);

                                    // メス遺伝のための重みづけ
                                    // メスのコスト

                                    calc_female_sum(gsum, i, j, femaleT, femaleP);
                                    // printf("calc_female_sum OK");

                                    rnd4 = genrand_real2_mt(rng);
                                    // オス遺伝
                                    if (rnd4 < 0.5)
                                    {
                                        rnd5 = genrand_real2_mt(rng);
                                        if (rnd5 < 0.5)
                                        { // 次世代がオスの場合
                                            do
                                            {
                                                rnd2 = genrand_real2_mt(rng);
                                                genotype(sum, &maleT0, &maleP0, rng);
                                                // printf("genotype OK");
                                                if (maleT0 == 1)
                                                    break;
                                                if (maleT0 == 2 && rnd2 < u)
                                                    continue;
                                                break;
                                            } while (1);

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
                                                if (maleP0 == 2 && rnd2 < K)
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
                                                if (femaleT0 == 2 && rnd2 < u)
                                                {
                                                    do
                                                    {
                                                        rnd2 = genrand_real2_mt(rng);
                                                        genotype(gsum, &femaleT0, &femaleP0, rng);
                                                        // printf("genotype OK");
                                                        if (femaleT0 == 1)
                                                            break;
                                                        if (femaleT0 == 2 && rnd2 < u)
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
                                                if (femaleP0 == 2 && rnd2 < K)
                                                {
                                                    do
                                                    {
                                                        rnd2 = genrand_real2_mt(rng);
                                                        genotype(gsum, &femaleT0, &femaleP0, rng);
                                                        // printf("genotype OK");
                                                        if (femaleP0 == 1)
                                                            break;
                                                        if (femaleP0 == 2 && rnd2 < K)
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
                        // if (t % 100 == 0 && fabs(initT2P1 - mapinitP) < 1e-9 && t < 40000)
                        // {
                        //     sprintf(snapshot_file2, "Two_intmap_t_%d_K_%f_u_%f_initT2P1_%g.dat", t, K, u, initT2P1);

                        //     mgenotype = fgenotype = 0;

                        //     for (i = 0; i < LH; i++)
                        //     {
                        //         for (j = 0; j < LV; j++)
                        //         {
                        //             if (maleT[i][j] == 1 && maleP[i][j] == 1)
                        //                 mgenotype = 1;
                        //             else if (maleT[i][j] == 1 && maleP[i][j] == 2)
                        //                 mgenotype = 2;
                        //             else if (maleT[i][j] == 2 && maleP[i][j] == 1)
                        //                 mgenotype = 3;
                        //             else if (maleT[i][j] == 2 && maleP[i][j] == 2)
                        //                 mgenotype = 4;
                        //             if (femaleT[i][j] == 1 && femaleP[i][j] == 1)
                        //                 fgenotype = 1;
                        //             else if (femaleT[i][j] == 1 && femaleP[i][j] == 2)
                        //                 fgenotype = 2;
                        //             else if (femaleT[i][j] == 2 && femaleP[i][j] == 1)
                        //                 fgenotype = 3;
                        //             else if (femaleT[i][j] == 2 && femaleP[i][j] == 2)
                        //                 fgenotype = 4;
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
                        //     Map("male", snapshot_file2, K, mapinitP, t);
                        //     Map("female", snapshot_file2, K, mapinitP, t);
                        // }

                        // 遺伝子型の割合出力
                        sum1 = sum2 = sum3 = sum4 = 0.0;
                        // data7=fopen(data_file7,"a");
                        for (i = 0; i < LH; i++)
                        {
                            for (j = 0; j < LV; j++)
                            {
                                if ((maleT[i][j] == 1 && maleP[i][j] == 1))
                                    sum1 += 1;
                                else if ((maleT[i][j] == 1 && maleP[i][j] == 2))
                                    sum2 += 1;
                                else if ((maleT[i][j] == 2 && maleP[i][j] == 1))
                                    sum3 += 1;
                                else if ((maleT[i][j] == 2 && maleP[i][j] == 2))
                                    sum4 += 1;
                                calc_pair(sum_pair, i, j, maleT, maleP);
                            }
                        }
                        genorepo[geno_count].t = t;
                        genorepo[geno_count].sum1 = (double)sum1 / (double)(LH * LV);
                        genorepo[geno_count].sum2 = (double)sum2 / (double)(LH * LV);
                        genorepo[geno_count].sum3 = (double)sum3 / (double)(LH * LV);
                        genorepo[geno_count].sum4 = (double)sum4 / (double)(LH * LV);
                        geno_count++;

                        if (t % 100 == 0)
                        {
                            Pair[pair_count].t = t;
                            if (sum1 == 0)
                            {
                                Pair[pair_count].x11 = 0.0;
                                Pair[pair_count].x21 = 0.0;
                                Pair[pair_count].x31 = 0.0;
                            }
                            if (sum3 == 0)
                            {
                                Pair[pair_count].x22 = 0.0;
                                Pair[pair_count].x32 = 0.0;
                            }
                            if (sum4 == 0)
                            {
                                Pair[pair_count].x33 = 0.0;
                            }

                            Pair[pair_count].x11 = sum_pair[0] / ((double)4.0 * sum1);
                            Pair[pair_count].x21 = sum_pair[1] / ((double)4.0 * sum1);
                            Pair[pair_count].x31 = sum_pair[2] / ((double)4.0 * sum1);
                            Pair[pair_count].x22 = sum_pair[3] / ((double)4.0 * sum3);
                            Pair[pair_count].x32 = sum_pair[4] / ((double)4.0 * sum3);
                            Pair[pair_count].x33 = sum_pair[5] / ((double)4.0 * sum4);
                            Pair[pair_count].initT2P1 = initT2P1;
                            pair_count++;

                            buffer[buf_count].t = t;
                            buffer[buf_count].geno1 = (double)sum1 / (double)(LH * LV);
                            buffer[buf_count].geno2 = (double)sum2 / (double)(LH * LV);
                            buffer[buf_count].geno3 = (double)sum3 / (double)(LH * LV);
                            buffer[buf_count].geno4 = (double)sum4 / (double)(LH * LV);
                            buffer[buf_count].initT2P1 = initT2P1;
                            buf_count++;
                        }
                        Pair_freq[freq_count].t = t;
                        if (sum1 == 0)
                        {
                            Pair_freq[freq_count].x11 = 0.0;
                            Pair_freq[freq_count].x21 = 0.0;
                            Pair_freq[freq_count].x31 = 0.0;
                        }
                        if (sum3 == 0)
                        {
                            Pair_freq[freq_count].x22 = 0.0;
                            Pair_freq[freq_count].x32 = 0.0;
                        }
                        if (sum4 == 0)
                        {
                            Pair_freq[freq_count].x33 = 0.0;
                        }

                        Pair_freq[freq_count].x11 = sum_pair[0] / ((double)4.0 * sum1);
                        Pair_freq[freq_count].x21 = sum_pair[1] / ((double)4.0 * sum1);
                        Pair_freq[freq_count].x31 = sum_pair[2] / ((double)4.0 * sum1);
                        Pair_freq[freq_count].x22 = sum_pair[3] / ((double)4.0 * sum3);
                        Pair_freq[freq_count].x32 = sum_pair[4] / ((double)4.0 * sum3);
                        Pair_freq[freq_count].x33 = sum_pair[5] / ((double)4.0 * sum4);
                        Pair_freq[freq_count].initT2P1 = initT2P1;
                        freq_count++;
                        for (i = 0; i < 6; i++)
                            sum_pair[i] = 0.0;
                    }
                }
            }

            data1 = fopen(data_file1, "w");
            for (int n = 0; n < buf_count; n++)
            {
                fprintf(data1, "%d\t%f\t%f\t%f\t%f\t%f\n",
                        buffer[n].t, buffer[n].geno1, buffer[n].geno2, buffer[n].geno3, buffer[n].geno4, buffer[n].initT2P1);
            }
            fclose(data1);

            data7 = fopen(data_file7, "w");
            for (int n = 0; n < geno_count; n++)
            {
                fprintf(data7, "%d\t%lf\t%lf\t%lf\t%lf\n",
                        genorepo[n].t, genorepo[n].sum1, genorepo[n].sum2, genorepo[n].sum3, genorepo[n].sum4);
            }
            fclose(data7);

            data4 = fopen(data_file4, "w");
            for (int n = 0; n < pair_count; n++)
            {
                fprintf(data4, "%d\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\n",
                        Pair[n].t, Pair[n].x11, Pair[n].x21, Pair[n].x31, Pair[n].x22,
                        Pair[n].x32, Pair[n].x33,
                        genorepo[n].sum1, genorepo[n].sum3, genorepo[n].sum4, Pair[n].initT2P1);
            }
            fclose(data4);

            data8 = fopen(data_file8, "w");
            for (int n = 0; n < freq_count; n++)
            {
                fprintf(data8, "%d\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\n",
                        Pair_freq[n].t, Pair_freq[n].x11, Pair_freq[n].x21,
                        Pair_freq[n].x31, Pair_freq[n].x22,
                        Pair_freq[n].x32, Pair_freq[n].x33,
                        Pair_freq[n].initT2P1);
            }
            fclose(data8);

            // T1P1-T2P2-T2P1図
            data_file3 = malloc(100);
            sprintf(data_file3, "Two_env_2dime_final_T2P2_K_%f_u_%f.dat", K, u);
            gp = fopen(data_file1, "r");
            data3 = fopen(data_file3, "w");
            while (fscanf(gp, "%d %lf %lf %lf %lf %lf", &x1, &gen1, &gen2, &gen3, &gen4, &init) == 6)
            {
                if (x1 == (tend))
                {

                    fprintf(data3, "%d\t%f\t%f\t%f\t%f\n", x1, gen1, gen2, gen3, gen4);
                }
            }
            fclose(data3);
            fclose(gp);

            data3 = fopen(data_file3, "r");
            if (!data3)
            {
                fprintf(stderr, "Error: Cannot open file\n");
                return 1;
            }
            if (fscanf(data3, "%d %lf %lf %lf %lf", &x1, &gen1, &gen2, &gen3, &gen4) != 5)
                return 1;
            // fscanf(data3, "%d %lf %lf %lf %lf", &x1, &gen1, &gen2, &gen3, &gen4);
            if ((fabs(gen1 - 0.0) < 1e-6 && fabs(gen2 - 0.0) < 1e-6&& fabs(gen3 - 0.0) < 1e-6) || (fabs(gen1 - 0.0) < 1e-6 && fabs(gen2 - 0.0) < 1e-6&& fabs(gen4 - 0.0) < 1e-6))
            {
                // 全部同じ遺伝子型になった場合、流れを出力しない
                situ = 1;
            }
            else if (fabs(gen1 - 1.0) < 1e-6)
            {
                // 全部同じ遺伝子型になった場合、流れを出力しない
                situ = 3;
            }
            else
            {
                situ = 2;
            }
            fclose(data3);

            data1 = fopen(data_file1, "r");
            data2 = fopen(data_file2, "w");
            if (fscanf(data1, "%d %lf %lf %lf %lf %lf", &x1, &gen1, &gen2, &gen3, &gen4, &init) != 6)
                return 1;
            while (fscanf(data1, "%d %lf %lf %lf %lf %lf", &x2, &geno1, &geno2, &geno3, &geno4, &init1) == 6)
            {
                if (fabs(init - init1) < 1e-9)
                {
                    fprintf(data2, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", gen1, gen3, gen4, geno1, geno3, geno4);
                }
                x1 = x2;
                gen1 = geno1;
                gen2 = geno2;
                gen3 = geno3;
                gen4 = geno4;
                init = init1;
            }
            fclose(data1);
            fclose(data2);

            // // T1P1-T2P1図
            // gp = popen("gnuplot -persist", "w");
            // fprintf(gp, "set terminal png\n");
            // fprintf(gp, "set term pngcairo size 1000,700\n");
            // fprintf(gp, "set output 'Genotype_Twoalleles/K_%f_a1_%f_T1P1_T2P1.png'\n", K, a1);
            // fprintf(gp, "set xrange [0:%f]\n", 1.0);
            // fprintf(gp, "set xlabel 'T1P1'\n");
            // fprintf(gp, "set yrange [0:%f]\n", 1.0);
            // fprintf(gp, "set ylabel 'T2P1'\n");
            // fprintf(gp, "plot \'%s\' using 2:4 with points pointtype 7 lc rgb 'blue' title \"survivalrateK=%f\",\'%s\' using 1:2:($4-$1):($5-$2) with vectors head filled lc rgb 'blue',\'%s\' using 2:4 with points pointtype 7 lc rgb 'red' title \"finalarrival\"\n", data_file1, K, data_file2, data_file3);
            // pclose(gp);

            // // T1P1-T2P2図
            // gp = popen("gnuplot -persist", "w");
            // fprintf(gp, "set terminal png\n");
            // fprintf(gp, "set term pngcairo size 1000,700\n");
            // fprintf(gp, "set output 'Genotype_Twoalleles/K_%f_a1_%f_T1P1_T2P2.png'\n", K, a1);
            // fprintf(gp, "set xrange [0:%f]\n", 1.0);
            // fprintf(gp, "set xlabel 'T1P1'\n");
            // fprintf(gp, "set yrange [0:%f]\n", 1.0);
            // fprintf(gp, "set ylabel 'T2P2'\n");
            // fprintf(gp, "plot \'%s\' using 2:5 with points pointtype 7 lc rgb 'blue' title \"survivalrateK=%f\",\'%s\' using 1:3:($4-$1):($6-$3) with vectors head filled lc rgb 'blue',\'%s\' using 2:5 with points pointtype 7 lc rgb 'red' title \"finalarrival\"\n", data_file1, K, data_file2, data_file3);
            // pclose(gp);

            // for (i = 1; i <= 9; i++)
            // {
            //     initT2P1 = (double)0.1 * i;
            //     gp = popen("gnuplot -persist", "w");
            //     fprintf(gp, "set terminal png\n");
            //     fprintf(gp, "set output 'Genotype_twoalleles_genoport/Two_env_genoport_K_%f_a1_%f_initT2P1_%f.png'\n", K, a1, initT2P1);
            //     fprintf(gp, "set xrange [0:%d]\n", tend);
            //     fprintf(gp, "set xlabel 't'\n");
            //     fprintf(gp, "set yrange [0:%f]\n", 1.0);
            //     fprintf(gp, "set ylabel 'genotype_frequency'\n");
            //     fprintf(gp, "titles='T1P1 T1P2 T2P1 T2P2'\n");
            //     fprintf(gp, "set style line 1 lc rgb \"#0000FF\" lw 2\n");
            //     fprintf(gp, "set style line 2 lc rgb \"#00CC00\" lw 2\n");
            //     fprintf(gp, "set style line 3 lc rgb \"#FF8800\" lw 2\n");
            //     fprintf(gp, "set style line 4 lc rgb \"#FF0000\" lw 2\n");
            //     fprintf(gp, "plot for [j=2:5] \'%s\' every ::%d::%d using 1:j with lines ls (j-1) title word(titles, j-1)\n", data_file7, (i - 1) * (tend + 1), i * (tend + 1) - 1);
            //     pclose(gp);
            // }

            // pairのグラフ
            gp = fopen(data_file4, "r");
            data6 = fopen(data_file6, "w");
            while (fscanf(gp, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &x1, &x11, &x21, &x31, &x22, &x32, &x33, &sum1, &sum3, &sum4, &init) == 11)
            {
                if (x1 == (tend))
                {

                    fprintf(data6, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", x1, x11, x21, x31, x22, x32, x33, sum1, sum3, sum4);
                }
            }
            fclose(data6);
            fclose(gp);

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

            // for (i = 0; i < 9; i++)
            // {
            //     for (j = i + 1; j < 9; j++)
            //     {
            //         gp = popen("gnuplot -persist", "w");
            //         fprintf(gp, "set terminal png\n");
            //         fprintf(gp, "set term pngcairo size 1000,700\n");
            //         fprintf(gp, "set output 'Genotype_Twoalleles_pair/K_%f_a1_%f_%s_%s.png'\n", K, a1, Figaxis[i], Figaxis[j]);
            //         fprintf(gp, "set xrange [0:%f]\n", 1.0);
            //         fprintf(gp, "set xlabel \'%s\'\n", Figaxis[i]);
            //         fprintf(gp, "set yrange [0:%f]\n", 1.0);
            //         fprintf(gp, "set ylabel \'%s\'\n", Figaxis[j]);
            //         fprintf(gp, "plot \'%s\' using %d:%d with points pointtype 7 lc rgb 'blue' title \
            //             \"survivalrateK=%f\",\
            //             \'%s\' using %d:%d:($%d-$%d):($%d-$%d) with vectors head filled lc rgb 'blue',\
            //             \'%s\' using %d:%d with points pointtype 7 lc rgb 'red' title \"finalarrival\"\n",
            //                 data_file4, i + 2, j + 2,
            //                 K,
            //                 data_file5, i + 1, j + 1, i + 10, i + 1, j + 10, j + 1,
            //                 data_file6, i + 2, j + 2);
            //         pclose(gp);
            //     }
            // }

            for (i = 1; i <= 1; i++) // i = 1; i <= 4; i++
            {
                // initT2P1 = (double)0.2 * i;
                initT2P1 = (double)0.25 * i;
                gp = popen("gnuplot -persist", "w");
                fprintf(gp, "set terminal png\n");
                fprintf(gp, "set output 'Genotype_Twoalleles_pair/Two_env_pair_K_%f_u_%f_a1_%f_initT2P1_%f.png'\n", K, u, a1, initT2P1);
                fprintf(gp, "set xrange [0:%d]\n", tend);
                fprintf(gp, "set xlabel 't'\n");
                fprintf(gp, "set yrange [0:%f]\n", 1.0);
                fprintf(gp, "set ylabel 'pair_frequency'\n");
                fprintf(gp, "titles='x_1/1 x_2/1 x_3/1 x_2/2 x_3/2 x_3/3'\n");
                fprintf(gp, "set style line 1 lc rgb \"#0000FF\" lw 2\n");
                fprintf(gp, "set style line 2 lc rgb \"#00CC00\" lw 2\n");
                fprintf(gp, "set style line 3 lc rgb \"#FF8800\" lw 2\n");
                fprintf(gp, "set style line 4 lc rgb \"#FF0000\" lw 2\n");
                fprintf(gp, "set style line 5 lc rgb \"#FF00FF\" lw 2\n");
                fprintf(gp, "set style line 6 lc rgb \"#000000\" lw 2\n");
                fprintf(gp, "plot for [j=2:7] \'%s\' every ::%d::%d using 1:j with lines ls (j-1) title word(titles, j-1)\n", data_file8, (i - 1) * (tend + 1), i * (tend + 1) - 1);
                pclose(gp);
            }

            // free(snapshot_file1);
            // free(snapshot_file2);
            // free(snapshot_file3);
            free(data_file1);
            free(data_file2);
            free(data_file3);
            free(data_file4);
            free(data_file5);
            free(data_file6);
            free(data_file7);
            free(data_file8);

            
            Threshold1[Threshold_count].k = K;
            Threshold1[Threshold_count].u = u;
            Threshold1[Threshold_count].situ = situ;
            Threshold1[Threshold_count].gen1 = genorepo[geno_count-1].sum1;
            Threshold1[Threshold_count].gen2 = genorepo[geno_count-1].sum2;
            Threshold1[Threshold_count].gen3 = genorepo[geno_count-1].sum3;
            Threshold1[Threshold_count].gen4 = genorepo[geno_count-1].sum4;
            Threshold_count = Threshold_count + 1;
            printf("K=%f u=%f situ=%d T1P1=%f T1P2=%f T2P1=%f T2P2=%f\n", K, u, situ,\
                 genorepo[geno_count-1].sum1, genorepo[geno_count-1].sum2,\
                  genorepo[geno_count-1].sum3, genorepo[geno_count-1].sum4);
        }
        
    }
    gp = fopen(data_file9, "w");
    fprintf(gp, "K,u,situ,freq_T1P1,freq_T1P2,freq_T2P1,freq_T2P2\n");
    for (i = 0; i < Threshold_count; i++)
    {
        fprintf(gp, "%f %f %d %f %f %f %f\n", Threshold1[i].k, Threshold1[i].u, Threshold1[i].situ, Threshold1[i].gen1, Threshold1[i].gen2, Threshold1[i].gen3, Threshold1[i].gen4);
    }
    fclose(gp);

    // gp = popen("gnuplot -persist", "w");
    // fprintf(gp, "set terminal png\n");
    // fprintf(gp, "set output 'Genotype_Twoalleles/Threshold_initT2P1_%f.png'\n", initT2P1);
    // fprintf(gp, "set xrange [0:%f]\n", 0.5);
    // fprintf(gp, "set xlabel 'Female cost(K)'\n");
    // fprintf(gp, "set yrange [0:%f]\n", 0.9);
    // fprintf(gp, "set ylabel 'Male cost(U)'\n");
    // fprintf(gp, "plot \'%s\' using 1:2 with lines lw 2 lc rgb 'blue' title 'Threshold_coexistance',\'%s\' using 1:3 with lines lw 2 lc rgb 'red' title 'Threshold_extinction'\n", data_file9, data_file9);
    // pclose(gp);

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
    free(buffer);
    free(genorepo);
    free(rng_states);
    free(Pair);
    free(Pair_freq);
    free(Threshold1);

    clock_gettime(CLOCK_MONOTONIC, &end);
    double elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) * 1e-9;
    printf("simulation elapsed: %.3f s\n", elapsed);
}
