/* nearest neighbor interaction */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "MT.h"

#define LH 1000 // 10000
#define LV 1000 // 10000
// #define K 0.075 //P3メスのコスト
// #define V 0.074 //P2メスのコスト(0<V<K)
#define u 0.3 // T3オスのコスト0.3
// #define l 0.15  //T2オスのコスト(0<l<u)
#define a1 3.0 // P2メスがT2オスを選好する倍率3.0
// #define a2 6.0    // P3メスがT3オスを選好する倍率
#define tend 100000 // 4000 80000 10000
#define mapinitP 0.6
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
            fprintf(gp, "set output 'Twoalleles_hyper_map/MapSt_%s_two_%g_initP_%g.png'\n", sex, K, initP);
        else if (strstr(filename, "intmap"))
            fprintf(gp, "set output 'Twoalleles_hyper_map/MapInt_%s_two_%g_initP_%g_t_%06d.png'\n", sex, K, initP, t);
        else if (strstr(filename, "finmap"))
            fprintf(gp, "set output 'Twoalleles_hyper_map/MapFin_%s_two_%g_initP_%g.png'\n", sex, K, initP);
        fprintf(gp, "unset key\n");
        fprintf(gp, "set size ratio -1\n");
        fprintf(gp, "set xrange [0:%d]\n", LH - 1);
        // fprintf(gp,"set xlabel 'T2'\n");
        fprintf(gp, "set yrange [0:%d]\n", LV - 1);
        fprintf(gp, "set palette defined(1 'blue',2 'green',3 'orange',4 'red')\n");
        // fprintf(gp,"set multiplot layout 1,2 title 'Genotype map (T×P: 0=T1P1, 1=T1P2, 2=T2P1, 3=T2P2)'\n");
        fprintf(gp, "set title 'Male map'\n");
        fprintf(gp, "unset xtics;unset ytics\n");
        // fprintf(gp,"unset yticks\n");
        fprintf(gp, "plot \'%s\' using 1:2:3 with image\n", filename);
    }
    else if (strcmp(sex, "female") == 0)
    {
        if (strstr(filename, "stmap"))
            fprintf(gp, "set output 'Twoalleles_hyper_map/MapSt_%s_two_%g_initP_%g.png'\n", sex, K, initP);
        else if (strstr(filename, "intmap"))
            fprintf(gp, "set output 'Twoalleles_hyper_map/MapInt_%s_two_%g_initP_%g_t_%06d.png'\n", sex, K, initP, t);
        else if (strstr(filename, "finmap"))
            fprintf(gp, "set output 'Twoalleles_hyper_map/MapFin_%s_two_%g_initP_%g.png'\n", sex, K, initP);
        fprintf(gp, "unset key\n");
        fprintf(gp, "set size ratio -1\n");
        fprintf(gp, "set xrange [0:%d]\n", LH - 1);
        // fprintf(gp,"set xlabel 'T2'\n");
        fprintf(gp, "set yrange [0:%d]\n", LV - 1);
        fprintf(gp, "set palette defined(1 'blue',2 'green',3 'orange',4 'red')\n");
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
    double initT2, initP2, initT3, initP3;
    int k, k2, i, j, i2, j2, t, ok, x1, x2, a, b, n;
    int maleI, maleJ, femaleI, femaleJ;
    int numMT1, numMT2, numMT3, numMP1, numMP2, numMP3;
    int numFT1, numFT2, numFT3, numFP1, numFP2, numFP3;
    double sum[4], gsum[4];
    double rnd, rnd2, rnd3, rnd4, rnd5, sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9, gen1, gen2, gen3, gen4, geno1, geno2, geno3, geno4, init, init1;
    double gsum1, gsum2, gsum3, gsum4, gsum5, gsum6, gsum7, gsum8, gsum9;
    int iK, iV, il, ia2, ia1;
    double K, V, l, a2;
    int maleT0, maleP0, femaleT0, femaleP0, mgenotype, fgenotype, count;
    FILE *gp, *data1, *data2, *data3, *data4, *data5, *data6, *data7;
    FILE *snapshot1, *snapshot2, *snapshot3, *snapshot4, *snapshot5, *snapshot6;
    char *data_file1, *data_file2, *data_file3, *data_file4, *data_file5, *data_file6, *data_file7;
    char *snapshot_file1, *snapshot_file2, *snapshot_file3, *snapshot_file4, *snapshot_file5, *snapshot_file6;
    // data1(T2P2頻度図の書き込み用)
    typedef struct
    {
        int t;
        double geno1, geno2, geno3, geno4, initP2;
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

    for (iK = 0; iK <= 2; iK++)
    {
        // K=(double)(iK*2-1)*0.00;
        // K = 0.1 + (double)iK * 0.01; // K=0.05~0.20まで0.01刻み
        K = 0.15 + (double)iK * 0.01;
        printf("K:%f\n", K);
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

        init_genrand(0);

        // snapshot_file1 = malloc(100);
        // sprintf(snapshot_file1, "Three_env_stmap_%f_initP_%g.dat", K, mapinitP);

        // snapshot_file2 = malloc(100);
        // // sprintf(snapshot2, "efcs_intmap_%f.dat", K);

        // snapshot_file3 = malloc(100);
        // sprintf(snapshot_file3, "Three_env_finmap_%f_initP_%g.dat", K, mapinitP);

        data_file1 = malloc(100);
        sprintf(data_file1, "Two_env_3dime_T2P2_K_%f.dat", K);

        data_file2 = malloc(100);
        sprintf(data_file2, "Two_env_3dime_T2P2_flow_K_%f.dat", K);

        data_file7 = malloc(100);
        sprintf(data_file7, "Two_env_genoport_K_%f.dat", K);

        snapshot_file2 = malloc(100);

        buf_count = 0;
        geno_count = 0;
        for (k = 1; k <= 1; k++)
        {
            initT2 = 0.1 * initialT;
            for (k2 = 1; k2 <= 9; k2++) // k2=1; k2<=9; k2++
            {
                // data1 = fopen(data_file1, "a");
                // data7 = fopen(data_file7, "a");

                initP2 = 0.1 * k2;
                for (i = 0; i < LH; i++)
                {
                    for (j = 0; j < LV; j++)
                    {
                        // 通常
                        rnd = genrand_real2();
                        if (rnd < initT2)
                            maleT[i][j] = 2;
                        else
                            maleT[i][j] = 1;
                        rnd = genrand_real2();
                        if (rnd < initP2)
                            maleP[i][j] = 2;
                        else
                            maleP[i][j] = 1;
                        rnd = genrand_real2();
                        if (rnd < initT2)
                            femaleT[i][j] = 2;
                        else
                            femaleT[i][j] = 1;
                        rnd = genrand_real2();
                        if (rnd < initP2)
                            femaleP[i][j] = 2;
                        else
                            femaleP[i][j] = 1;

                        // T1P1vsT2P1

                        // T2P1vsT2P2

                        // T1P1vsT2P2
                    }
                }

                numMT1 = numMT2 = numFT1 = numFT2 = numMP1 = numMP2 = numFP1 = numFP2 = 0;
                for (i = 0; i < LH; i++)
                {
                    for (j = 0; j < LV; j++)
                    {
                        if (maleT[i][j] == 1)
                            numMT1++;
                        else if (maleT[i][j] == 2)
                            numMT2++;
                        if (femaleT[i][j] == 1)
                            numFT1++;
                        else if (femaleT[i][j] == 2)
                            numFT2++;
                        if (maleP[i][j] == 1)
                            numMP1++;
                        else if (maleP[i][j] == 2)
                            numMP2++;
                        if (femaleP[i][j] == 1)
                            numFP1++;
                        else if (femaleP[i][j] == 2)
                            numFP2++;
                    }
                }
                // printf("0 T2 P2:%f %f \n", (double)numMT2 / (double)(LH * LV), (double)numMP2 / (double)(LH * LV));

                // 最初の割合を出力
                sum1 = sum2 = sum3 = sum4 = 0.0;
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
                // fclose(data7);

                // 初期の図
                // if (initP2 == mapinitP)
                // {
                //     snapshot1 = fopen(snapshot_file1, "w");
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
                //             fprintf(snapshot1, "%d\t%d\t%d\t%d\n", i, j, mgenotype, fgenotype);
                //         }
                //     }
                //     fclose(snapshot1);
                // }

                // t = 0;

                // --- 追加: dummy 配列を初期化（未書き込み領域を防ぐ） ---

                for (t = 0; t < tend; t++) // while(t<tend)
                {

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

                    // data1 = fopen(data_file1, "a");

                    if (t % 10 == 0)
                    {
                        buffer[buf_count].t = t;
                        buffer[buf_count].geno1 = (double)sum1 / (double)(LH * LV);
                        buffer[buf_count].geno2 = (double)sum2 / (double)(LH * LV);
                        buffer[buf_count].geno3 = (double)sum3 / (double)(LH * LV);
                        buffer[buf_count].geno4 = (double)sum4 / (double)(LH * LV);
                        buffer[buf_count].initP2 = initP2;
                        buf_count++;
                    }
                    // fprintf(data1, "%d\t%f\t%f\t%f\n", t, (double)numMT2 / (double)(LH * LV), (double)numMP2 / (double)(LH * LV),initP2);
                    // fclose(data1);

                    // #pragma omp parallel
                    // {
                    //     #pragma omp single
                    //     printf("OpenMP threads = %d\n", omp_get_num_threads());
                    // }
                    // --- 追加: dummy 配列を初期化（未書き込み領域を防ぐ） ---
                     // #pragma omp parallel for collapse(2) schedule(static) default(none) shared(maleT, maleP, femaleT, femaleP,                         \
                                                                               maleTdummy, malePdummy, femaleTdummy, femalePdummy) \
    private(i, j)

#pragma omp parallel for collapse(2) schedule(static) private(sum, gsum, rnd2, rnd4, rnd5, maleT0, maleP0, femaleT0, femaleP0)
                    // #pragma omp parallel for schedule(static)
                    for (i = 0; i < LH; i++)
                        for (j = 0; j < LV; j++)
                        {
                            int tid = omp_get_thread_num();
                            // if (i == 0 && j < 10) {
                            //     printf("Thread %d is working on (0,%d)\n", omp_get_thread_num(), j);
                            // }

                            calc_male_sum(sum, i, j, maleT, maleP, femaleP[i][j]);

                            // calc_male_sum(sum,i,j,maleT,maleP,femaleP[i][j]);
                            // for(n=0;n<9;n++){printf("sum[%d]:%f",n,sum[n]);}
                            // printf("calc_male_sum OK");
                            // メス遺伝のための重みづけ
                            // メスのコスト

                            calc_female_sum(gsum, i, j, femaleT, femaleP);
                            // printf("calc_female_sum OK");

                            rnd4 = genrand_real2_mt(&rng_states[tid]);
                            // オス遺伝
                            if (rnd4 < 0.5)
                            {
                                rnd5 = genrand_real2_mt(&rng_states[tid]);
                                if (rnd5 < 0.5)
                                { // 次世代がオスの場合
                                    do
                                    {
                                        rnd2 = genrand_real2_mt(&rng_states[tid]);
                                        genotype(sum, &maleT0, &maleP0, &rng_states[tid]);
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
                                        rnd2 = genrand_real2_mt(&rng_states[tid]);
                                        genotype(sum, &maleT0, &maleP0, &rng_states[tid]);
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
                                rnd5 = genrand_real2_mt(&rng_states[tid]);
                                if (rnd5 < 0.5)
                                { // 次世代がオス
                                    femaleT0 = femaleT[i][j];
                                    femaleP0 = femaleP[i][j];
                                    if (femaleT0 != 1)
                                    {
                                        rnd2 = genrand_real2_mt(&rng_states[tid]);
                                        if (femaleT0 == 2 && rnd2 < u)
                                        {
                                            do
                                            {
                                                rnd2 = genrand_real2_mt(&rng_states[tid]);
                                                genotype(gsum, &femaleT0, &femaleP0, &rng_states[tid]);
                                                // printf("genotype OK");
                                                if (femaleT0 == 1)
                                                    break;
                                                if (femaleT0 == 2 && rnd2 < u)
                                                    continue;
                                                break;
                                            } while (1);
                                        }
                                    }

                                    // do
                                    //     {
                                    //         rnd2 = genrand_real2_mt(&rng_states[tid]);
                                    //         genotype(gsum, &femaleT0, &femaleP0, &rng_states[tid]);
                                    //         // printf("genotype OK");
                                    //         if (femaleT0 == 1)
                                    //             break;
                                    //         if (femaleT0 == 2 && rnd2 < u)
                                    //             continue;
                                    //         break;
                                    //     } while (1);
                                    maleTdummy[i][j] = femaleT0;
                                    malePdummy[i][j] = femaleP0;
                                }
                                else
                                { // 次世代がメス
                                    femaleT0 = femaleT[i][j];
                                    femaleP0 = femaleP[i][j];
                                    if (femaleP0 != 1)
                                    {
                                        rnd2 = genrand_real2_mt(&rng_states[tid]);
                                        if (femaleP0 == 2 && rnd2 < K)
                                        {
                                            do
                                            {
                                                rnd2 = genrand_real2_mt(&rng_states[tid]);
                                                genotype(gsum, &femaleT0, &femaleP0, &rng_states[tid]);
                                                // printf("genotype OK");
                                                if (femaleP0 == 1)
                                                    break;
                                                if (femaleP0 == 2 && rnd2 < K)
                                                    continue;
                                                break;
                                            } while (1);
                                        }
                                    }
                                    // do
                                    //     {
                                    //         rnd2 = genrand_real2_mt(&rng_states[tid]);
                                    //         genotype(gsum, &femaleT0, &femaleP0,&rng_states[tid]);
                                    //         // printf("genotype OK");
                                    //         if (femaleP0 == 1)
                                    //             break;
                                    //         if (femaleP0 == 2 && rnd2 < K)
                                    //             continue;
                                    //         break;

                                    //     } while (1);
                                    femaleTdummy[i][j] = femaleT0;
                                    femalePdummy[i][j] = femaleP0;
                                }
                            }
                        }

                    // --- 追加: dummy 配列を初期化（未書き込み領域を防ぐ） ---
#pragma omp parallel for collapse(2) schedule(static) default(none) shared(maleT, maleP, femaleT, femaleP,                         \
                                                                               maleTdummy, malePdummy, femaleTdummy, femalePdummy) \
    private(i, j)
                    for (i = 0; i < LH; i++)
                    {
                        for (j = 0; j < LV; j++)
                        {
                            maleT[i][j] = maleTdummy[i][j];
                            maleP[i][j] = malePdummy[i][j];
                            femaleT[i][j] = femaleTdummy[i][j];
                            femaleP[i][j] = femalePdummy[i][j];
                        }
                    }
                    numMT1 = numMT2 = numFT1 = numFT2 = numMP1 = numMP2 = numFP1 = numFP2 = 0;
                    for (i = 0; i < LH; i++)
                    {
                        for (j = 0; j < LV; j++)
                        {
                            if (maleT[i][j] == 1)
                                numMT1++;
                            else if (maleT[i][j] == 2)
                                numMT2++;
                            if (femaleT[i][j] == 1)
                                numFT1++;
                            else if (femaleT[i][j] == 2)
                                numFT2++;
                            if (maleP[i][j] == 1)
                                numMP1++;
                            else if (maleP[i][j] == 2)
                                numMP2++;
                            if (femaleP[i][j] == 1)
                                numFP1++;
                            else if (femaleP[i][j] == 2)
                                numFP2++;
                        }
                    }

                    // 途中の図
                    // 途中の図
                    if (t % 100 == 0 && fabs(initP2 - mapinitP) < 1e-12 && t < 40000)
                    {
                        sprintf(snapshot_file2, "Two_intmap_t_%d_cost_%f_initP_%g.dat", t, K, mapinitP);

                        mgenotype = fgenotype = 0;

                        for (i = 0; i < LH; i++)
                        {
                            for (j = 0; j < LV; j++)
                            {
                                if (maleT[i][j] == 1 && maleP[i][j] == 1)
                                    mgenotype = 1;
                                else if (maleT[i][j] == 1 && maleP[i][j] == 2)
                                    mgenotype = 2;
                                else if (maleT[i][j] == 2 && maleP[i][j] == 1)
                                    mgenotype = 3;
                                else if (maleT[i][j] == 2 && maleP[i][j] == 2)
                                    mgenotype = 4;
                                if (femaleT[i][j] == 1 && femaleP[i][j] == 1)
                                    fgenotype = 1;
                                else if (femaleT[i][j] == 1 && femaleP[i][j] == 2)
                                    fgenotype = 2;
                                else if (femaleT[i][j] == 2 && femaleP[i][j] == 1)
                                    fgenotype = 3;
                                else if (femaleT[i][j] == 2 && femaleP[i][j] == 2)
                                    fgenotype = 4;
                                recomap[i * LH + j].i = i;
                                recomap[i * LH + j].j = j;
                                recomap[i * LH + j].mgenotype = mgenotype;
                                recomap[i * LH + j].fgenotype = fgenotype;
                            }
                        }
                        snapshot2 = fopen(snapshot_file2, "w");
                        for (n = 0; n < LH * LV; n++)
                        {
                            fprintf(snapshot2, "%d\t%d\t%d\t%d\n", recomap[n].i, recomap[n].j, recomap[n].mgenotype, recomap[n].fgenotype);
                        }
                        fclose(snapshot2);
                        Map("male", snapshot_file2, K, mapinitP, t);
                        Map("female", snapshot_file2, K, mapinitP, t);
                    }

                    // 遺伝子型の割合出力
                    sum1 = sum2 = sum3 = sum4 = 0.0;
                    // data7=fopen(data_file7,"a");
                    for (i = 0; i < LH; i++)
                    {
                        for (j = 0; j < LV; j++)
                        {
                            // if ((maleT[i][j] == 1 && maleP[i][j] == 1) || (femaleT[i][j] == 1 && femaleP[i][j] == 1))
                            //     sum1 += 1;
                            // else if ((maleT[i][j] == 1 && maleP[i][j] == 2) || (femaleT[i][j] == 1 && femaleP[i][j] == 2))
                            //     sum2 += 1;
                            // else if ((maleT[i][j] == 2 && maleP[i][j] == 1) || (femaleT[i][j] == 2 && femaleP[i][j] == 1))
                            //     sum3 += 1;
                            // else if ((maleT[i][j] == 2 && maleP[i][j] == 2) || (femaleT[i][j] == 2 && femaleP[i][j] == 2))
                            //     sum4 += 1;
                            if ((maleT[i][j] == 1 && maleP[i][j] == 1))
                                sum1 += 1;
                            else if ((maleT[i][j] == 1 && maleP[i][j] == 2))
                                sum2 += 1;
                            else if ((maleT[i][j] == 2 && maleP[i][j] == 1))
                                sum3 += 1;
                            else if ((maleT[i][j] == 2 && maleP[i][j] == 2))
                                sum4 += 1;
                        }
                    }
                    genorepo[geno_count].t = t;
                    genorepo[geno_count].sum1 = (double)sum1 / (double)(LH * LV);
                    genorepo[geno_count].sum2 = (double)sum2 / (double)(LH * LV);
                    genorepo[geno_count].sum3 = (double)sum3 / (double)(LH * LV);
                    genorepo[geno_count].sum4 = (double)sum4 / (double)(LH * LV);
                    geno_count++;
                    // fclose(data7);

                    // int **tmp;

                    // tmp = maleT;
                    // maleT = maleTdummy;
                    // maleTdummy = tmp;

                    // tmp = maleP;
                    // maleP = malePdummy;
                    // malePdummy = tmp;

                    // tmp = femaleT;
                    // femaleT = femaleTdummy;
                    // femaleTdummy = tmp;

                    // tmp = femaleP;
                    // femaleP = femalePdummy;
                    // femalePdummy = tmp;
                }

                // 最後の図
                // if (initP2 == mapinitP)
                // {
                //     snapshot3 = fopen(snapshot_file3, "w");
                //     mgenotype = fgenotype = 0;
                //     for (i = 0; i < LH; i++)
                //     {
                //         for (j = 0; j < LV; j++)
                //         {
                //             if (maleT[i][j] == 1 && maleP[i][j] == 1)
                //                 mgenotype = 1;
                //             else if (maleT[i][j] == 1 && maleP[i][j] == 2)
                //                 mgenotype = 2;
                //             else if (maleT[i][j] == 1 && maleP[i][j] == 2)
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
                //             fprintf(snapshot3, "%d\t%d\t%d\t%d\n", i, j, mgenotype, fgenotype);
                //         }
                //     }
                //     fclose(snapshot3);
                // }
                // fclose(data1);
                // fclose(data7);

            } // omp parallel
        }

        data1 = fopen(data_file1, "w");
        for (int n = 0; n < buf_count; n++)
        {
            fprintf(data1, "%d\t%f\t%f\t%f\t%f\t%f\n",
                    buffer[n].t, buffer[n].geno1, buffer[n].geno2, buffer[n].geno3, buffer[n].geno4, buffer[n].initP2);
        }
        fclose(data1);

        data7 = fopen(data_file7, "w");
        for (int n = 0; n < geno_count; n++)
        {
            fprintf(data7, "%d\t%lf\t%lf\t%lf\t%lf\n",
                    genorepo[n].t, genorepo[n].sum1, genorepo[n].sum2, genorepo[n].sum3, genorepo[n].sum4);
        }
        fclose(data7);

        printf("Ok\n");
        // T1P1-T2P2-T2P1図
        data_file3 = malloc(100);
        sprintf(data_file3, "Two_env_3dime_final_T2P2_K_%f.dat", K);
        // data_file3 = "final_env.dat";
        gp = fopen(data_file1, "r");
        data3 = fopen(data_file3, "w");
        while (fscanf(gp, "%d %lf %lf %lf %lf %lf", &x1, &gen1, &gen2, &gen3, &gen4, &init) == 6)
        {
            if (x1 == (tend - 10))
            {

                fprintf(data3, "%d\t%f\t%f\t%f\t%f\n", x1, gen1, gen2, gen3, gen4);
            }
        }
        fclose(data3);
        fclose(gp);

        data1 = fopen(data_file1, "r");
        data2 = fopen(data_file2, "w");
        if (fscanf(gp, "%d %lf %lf %lf %lf %lf", &x1, &gen1, &gen2, &gen3, &gen4, &init) != 6)
            return 1;
        while (fscanf(gp, "%d %lf %lf %lf %lf %lf", &x2, &geno1, &geno2, &geno3, &geno4, &init1) == 6)
        {
            if (fabs(init - init1) < 1e-12)
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

        gp = popen("gnuplot -persist", "w");
        fprintf(gp, "set terminal png\n");
        fprintf(gp, "set term pngcairo size 1000,700\n");
        fprintf(gp, "set output 'Genotype_Twoalleles_3dime/K_%f_a1_%f.png'\n", K, a1);
        fprintf(gp, "set xrange [0:%f]\n", 1.0);
        fprintf(gp, "set xlabel 'T1P1'\n");
        fprintf(gp, "set yrange [0:%f]\n", 1.0);
        fprintf(gp, "set ylabel 'T2P1'\n");
        fprintf(gp, "set zrange [0:%f]\n", 1.0);
        fprintf(gp, "set zlabel 'T2P2'\n");

        // fprintf(gp, "plot \'%s\' using 2:3 with points pointtype 7 lc rgb 'blue' title \"survivalrateK=%f\",\'%s\' using 1:2:($3-$1):($4-$2) with vectors head filled lc rgb 'blue',\'%s\' using 2:3 with points pointtype 7 lc rgb 'red' title \"finalarrival\"\n", data_file1, K, data_file2, data_file3);
        fprintf(gp, "splot \'%s\' using 2:4:5 with points pointtype 7 lc rgb 'blue' title \"survivalrateK=%f\",\'%s\' using 1:2:3:($4-$1):($5-$2):($6-$3) with vectors head filled lc rgb 'blue',\'%s\' using 2:4:5 with points pointtype 7 lc rgb 'red' title \"finalarrival\"\n", data_file1, K, data_file2, data_file3);

        pclose(gp);

        for (i = 1; i <= 9; i++)
        {
            initP2 = (double)0.1 * i;
            gp = popen("gnuplot -persist", "w");
            fprintf(gp, "set terminal png\n");
            fprintf(gp, "set output 'Twoalleles_hyper_genoportion_new/Two_env_genoport_K_%f_a1_%f_initP2_%f.png'\n", K, a1, initP2);
            fprintf(gp, "set xrange [0:%d]\n", tend);
            fprintf(gp, "set xlabel 't'\n");
            fprintf(gp, "set yrange [0:%f]\n", 1.0);
            fprintf(gp, "set ylabel 'genotype_frequency'\n");
            fprintf(gp, "titles='T1P1 T1P2 T2P1 T2P2'\n");

            fprintf(gp, "set style line 1 lc rgb \"#0000FF\" lw 2\n");
            fprintf(gp, "set style line 2 lc rgb \"#FFFF00\" lw 2\n");
            fprintf(gp, "set style line 3 lc rgb \"#FF0000\" lw 2\n");
            fprintf(gp, "set style line 4 lc rgb \"#000000\" lw 2\n");

            // fprintf(gp, "plot \'%s\' using 2:3 with points pointtype 7 lc rgb 'blue' title \"survivalrateV=%f\",\'%s\' using 1:2:($3-$1):($4-$2) with vectors head filled lc rgb 'blue',\'%s\' using 2:3 with points pointtype 7 lc rgb 'red' title \"finalarrival\"\n", data_file4, V, data_file5, data_file6);
            fprintf(gp, "plot for [j=2:5] \'%s\' every ::%d::%d using 1:j with lines ls (j-1) title word(titles, j-1)\n", data_file7, (i - 1) * (tend + 1), i * (tend + 1) - 1);
            pclose(gp);
        }

        // Map("male", snapshot_file1, initP2, 0);
        // Map("female", snapshot_file1, initP2, 0);
        // Map("male", snapshot_file3, initP2, 0);
        // Map("female", snapshot_file3, initP2, 0);

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
        // free(snapshot_file1);
        // free(snapshot_file2);
        // free(snapshot_file3);
        free(data_file1);
        free(data_file2);
        free(data_file3);
        free(data_file7);
    }
    free(buffer);
    free(genorepo);
    free(rng_states);
}
