/* nearest neighbor interaction */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MT.h"

#define LH 1000 // 10000
#define LV 1000 // 10000
#define K 0.075
#define u 0.3
#define a1 3.0    // P2メスがT2オスを選好する倍率
#define a2 3.0    // P3メスがT3オスを選好する倍率
#define tend 1000 // 4000 80000 10000
#define mapinitP 0.5
#define initialP 5
#define initialT 5

void Map(const char *sex, const char *filename, double initP, int t)
{
    FILE *gp;
    gp = popen("gnuplot -persist", "w");
    fprintf(gp, "set term pngcairo size 1000,1000\n");
    //  fprintf(gp,"set terminal png\n");
    if (strcmp(sex, "male") == 0)
    {
        if (strstr(filename, "stmap"))
            fprintf(gp, "set output 'Threealleles/ThreeMapSt_%s_noenv_cost_%g_initP_%g.png'\n", sex, K, initP);
        else if (strstr(filename, "intmap"))
            fprintf(gp, "set output 'Threealleles/ThreeMapInt_%s_noenv_cost_%g_initP_%g_t_%d.png'\n", sex, K, initP, t);
        else if (strstr(filename, "finmap"))
            fprintf(gp, "set output 'Threealleles/ThreeMapFin_%s_noenv_cost_%g_initP_%g.png'\n", sex, K, initP);
        fprintf(gp, "unset key\n");
        fprintf(gp, "set size ratio -1\n");
        fprintf(gp, "set xrange [0:%d]\n", LH - 1);
        // fprintf(gp,"set xlabel 'T2'\n");
        fprintf(gp, "set yrange [0:%d]\n", LV - 1);
        fprintf(gp, "set palette defined(1 'blue',2 'green',3 'orange',4 'red')\n");
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
            fprintf(gp, "set output 'Threealleles/ThreeMapSt_%s_noenv_cost_%g_initP_%g.png'\n", sex, K, initP);
        else if (strstr(filename, "intmap"))
            fprintf(gp, "set output 'Threealleles/ThreeMapInt_%s_noenv_cost_%g_initP_%g_t_%d.png'\n", sex, K, initP, t);
        else if (strstr(filename, "finmap"))
            fprintf(gp, "set output 'Threealleles/ThreeMapFin_%s_noenv_cost_%g_initP_%g.png'\n", sex, K, initP);
        fprintf(gp, "unset key\n");
        fprintf(gp, "set size ratio -1\n");
        fprintf(gp, "set xrange [0:%d]\n", LH - 1);
        // fprintf(gp,"set xlabel 'T2'\n");
        fprintf(gp, "set yrange [0:%d]\n", LV - 1);
        fprintf(gp, "set palette defined(1 'blue',2 'green',3 'orange',4 'red')\n");
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
    double initT2, initP2, initT3, initP3;
    int k, k2, i, j, i2, j2, t, ok, x1, x2;
    int maleI, maleJ, femaleI, femaleJ;
    int numMT1, numMT2, numMT3, numMP1, numMP2, numMP3;
    int numFT1, numFT2, numFT3, numFP1, numFP2, numFP3;
    double rnd, rnd2, sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9, y1, z1, y2, z2;
    int maleT0, maleP0, femaleT0, femaleP0, mgenotype, fgenotype;
    FILE *gp, *data1, *data2, *data3, *data4, *data5, *data6;
    FILE *snapshot1, *snapshot2, *snapshot3, *snapshot4, *snapshot5, *snapshot6;
    char *data_file1, *data_file2, *data_file3, *data_file4, *data_file5, *data_file6;
    char *snapshot_file1, *snapshot_file2, *snapshot_file3, *snapshot_file4, *snapshot_file5, *snapshot_file6;

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

    snapshot_file1 = malloc(100);
    sprintf(snapshot_file1, "Three_noenv_stmap_%f_initP_%g.dat", K, mapinitP);

    snapshot_file2 = malloc(100);
    // sprintf(snapshot2, "efcs_intmap_%f.dat", K);

    snapshot_file3 = malloc(100);
    sprintf(snapshot_file3, "Three_noenv_finmap_%f_initP_%g.dat", K, mapinitP);

    data_file1 = malloc(100);
    sprintf(data_file1, "Three_noenv_T3P3_%f.dat", K);

    data_file2 = malloc(100);
    sprintf(data_file2, "Three_noenv_T3P3_flow_%f.dat", K);

    data_file4 = malloc(100);
    sprintf(data_file4, "Three_noenv_T2P2_%f.dat", K);

    data_file5 = malloc(100);
    sprintf(data_file5, "Three_noenv_T2P2_flow_%f.dat", K);

    for (k = 1; k <= 1; k++)
    {
        initT2 = 0.1 * initialT;
        initT3 = (1.0 - initT2) / 2;
        for (k2 = initialP; k2 <= initialP; k2++) // k2=1; k2<=9; k2++
        {
            initP2 = 0.1 * k2;
            initP3 = (1.0 - initP2) / 2;
            for (i = 0; i < LH; i++)
            {
                for (j = 0; j < LV; j++)
                {
                    // 通常
                    rnd = genrand_real2();
                    if (rnd < initT3)
                        maleT[i][j] = 3;
                    else if (rnd < initT3 + initT2)
                        maleT[i][j] = 2;
                    else
                        maleT[i][j] = 1;
                    rnd = genrand_real2();
                    if (rnd < initP3)
                        maleP[i][j] = 3;
                    else if (rnd < initP3 + initP2)
                        maleP[i][j] = 2;
                    else
                        maleP[i][j] = 1;
                    rnd = genrand_real2();
                    if (rnd < initT3)
                        femaleT[i][j] = 3;
                    else if (rnd < initT3 + initT2)
                        femaleT[i][j] = 2;
                    else
                        femaleT[i][j] = 1;
                    rnd = genrand_real2();
                    if (rnd < initP3)
                        femaleP[i][j] = 3;
                    else if (rnd < initP3 + initP2)
                        femaleP[i][j] = 2;
                    else
                        femaleP[i][j] = 1;

                    // T1P1vsT2P1

                    // T2P1vsT2P2

                    // T1P1vsT2P2
                }
            }

            numMT1 = numMT2 = numMT3 = numFT1 = numFT2 = numFT3 = numMP1 = numMP2 = numMP3 = numFP1 = numFP2 = numFP3 = 0;
            for (i = 0; i < LH; i++)
            {
                for (j = 0; j < LV; j++)
                {
                    if (maleT[i][j] == 1)
                        numMT1++;
                    else if (maleT[i][j] == 2)
                        numMT2++;
                    else if (maleT[i][j] == 3)
                        numMT3++;
                    if (femaleT[i][j] == 1)
                        numFT1++;
                    else if (femaleT[i][j] == 2)
                        numFT2++;
                    else if (femaleT[i][j] == 3)
                        numFT3++;
                    if (maleP[i][j] == 1)
                        numMP1++;
                    else if (maleP[i][j] == 2)
                        numMP2++;
                    else if (maleP[i][j] == 3)
                        numMP3++;
                    if (femaleP[i][j] == 1)
                        numFP1++;
                    else if (femaleP[i][j] == 2)
                        numFP2++;
                    else if (femaleP[i][j] == 3)
                        numFP3++;
                }
            }
            printf("0 T3P3:%f %f \n", (double)numMT3 / (double)(LH * LV), (double)numMP3 / (double)(LH * LV));
            // 初期の図
            if (initP2 == mapinitP)
            {
                snapshot1 = fopen(snapshot_file1, "w");
                mgenotype = fgenotype = 0;
                for (i = 0; i < LH; i++)
                {
                    for (j = 0; j < LV; j++)
                    {
                        if (maleT[i][j] == 1 && maleP[i][j] == 1)
                            mgenotype = 1;
                        else if (maleT[i][j] == 1 && maleP[i][j] == 2)
                            mgenotype = 2;
                        else if (maleT[i][j] == 1 && maleP[i][j] == 2)
                            mgenotype = 3;
                        else if (maleT[i][j] == 2 && maleP[i][j] == 1)
                            mgenotype = 4;
                        else if (maleT[i][j] == 2 && maleP[i][j] == 2)
                            mgenotype = 5;
                        else if (maleT[i][j] == 2 && maleP[i][j] == 3)
                            mgenotype = 6;
                        else if (maleT[i][j] == 3 && maleP[i][j] == 1)
                            mgenotype = 7;
                        else if (maleT[i][j] == 3 && maleP[i][j] == 2)
                            mgenotype = 8;
                        else if (maleT[i][j] == 3 && maleP[i][j] == 3)
                            mgenotype = 9;
                        if (femaleT[i][j] == 1 && femaleP[i][j] == 1)
                            fgenotype = 1;
                        else if (femaleT[i][j] == 1 && femaleP[i][j] == 2)
                            fgenotype = 2;
                        else if (femaleT[i][j] == 1 && femaleP[i][j] == 3)
                            fgenotype = 3;
                        else if (femaleT[i][j] == 2 && femaleP[i][j] == 1)
                            fgenotype = 4;
                        else if (femaleT[i][j] == 2 && femaleP[i][j] == 2)
                            fgenotype = 5;
                        else if (femaleT[i][j] == 2 && femaleP[i][j] == 3)
                            fgenotype = 6;
                        else if (femaleT[i][j] == 3 && femaleP[i][j] == 1)
                            fgenotype = 7;
                        else if (femaleT[i][j] == 3 && femaleP[i][j] == 2)
                            fgenotype = 8;
                        else if (femaleT[i][j] == 3 && femaleP[i][j] == 3)
                            fgenotype = 9;
                        fprintf(snapshot1, "%d\t%d\t%d\t%d\n", i, j, mgenotype, fgenotype);
                    }
                }
                fclose(snapshot1);
            }

            t = 0;
            while (t < tend)
            {
                data1 = fopen(data_file1, "a");
                if (t % 10 == 0)
                    fprintf(data1, "%d\t%f\t%f\n", t, (double)numMT3 / (double)(LH * LV), (double)numMP3 / (double)(LH * LV));
                fclose(data1);

                data4 = fopen(data_file1, "a");
                if (t % 10 == 0)
                    fprintf(data4, "%d\t%f\t%f\n", t, (double)numMT2 / (double)(LH * LV), (double)numMP2 / (double)(LH * LV));
                fclose(data4);

                for (i = 0; i < LH; i++)
                    for (j = 0; j < LV; j++)
                    {
                        sum1 = sum2 = sum3 = sum4 = sum5 = sum6 = sum7 = sum8 = sum9 = 0.0;
                        if (maleT[i][j] == 1 && maleP[i][j] == 1)
                            sum1 += 1.0;
                        else if (maleT[i][j] == 1 && maleP[i][j] == 2)
                            sum2 += 1.0;
                        else if (maleT[i][j] == 1 && maleP[i][j] == 3)
                            sum3 += 1.0;
                        else if (maleT[i][j] == 2 && maleP[i][j] == 1)
                            sum4 += a1;
                        else if (maleT[i][j] == 2 && maleP[i][j] == 2)
                            sum5 += a1;
                        else if (maleT[i][j] == 2 && maleP[i][j] == 3)
                            sum6 += a1;
                        else if (maleT[i][j] == 3 && maleP[i][j] == 1)
                            sum7 += a2;
                        else if (maleT[i][j] == 3 && maleP[i][j] == 2)
                            sum8 += a2;
                        else if (maleT[i][j] == 3 && maleP[i][j] == 3)
                            sum9 += a2;
                    }
                // ここまでやった2025/10/20
                for (i = 0; i < LH; i++)
                    for (j = 0; j < LV; j++)
                    {
                        if (femaleP[i][j] == 1)
                        {
                            // sum1 = sum2 = sum3 = sum4 =sum5 = sum6= sum7 = sum8 =sum9= 0.0;
                            // i2 = (i - 1 + LH) % LH;
                            // j2 = j;
                            // if (maleT[i2][j2] == 1 && maleP[i2][j2] == 1)
                            //     sum1 += 1.0;
                            // else if (maleT[i2][j2] == 1 && maleP[i2][j2] == 2)
                            //     sum2 += 1.0;
                            // else if (maleT[i2][j2] == 2 && maleP[i2][j2] == 1)
                            //     sum3 += 1.0;
                            // else if (maleT[i2][j2] == 2 && maleP[i2][j2] == 2)
                            //     sum4 += 1.0;
                            // i2 = (i + 1 + LH) % LH;
                            // j2 = j;
                            // if (maleT[i2][j2] == 1 && maleP[i2][j2] == 1)
                            //     sum1 += 1.0;
                            // else if (maleT[i2][j2] == 1 && maleP[i2][j2] == 2)
                            //     sum2 += 1.0;
                            // else if (maleT[i2][j2] == 2 && maleP[i2][j2] == 1)
                            //     sum3 += 1.0;
                            // else if (maleT[i2][j2] == 2 && maleP[i2][j2] == 2)
                            //     sum4 += 1.0;
                            // i2 = i;
                            // j2 = (j - 1 + LV) % LV;
                            // if (maleT[i2][j2] == 1 && maleP[i2][j2] == 1)
                            //     sum1 += 1.0;
                            // else if (maleT[i2][j2] == 1 && maleP[i2][j2] == 2)
                            //     sum2 += 1.0;
                            // else if (maleT[i2][j2] == 2 && maleP[i2][j2] == 1)
                            //     sum3 += 1.0;
                            // else if (maleT[i2][j2] == 2 && maleP[i2][j2] == 2)
                            //     sum4 += 1.0;
                            // i2 = i;
                            // j2 = (j + 1 + LV) % LV;
                            // if (maleT[i2][j2] == 1 && maleP[i2][j2] == 1)
                            //     sum1 += 1.0;
                            // else if (maleT[i2][j2] == 1 && maleP[i2][j2] == 2)
                            //     sum2 += 1.0;
                            // else if (maleT[i2][j2] == 2 && maleP[i2][j2] == 1)
                            //     sum3 += 1.0;
                            // else if (maleT[i2][j2] == 2 && maleP[i2][j2] == 2)
                            //     sum4 += 1.0;

                            // // 新たに挿入
                            // i2 = i;
                            // j2 = j;
                            // if (maleT[i2][j2] == 1 && maleP[i2][j2] == 1)
                            //     sum1 += 1.0;
                            // else if (maleT[i2][j2] == 1 && maleP[i2][j2] == 2)
                            //     sum2 += 1.0;
                            // else if (maleT[i2][j2] == 2 && maleP[i2][j2] == 1)
                            //     sum3 += 1.0;
                            // else if (maleT[i2][j2] == 2 && maleP[i2][j2] == 2)
                            //     sum4 += 1.0;
                        }
                        else
                        {
                            // sum1 = sum2 = sum3 = sum4 = 0.0;
                            // i2 = (i - 1 + LH) % LH;
                            // j2 = j;
                            // if (maleT[i2][j2] == 1 && maleP[i2][j2] == 1)
                            //     sum1 += 1.0;
                            // else if (maleT[i2][j2] == 1 && maleP[i2][j2] == 2)
                            //     sum2 += 1.0;
                            // else if (maleT[i2][j2] == 2 && maleP[i2][j2] == 1)
                            //     sum3 += a;
                            // else if (maleT[i2][j2] == 2 && maleP[i2][j2] == 2)
                            //     sum4 += a;
                            // i2 = (i + 1 + LH) % LH;
                            // j2 = j;
                            // if (maleT[i2][j2] == 1 && maleP[i2][j2] == 1)
                            //     sum1 += 1.0;
                            // else if (maleT[i2][j2] == 1 && maleP[i2][j2] == 2)
                            //     sum2 += 1.0;
                            // else if (maleT[i2][j2] == 2 && maleP[i2][j2] == 1)
                            //     sum3 += a;
                            // else if (maleT[i2][j2] == 2 && maleP[i2][j2] == 2)
                            //     sum4 += a;
                            // i2 = i;
                            // j2 = (j - 1 + LV) % LV;
                            // if (maleT[i2][j2] == 1 && maleP[i2][j2] == 1)
                            //     sum1 += 1.0;
                            // else if (maleT[i2][j2] == 1 && maleP[i2][j2] == 2)
                            //     sum2 += 1.0;
                            // else if (maleT[i2][j2] == 2 && maleP[i2][j2] == 1)
                            //     sum3 += a;
                            // else if (maleT[i2][j2] == 2 && maleP[i2][j2] == 2)
                            //     sum4 += a;
                            // i2 = i;
                            // j2 = (j + 1 + LV) % LV;
                            // if (maleT[i2][j2] == 1 && maleP[i2][j2] == 1)
                            //     sum1 += 1.0;
                            // else if (maleT[i2][j2] == 1 && maleP[i2][j2] == 2)
                            //     sum2 += 1.0;
                            // else if (maleT[i2][j2] == 2 && maleP[i2][j2] == 1)
                            //     sum3 += a;
                            // else if (maleT[i2][j2] == 2 && maleP[i2][j2] == 2)
                            //     sum4 += a;
                            // // 新たに挿入
                            // i2 = i;
                            // j2 = j;
                            // if (maleT[i2][j2] == 1 && maleP[i2][j2] == 1)
                            //     sum1 += 1.0;
                            // else if (maleT[i2][j2] == 1 && maleP[i2][j2] == 2)
                            //     sum2 += 1.0;
                            // else if (maleT[i2][j2] == 2 && maleP[i2][j2] == 1)
                            //     sum3 += a;
                            // else if (maleT[i2][j2] == 2 && maleP[i2][j2] == 2)
                            //     sum4 += a;
                        }

                        do
                        {
                            rnd = genrand_real2();
                            if (rnd < sum1 / (sum1 + sum2 + sum3 + sum4))
                            {
                                maleT0 = 1;
                                maleP0 = 1;
                            }
                            else if (rnd < (sum1 + sum2) / (sum1 + sum2 + sum3 + sum4))
                            {
                                maleT0 = 1;
                                maleP0 = 2;
                            }
                            else if (rnd < (sum1 + sum2 + sum3) / (sum1 + sum2 + sum3 + sum4))
                            {
                                maleT0 = 2;
                                maleP0 = 1;
                            }
                            else
                            {
                                maleT0 = 2;
                                maleP0 = 2;
                            }
                            if (maleT0 == 1)
                                rnd2 = 1.0;
                            else
                                rnd2 = genrand_real2();
                        } while (rnd2 < u);
                        // メスのコスト
                        sum1 = sum2 = sum3 = sum4 = 0.0;
                        i2 = (i - 1 + LH) % LH;
                        j2 = j;
                        if (femaleT[i2][j2] == 1 && femaleP[i2][j2] == 1)
                            sum1 += 1.0;
                        else if (femaleT[i2][j2] == 1 && femaleP[i2][j2] == 2)
                            sum2 += 1.0;
                        else if (femaleT[i2][j2] == 2 && femaleP[i2][j2] == 1)
                            sum3 += 1.0;
                        else if (femaleT[i2][j2] == 2 && femaleP[i2][j2] == 2)
                            sum4 += 1.0;
                        i2 = (i + 1 + LH) % LH;
                        j2 = j;
                        if (femaleT[i2][j2] == 1 && femaleP[i2][j2] == 1)
                            sum1 += 1.0;
                        else if (femaleT[i2][j2] == 1 && femaleP[i2][j2] == 2)
                            sum2 += 1.0;
                        else if (femaleT[i2][j2] == 2 && femaleP[i2][j2] == 1)
                            sum3 += 1.0;
                        else if (femaleT[i2][j2] == 2 && femaleP[i2][j2] == 2)
                            sum4 += 1.0;
                        i2 = i;
                        j2 = (j - 1 + LV) % LV;
                        if (femaleT[i2][j2] == 1 && femaleP[i2][j2] == 1)
                            sum1 += 1.0;
                        else if (femaleT[i2][j2] == 1 && femaleP[i2][j2] == 2)
                            sum2 += 1.0;
                        else if (femaleT[i2][j2] == 2 && femaleP[i2][j2] == 1)
                            sum3 += 1.0;
                        else if (femaleT[i2][j2] == 2 && femaleP[i2][j2] == 2)
                            sum4 += 1.0;
                        i2 = i;
                        j2 = (j + 1 + LV) % LV;
                        if (femaleT[i2][j2] == 1 && femaleP[i2][j2] == 1)
                            sum1 += 1.0;
                        else if (femaleT[i2][j2] == 1 && femaleP[i2][j2] == 2)
                            sum2 += 1.0;
                        else if (femaleT[i2][j2] == 2 && femaleP[i2][j2] == 1)
                            sum3 += 1.0;
                        else if (femaleT[i2][j2] == 2 && femaleP[i2][j2] == 2)
                            sum4 += 1.0;
                        // 新たに挿入
                        i2 = i;
                        j2 = j;
                        if (femaleT[i2][j2] == 1 && femaleP[i2][j2] == 1)
                            sum1 += 1.0;
                        else if (femaleT[i2][j2] == 1 && femaleP[i2][j2] == 2)
                            sum2 += 1.0;
                        else if (femaleT[i2][j2] == 2 && femaleP[i2][j2] == 1)
                            sum3 += 1.0;
                        else if (femaleT[i2][j2] == 2 && femaleP[i2][j2] == 2)
                            sum4 += 1.0;

                        do
                        {
                            rnd = genrand_real2();
                            if (rnd < sum1 / (sum1 + sum2 + sum3 + sum4))
                            {
                                femaleT0 = 1;
                                femaleP0 = 1;
                            }
                            else if (rnd < (sum1 + sum2) / (sum1 + sum2 + sum3 + sum4))
                            {
                                femaleT0 = 1;
                                femaleP0 = 2;
                            }
                            else if (rnd < (sum1 + sum2 + sum3) / (sum1 + sum2 + sum3 + sum4))
                            {
                                femaleT0 = 2;
                                femaleP0 = 1;
                            }
                            else
                            {
                                femaleT0 = 2;
                                femaleP0 = 2;
                            }
                            if (femaleP0 == 1)
                                rnd2 = 1.0;
                            else
                                rnd2 = genrand_real2();
                        } while (rnd2 < K);

                        rnd = genrand_real2();
                        if (rnd < 0.5)
                        {
                            maleTdummy[i][j] = maleT0;
                            malePdummy[i][j] = maleP0;
                        }
                        else
                        {
                            maleTdummy[i][j] = femaleT0;
                            malePdummy[i][j] = femaleP0;
                        }
                        rnd = genrand_real2();
                        if (rnd < 0.5)
                        {
                            femaleTdummy[i][j] = maleT0;
                            femalePdummy[i][j] = maleP0;
                        }
                        else
                        {
                            femaleTdummy[i][j] = femaleT0;
                            femalePdummy[i][j] = femaleP0;
                        }
                    }
                for (i = 0; i < LH; i++)
                    for (j = 0; j < LV; j++)
                    {
                        maleT[i][j] = maleTdummy[i][j];
                        maleP[i][j] = malePdummy[i][j];
                        femaleT[i][j] = femaleTdummy[i][j];
                        femaleP[i][j] = femalePdummy[i][j];
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
                t++;

                // 途中の図
                if (t % 10 == 0 && t < tend && fabs(initP2 - mapinitP) < 1e-12)
                { // t<tend
                    sprintf(snapshot_file2, "Three_noenv_intmap_t_%d_cost_%f_initP_%g.dat", t, K, mapinitP);
                    snapshot2 = fopen(snapshot_file2, "w");
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
                            fprintf(snapshot2, "%d\t%d\t%d\t%d\n", i, j, mgenotype, fgenotype);
                        }
                    }
                    fclose(snapshot2);
                    Map("male", snapshot_file2, mapinitP, t);
                    Map("female", snapshot_file2, mapinitP, t);
                }
            }
            // 最後の図
            if (initP2 == mapinitP)
            {
                snapshot3 = fopen(snapshot_file3, "w");
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
                        fprintf(snapshot3, "%d\t%d\t%d\t%d\n", i, j, mgenotype, fgenotype);
                    }
                }
                fclose(snapshot3);
            }
        }
    }
    data_file3 = malloc(100);
    sprintf(data_file3, "Three_noenv_final_%f.dat", K);
    // data_file3 = "final_env.dat";
    gp = fopen(data_file1, "r");
    while (fscanf(gp, "%d %lf %lf", &x1, &y1, &z1) == 3)
    {
        if (x1 == (tend - 10))
        {
            data3 = fopen(data_file3, "a");
            fprintf(data3, "%d\t%f\t%f\n", x1, y1, z1);
            fclose(data3);
        }
    }
    fclose(gp);

    data2 = fopen(data_file2, "w");
    data1 = fopen(data_file1, "r");
    if (fscanf(data1, "%d %lf %lf", &x1, &y1, &z1) != 3)
        return 1;
    while (fscanf(data1, "%d %lf %lf", &x2, &y2, &z2) == 3)
    {
        fprintf(data2, "%lf\t%lf\t%lf\t%lf\n", y1, z1, y2, z2);
        x1 = x2;
        y1 = y2;
        z1 = z2;
    }
    fclose(data1);
    fclose(data2);

    gp = popen("gnuplot -persist", "w");
    fprintf(gp, "set terminal png\n");
    fprintf(gp, "set output 'Threealleles/Three_noenv_%f.png'\n", K);
    fprintf(gp, "set xrange [0:%f]\n", 1.0);
    fprintf(gp, "set xlabel 'T2'\n");
    fprintf(gp, "set yrange [0:%f]\n", 1.0);
    fprintf(gp, "set ylabel 'P2'\n");
    // fprintf(gp, "u=%f\n", u);
    // fprintf(gp, "a=%f\n", a);
    // fprintf(gp, "f(x)=u*(x*(a*(1-u)-1)+1)/((1-u)*(a-1))\n");
    // fprintf(gp,"plot \'%s\' using 1:2 with lines linetype 1 title \"a= %f  S \",\'%s\' using 1:3 with lines linetype 3 title \"I \",\'%s\' using 1:4 with lines linetype 4 title \"0\"\n",data_file3,A_list[h],data_file3,data_file3);
    //  fprintf(gp,"plot \'%s\' using 1:2 with points pointtype 7 lc rgb 'blue' title \"point \"\n",data_file1);
    fprintf(gp, "plot \'%s\' using 2:3 with points pointtype 7 lc rgb 'blue' title \"survivalrateK=%f\",\'%s\' using 1:2:3:4 with vectors head filled lc rgb 'blue',\'%s\' using 2:3 with points pointtype 7 lc rgb 'red' title \"finalarrival\"\n", data_file1, K, data_file2, data_file3);
    // fprintf(gp,"plot f(x) with lines linetype 1 lc rgb 'red' title \"equilibria line\"\n");
    pclose(gp);

    Map("male", snapshot_file1, mapinitP, 0);
    Map("female", snapshot_file1, mapinitP, 0);
    // Map("male",snapshot_file2);
    // Map("female",snapshot_file2);
    Map("male", snapshot_file3, mapinitP, 0);
    Map("female", snapshot_file3, mapinitP, 0);
    // gp=popen("gnuplot -persist","w");
    // // fprintf(gp,"set term pngcairo size 2000,1000\n");
    //  fprintf(gp,"set terminal png\n");
    // fprintf(gp,"set output 'Map_male_env_cost_self_%f.png'\n",K);
    // fprintf(gp,"unset key\n");
    // fprintf(gp,"set size ratio -1\n");
    // fprintf(gp,"set xrange [0:%d]\n",LH-1);
    // // fprintf(gp,"set xlabel 'T2'\n");
    // fprintf(gp,"set yrange [0:%d]\n",LV-1);
    // fprintf(gp,"set palette defined(1 'blue',2 'green',3 'orange',4 'red')\n");
    // // fprintf(gp,"set multiplot layout 1,2 title 'Genotype map (T×P: 0=T1P1, 1=T1P2, 2=T2P1, 3=T2P2)'\n");
    // fprintf(gp,"set title 'Male map at first'\n");
    // fprintf(gp,"unset xtics;unset ytics\n");
    // // fprintf(gp,"unset yticks\n");
    // fprintf(gp,"plot \'%s\' using 1:2:3 with image\n",snapshot_file1);
    // // fprintf(gp,"set ylabel 'P2'\n");
    // // fprintf(gp,"u=%f\n",u);
    // // fprintf(gp,"a=%f\n",a);
    // // fprintf(gp,"f(x)=u*(x*(a*(1-u)-1)+1)/((1-u)*(a-1))\n");
    // //fprintf(gp,"plot \'%s\' using 1:2 with lines linetype 1 title \"a= %f  S \",\'%s\' using 1:3 with lines linetype 3 title \"I \",\'%s\' using 1:4 with lines linetype 4 title \"0\"\n",data_file3,A_list[h],data_file3,data_file3);
    // // fprintf(gp,"plot \'%s\' using 1:2 with points pointtype 7 lc rgb 'blue' title \"point \"\n",data_file1);
    // // fprintf(gp,"plot \'%s\' using 2:3 with points pointtype 7 lc rgb 'blue' title \"survivalrateK=%f\",\'%s\' using 2:3 with points pointtype 7 lc rgb 'red' title \"finalarrival\",f(x) with lines linetype 1 lc rgb 'red' title \"equilibria line\"\n",data_file1,K,data_file3);
    // // fprintf(gp,"plot f(x) with lines linetype 1 lc rgb 'red' title \"equilibria line\"\n");
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
    free(snapshot_file1);
    free(snapshot_file2);
    free(snapshot_file3);
    free(data_file1);
    free(data_file2);
    free(data_file3);
}
