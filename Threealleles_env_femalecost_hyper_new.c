/* nearest neighbor interaction */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MT.h"

#define LH 1000 // 10000
#define LV 1000 // 10000
// #define K 0.075 //P3メスのコスト
// #define V 0.074 //P2メスのコスト(0<V<K)
#define u 0.3   //T3オスのコスト0.3
// #define l 0.15  //T2オスのコスト(0<l<u)
// #define a1 3.0    // P2メスがT2オスを選好する倍率
// #define a2 6.0    // P3メスがT3オスを選好する倍率
#define tend 1000 // 4000 80000 10000
#define mapinitP 0.3
#define initialP 3
#define initialT 3

// void Map(const char *sex, const char *filename, double initP, int t)
// {
//     FILE *gp;
//     gp = popen("gnuplot -persist", "w");
//     fprintf(gp, "set term pngcairo size 1000,1000\n");
//     //  fprintf(gp,"set terminal png\n");
//     if (strcmp(sex, "male") == 0)
//     {
//         if (strstr(filename, "stmap"))
//             fprintf(gp, "set output 'Threealleles/ThreeMapSt_%s_env_K_%g_V_%g_l_%f_a1_%f_a2_%f_initP_%g.png'\n", sex, K,V,l,a1,a2, initP);
//         else if (strstr(filename, "intmap"))
//             fprintf(gp, "set output 'Threealleles/ThreeMapInt_%s_env_K_%g_V_%g_l_%f_a1_%f_a2_%f_initP_%g_t_%d.png'\n", sex, K,V,l,a1,a2, initP, t);
//         else if (strstr(filename, "finmap"))
//             fprintf(gp, "set output 'Threealleles/ThreeMapFin_%s_env_K_%g_V_%g_l_%f_a1_%f_a2_%f_initP_%g.png'\n", sex, K,V,l,a1,a2, initP);
//         fprintf(gp, "unset key\n");
//         fprintf(gp, "set size ratio -1\n");
//         fprintf(gp, "set xrange [0:%d]\n", LH - 1);
//         // fprintf(gp,"set xlabel 'T2'\n");
//         fprintf(gp, "set yrange [0:%d]\n", LV - 1);
//         fprintf(gp, "set palette defined(1 \"#FFFFFF\", 2 \"#02befcff\", 3 \"#0000FF\", 4 \"#00FF00\", \
//              5 \"#FFFF00\", 6 \"#FFA500\", 7 \"#FF0000\", 8 \"#FF69B4\", 9 \"#000000\")\n");

//         fprintf(gp,"set cbtics ('T1P1' 1, 'T1P2' 2, 'T1P3' 3, 'T2P1' 4, \
//             'T2P2' 5, 'T2P3' 6, 'T3P1' 7, 'T3P2' 8, 'T3P3' 9)\n");
//         fprintf(gp, "unset autoscale cb\n");
//         fprintf(gp, "set cbrange [1:9]\n");
//         // fprintf(gp,"set multiplot layout 1,2 title 'Genotype map (T×P: 0=T1P1, 1=T1P2, 2=T2P1, 3=T2P2)'\n");
//         fprintf(gp, "set title 'Male map'\n");
//         fprintf(gp, "unset xtics;unset ytics\n");
//         // fprintf(gp,"unset yticks\n");
//         fprintf(gp, "plot \'%s\' using 1:2:3 with image\n", filename);
//     }
//     else if (strcmp(sex, "female") == 0)
//     {
//         if (strstr(filename, "stmap"))
//             fprintf(gp, "set output 'Threealleles/ThreeMapSt_%s_env_K_%g_V_%g_l_%f_a1_%f_a2_%f_initP_%g.png'\n", sex, K,V,l,a1,a2, initP);
//         else if (strstr(filename, "intmap"))
//             fprintf(gp, "set output 'Threealleles/ThreeMapInt_%s_env_K_%g_V_%g_l_%f_a1_%f_a2_%f_initP_%g_t_%d.png'\n", sex, K,V,l,a1,a2, initP, t);
//         else if (strstr(filename, "finmap"))
//             fprintf(gp, "set output 'Threealleles/ThreeMapFin_%s_env_K_%g_V_%g_l_%f_a1_%f_a2_%f_initP_%g.png'\n", sex, K,V,l,a1,a2, initP);
//         fprintf(gp, "unset key\n");
//         fprintf(gp, "set size ratio -1\n");
//         fprintf(gp, "set xrange [0:%d]\n", LH - 1);
//         // fprintf(gp,"set xlabel 'T2'\n");
//         fprintf(gp, "set yrange [0:%d]\n", LV - 1);
//         fprintf(gp, "set palette defined(1 \"#FFFFFF\", 2 \"#02befcff\", 3 \"#0000FF\", 4 \"#00FF00\", \
//              5 \"#FFFF00\", 6 \"#FFA500\", 7 \"#FF0000\", 8 \"#FF69B4\", 9 \"#000000\")\n");

//         fprintf(gp,"set cbtics ('T1P1' 1, 'T1P2' 2, 'T1P3' 3, 'T2P1' 4, \
//             'T2P2' 5, 'T2P3' 6, 'T3P1' 7, 'T3P2' 8, 'T3P3' 9)\n");
//         fprintf(gp, "unset autoscale cb\n");
//         fprintf(gp, "set cbrange [1:9]\n");
//         // fprintf(gp,"set multiplot layout 1,2 title 'Genotype map (T×P: 0=T1P1, 1=T1P2, 2=T2P1, 3=T2P2)'\n");
//         fprintf(gp, "set title 'Male map'\n");
//         fprintf(gp, "unset xtics;unset ytics\n");
//         // fprintf(gp,"unset yticks\n");
//         fprintf(gp, "plot \'%s\' using 1:2:3 with image\n", filename);
//     }
//     else
//         printf("それはだめよ");

//     pclose(gp);
// }

void calc_male_sum(double *sum, int i,int j,int **maleT,int **maleP,int female,double a1,double a2){
    int di[5]={-1,0,1,0,0};
    int dj[5]={0,1,0,-1,0};
    int i2,j2,n,t,p,id;
    double w;
    for(n=0;n<9;n++)sum[n]=0.0;
    for(n=0;n<5;n++){
        i2=(i+di[n]+LH)%LH;
        j2=(j+dj[n]+LV)%LV;
        t=maleT[i2][j2];
        p=maleP[i2][j2];
        w=1.0;
        if(female==2&&t==2)w=a1;
        else if(female==3&&t==3)w=a2;
        if(t==0)printf("T=0");
        id=3*(t-1)+p-1;
        sum[id]+=w;

    }
}

void calc_female_sum(double *sum, int i,int j,int **femaleT,int **femaleP){
    int di[5]={-1,0,1,0,0};
    int dj[5]={0,1,0,-1,0};
    int i2,j2,n,t,p,id;
    for(n=0;n<9;n++)sum[n]=0.0;
    for(n=0;n<5;n++){
        i2=(i+di[n]+LH)%LH;
        j2=(j+dj[n]+LV)%LV;
        t=femaleT[i2][j2];
        p=femaleP[i2][j2];
        if(t==0)printf("T=0");
        id=3*(t-1)+p-1;
        sum[id]+=1.0;

    }
}

void genotype(double *sum,int *sexT0,int *sexP0){
    int n;
    double rnd =genrand_real2();
    double total,acc;
    total=0.0;acc=0.0;
    for(n=0;n<9;n++)total+=sum[n];
    if(total==0.0){
        *sexT0=1;*sexP0=1;
        printf("fault");
        return;
    }
    // 
    if(rnd<sum[0]/total){ *sexT0=1;*sexP0=1;}
    else if(rnd<(sum[0]+sum[1])/total){ *sexT0=1;*sexP0=2;}
    else if(rnd<(sum[0]+sum[1]+sum[2])/total){*sexT0=1;*sexP0=3;}
    else if(rnd<(sum[0]+sum[1]+sum[2]+sum[3])/total){*sexT0=2;*sexP0=1;}
    else if(rnd<(sum[0]+sum[1]+sum[2]+sum[3]+sum[4])/total){*sexT0=2;*sexP0=2;}
    else if(rnd<(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5])/total){*sexT0=2;*sexP0=3;}
    else if(rnd<(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6])/total){*sexT0=3;*sexP0=1;}
    else if(rnd<(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7])/total){*sexT0=3;*sexP0=2;}
    else {*sexT0=3;*sexP0=3;}
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
    int k, k2, i, j, i2, j2, t, ok, x1, x2,a,b,n;
    int maleI, maleJ, femaleI, femaleJ;
    int numMT1, numMT2, numMT3, numMP1, numMP2, numMP3;
    int numFT1, numFT2, numFT3, numFP1, numFP2, numFP3;
    double sum[9],gsum[9];
    double rnd, rnd2,rnd3,rnd4,rnd5, sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9, y1, z1, y2, z2,init,init1;
    double  gsum1, gsum2, gsum3, gsum4, gsum5, gsum6, gsum7, gsum8, gsum9;
    int iK,iV,il,ia2,ia1;
    double K,V,l,a1,a2;
    int maleT0, maleP0, femaleT0, femaleP0, mgenotype, fgenotype,count;
    FILE *gp, *data1, *data2, *data3, *data4, *data5, *data6,*data7;
    FILE *snapshot1, *snapshot2, *snapshot3, *snapshot4, *snapshot5, *snapshot6;
    char *data_file1, *data_file2, *data_file3, *data_file4, *data_file5, *data_file6,*data_file7;
    char *snapshot_file1, *snapshot_file2, *snapshot_file3, *snapshot_file4, *snapshot_file5, *snapshot_file6;
    //data1(T2P2頻度図の書き込み用)
    typedef struct {
    int t;
    double mt2, mp2, initP2;
    } RecordT2P2;

    RecordT2P2 *buffer = malloc(sizeof(RecordT2P2) * 9*(tend + 1));
    RecordT2P2 *buft3p3 = malloc(sizeof(RecordT2P2) * 9*(tend + 1));
    int buf_count = 0;
    //data7(遺伝子頻度の書き込み用)
    typedef struct {
    int t;
    double sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9;
    } Recordgenoport;

    Recordgenoport *genorepo = malloc(sizeof(Recordgenoport) * 9*(tend + 1));
    int geno_count = 0;

    //snapshot2(空間構造の書き込み用)
    typedef struct {
    int i,j,mgenotype,fgenotype;
    
    } Recordmap;

    Recordmap *recomap = malloc(sizeof(Recordmap) * LH*LV);

    int num_threads = omp_get_num_procs(); // 最大利用可能スレッド数（論理コア数）
    mt_state *rng_states = malloc(sizeof(mt_state) * num_threads);
    if (!rng_states) { perror("malloc rng_states"); return 1; }

    for (int tid = 0; tid < num_threads; tid++) {
        init_genrand_mt(&rng_states[tid], 5489UL + (unsigned long)tid * 12345UL);
    }

    omp_set_num_threads(num_threads);
    printf("Using %d threads\n", num_threads);
    fflush(stdout);


    for(iK=1;iK<=3;iK++){
        K=(double)(iK*2-1)*0.01;
        for(iV=1;iV<=3;iV++){
            V=(double)(iV*2-1)*K/6.0;
            for(il=1;il<=3;il++){
                l=(double)(il*2-1)*u/6.0;
                
                for(ia1=2;ia1<=6;ia1++){
                    // if(ia1==3)continue;
                    // else a1=(double)ia1;
                    a1=(double)ia1;
                    for(ia2=2;ia2<=6;ia2++){
                        a2=(double)ia2;
    printf("K V l a2:%f %f %f %f\n",K,V,l,a2);
    RecordT2P2 *buffer = malloc(sizeof(RecordT2P2) * 9*(tend + 1));
    RecordT2P2 *buft3p3 = malloc(sizeof(RecordT2P2) * 9*(tend + 1));
    Recordmap *recomap = malloc(sizeof(Recordmap) * LH*LV);
    Recordgenoport *genorepo = malloc(sizeof(Recordgenoport) * 9*(tend + 1));
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
    sprintf(data_file1, "Three_env_T3P3_K_%f_V_%f_l_%f_a1_%f_a2_%f.dat", K,V,l,a1,a2);
    

    data_file2 = malloc(100);
    sprintf(data_file2, "Three_env_T3P3_flow_K_%f_V_%f_l_%f_a1_%f_a2_%f.dat", K,V,l,a1,a2);

    data_file4 = malloc(100);
    sprintf(data_file4, "Three_env_T2P2_K_%f_V_%f_l_%f_a1_%f_a2_%f.dat", K,V,l,a1,a2);
    

    data_file5 = malloc(100);
    sprintf(data_file5, "Three_env_T2P2_flow_K_%f_V_%f_l_%f_a1_%f_a2_%f.dat", K,V,l,a1,a2);

    data_file7 = malloc(100);
    sprintf(data_file7, "Three_env_genoport_K_%f_V_%f_l_%f_a1_%f_a2_%f.dat", K,V,l,a1,a2);
    

    buf_count=0;
    geno_count=0;

    for (k = 1; k <= 1; k++)
    {
        initT2 = 0.1 * initialT;
        initT3 = (1.0 - initT2) / 2;
        for (k2 = 1; k2 <= 9; k2++) // k2=1; k2<=9; k2++
        {
            // data1 = fopen(data_file1, "a");
            // data4 = fopen(data_file4, "a");
            // data7 = fopen(data_file7, "a");

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
            // printf("0 T2 P2:%f %f \n", (double)numMT2 / (double)(LH * LV), (double)numMP2 / (double)(LH * LV));

            //最初の割合を出力
            sum1 = sum2 = sum3 = sum4 =sum5 = sum6= sum7 = sum8 =sum9= 0.0;
            // data7=fopen(data_file7,"a");
            for(i=0;i<LH;i++){
                for(j=0;j<LV;j++){
                    if ((maleT[i][j] == 1 && maleP[i][j] == 1)||(femaleT[i][j] == 1 && femaleP[i][j] == 1))sum1+=1;
                    else if ((maleT[i][j] == 1 && maleP[i][j] == 2)||(femaleT[i][j] == 1 && femaleP[i][j] == 2))sum2+=1;
                    else if ((maleT[i][j] == 1 && maleP[i][j] == 3)||(femaleT[i][j] == 1 && femaleP[i][j] == 3))sum3+=1;
                    else if ((maleT[i][j] == 2 && maleP[i][j] == 1)||(femaleT[i][j] == 2 && femaleP[i][j] == 1))sum4+=1;
                    else if ((maleT[i][j] == 2 && maleP[i][j] == 2)||(femaleT[i][j] == 2 && femaleP[i][j] == 2))sum5+=1;
                    else if ((maleT[i][j] == 2 && maleP[i][j] == 3)||(femaleT[i][j] == 2 && femaleP[i][j] == 3))sum6+=1;
                    else if ((maleT[i][j] == 3 && maleP[i][j] == 1)||(femaleT[i][j] == 3 && femaleP[i][j] == 1))sum7+=1;
                    else if ((maleT[i][j] == 3 && maleP[i][j] == 2)||(femaleT[i][j] == 3 && femaleP[i][j] == 2))sum8+=1;
                    else if ((maleT[i][j] == 3 && maleP[i][j] == 3)||(femaleT[i][j] == 3 && femaleP[i][j] == 3))sum9+=1;
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
            //             fprintf(snapshot1, "%d\t%d\t%d\t%d\n", i, j, mgenotype, fgenotype);
            //         }
            //     }
            //     fclose(snapshot1);
            // }

            t = 0;
            for(t=0;t<tend;t++)
            {
                if (t % 10 == 0){
                    buffer[buf_count].t = t;
                    buffer[buf_count].mt2 = (double)numMT2 / (double)(LH * LV);
                    buffer[buf_count].mp2 = (double)numMP2 / (double)(LH * LV);
                    buffer[buf_count].initP2 = initP2;
                    buft3p3[buf_count].t = t;
                    buft3p3[buf_count].mt2 = (double)numMT3 / (double)(LH * LV);
                    buft3p3[buf_count].mp2 = (double)numMP3 / (double)(LH * LV);
                    buft3p3[buf_count].initP2 = initP2;
                    buf_count++;
                }
                // --- 追加: dummy 配列を初期化（未書き込み領域を防ぐ） ---
                #pragma omp parallel for collapse(2) schedule(static) \
                    default(none) shared(maleT, maleP, femaleT, femaleP, \
                    maleTdummy, malePdummy, femaleTdummy, femalePdummy) \
                    private(i,j)
                for (i = 0; i < LH; i++) {
                    for (j = 0; j < LV; j++) {
                        maleTdummy[i][j]   = maleT[i][j];
                        malePdummy[i][j]   = maleP[i][j];
                        femaleTdummy[i][j] = femaleT[i][j];
                        femalePdummy[i][j] = femaleP[i][j];
                    }
                }
                // // data1 = fopen(data_file1, "a");
                // if (t % 10 == 0)
                //     fprintf(data1, "%d\t%f\t%f\t%f\n", t, (double)numMT3 / (double)(LH * LV), (double)numMP3 / (double)(LH * LV),initP2);
                // // fclose(data1);

                // // data4 = fopen(data_file4, "a");
                // if (t % 10 == 0)
                //     fprintf(data4, "%d\t%f\t%f\t%f\n", t, (double)numMT2 / (double)(LH * LV), (double)numMP2 / (double)(LH * LV),initP2);
                // // fclose(data4);
                
                #pragma omp parallel for collapse(2) schedule(static) private(sum, gsum, rnd2, rnd4, rnd5, maleT0, maleP0, femaleT0, femaleP0)
                
                for (i = 0; i < LH; i++)
                    for (j = 0; j < LV; j++)
                    {   int tid = omp_get_thread_num();
                        calc_male_sum(sum,i,j,maleT,maleP,femaleP[i][j],a1,a2);
                        // for(n=0;n<9;n++){printf("sum[%d]:%f",n,sum[n]);}
                        // printf("calc_male_sum OK");
                        //メス遺伝のための重みづけ
                        // メスのコスト
                        
                        calc_female_sum(gsum,i,j,femaleT,femaleP);
                        // printf("calc_female_sum OK");

                        rnd4 = genrand_real2_mt(&rng_states[tid]);
                        //オス遺伝
                        if(rnd4<0.5){
                            rnd5 = genrand_real2_mt(&rng_states[tid]);
                            if(rnd5<0.5){
                                do{
                                    rnd2 = genrand_real2_mt(&rng_states[tid]);
                                    genotype(sum,&maleT0,&maleP0);
                                    // printf("genotype OK");
                                    
                                    if(maleT0==1)break;
                                    if(maleT0==2&&rnd2<l)continue;
                                    if(maleT0==3&&rnd2<u)continue;
                                    break;
                                }while(1);//次世代がオスの場合
                                    
                                
                                maleTdummy[i][j] = maleT0;
                                malePdummy[i][j] = maleP0;
                            }else{

                                do{
                                    rnd2 = genrand_real2_mt(&rng_states[tid]);
                                    genotype(sum,&maleT0,&maleP0);
                                        // printf("genotype OK");
                                        
                                    if(maleP0==1)break;
                                    if (maleP0 == 2&&rnd2 < V) continue;
                                    if (maleP0 == 3&&rnd2 < K) continue;
                                    break;    

                                }while(1);
                                 
                                    

                                femaleTdummy[i][j] = maleT0;
                                femalePdummy[i][j] = maleP0;
                                

                            }
                            
                            
                        }
                        else{//メス遺伝
                            rnd5=genrand_real2_mt(&rng_states[tid]);
                            if(rnd5<0.5){//次世代がオス
                                femaleT0=femaleT[i][j];
                                femaleP0=femaleP[i][j];
                                if(femaleT0 != 1){
                                    do{
                                        rnd2 = genrand_real2_mt(&rng_states[tid]);
                                        genotype(sum,&femaleT0,&femaleP0);
                                        // printf("genotype OK");
                                        
                                        if(femaleT0==1)break;
                                        if(femaleT0==2&&rnd2<l)continue;
                                        if(femaleT0==3&&rnd2<u)continue;
                                        break;
                                    }while(1);
                                    
                                
                                }
                                maleTdummy[i][j] = femaleT0;
                                malePdummy[i][j] = femaleP0;

                            }
                            else{//次世代がメス
                                femaleT0=femaleT[i][j];
                                femaleP0=femaleP[i][j];
                                if(femaleP0!=1){
                                   do{
                                        rnd2 = genrand_real2_mt(&rng_states[tid]);
                                        genotype(sum,&femaleT0,&femaleP0);
                                            // printf("genotype OK");
                                            
                                        if(femaleP0==1)break;
                                        if (femaleP0 == 2&&rnd2 < V) continue;
                                        if (femaleP0 == 3&&rnd2 < K) continue;
                                        break;    

                                    }while(1);
                                }
                                femaleTdummy[i][j] = femaleT0;
                                femalePdummy[i][j] = femaleP0;
                            }
                            

                        }

                        


                        

                        



                        
                    }
                // for (i = 0; i < LH; i++)
                //     for (j = 0; j < LV; j++)
                //     {
                //         maleT[i][j] = maleTdummy[i][j];
                //         maleP[i][j] = malePdummy[i][j];
                //         femaleT[i][j] = femaleTdummy[i][j];
                //         femaleP[i][j] = femalePdummy[i][j];
                //     }
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
                

                // 途中の図
                // if (t % 100 == 0 && t < 2000 )//fabs(initP2 - mapinitP) < 1e-12&&t<2000
                // { 
                //     sprintf(snapshot_file2, "Three_env_intmap_t_%d_K_%f_V_%f_l_%f_a1_%f_a2_%f_initP_%g.dat", t, K,V,l,a1,a2, initP2);
                //     snapshot2 = fopen(snapshot_file2, "w");
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
                //             fprintf(snapshot2, "%d\t%d\t%d\t%d\n", i, j, mgenotype, fgenotype);
                //         }
                //     }
                //     fclose(snapshot2);
                //     Map("male", snapshot_file2,initP2, t);
                //     Map("female", snapshot_file2, initP2, t);
                // }
                // else if (t % 200 == 0 && t>2000 && t < tend )//fabs(initP2 - mapinitP) < 1e-12
                // { 
                //     sprintf(snapshot_file2, "Three_env_intmap_t_%d_K_%f_V_%f_l_%f_a1_%f_a2_%f_initP_%g.dat", t, K,V,l,a1,a2 ,initP2);
                //     snapshot2 = fopen(snapshot_file2, "w");
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
                //             fprintf(snapshot2, "%d\t%d\t%d\t%d\n", i, j, mgenotype, fgenotype);
                //         }
                //     }
                //     fclose(snapshot2);
                //     Map("male", snapshot_file2,initP2, t);
                //     Map("female", snapshot_file2, initP2, t);
                // }
            //遺伝子型の割合出力
            sum1 = sum2 = sum3 = sum4 =sum5 = sum6= sum7 = sum8 =sum9= 0.0;
            // data7=fopen(data_file7,"a");
            for(i=0;i<LH;i++){
                for(j=0;j<LV;j++){
                    if ((maleT[i][j] == 1 && maleP[i][j] == 1)||(femaleT[i][j] == 1 && femaleP[i][j] == 1))sum1+=1;
                    else if ((maleT[i][j] == 1 && maleP[i][j] == 2)||(femaleT[i][j] == 1 && femaleP[i][j] == 2))sum2+=1;
                    else if ((maleT[i][j] == 1 && maleP[i][j] == 3)||(femaleT[i][j] == 1 && femaleP[i][j] == 3))sum3+=1;
                    else if ((maleT[i][j] == 2 && maleP[i][j] == 1)||(femaleT[i][j] == 2 && femaleP[i][j] == 1))sum4+=1;
                    else if ((maleT[i][j] == 2 && maleP[i][j] == 2)||(femaleT[i][j] == 2 && femaleP[i][j] == 2))sum5+=1;
                    else if ((maleT[i][j] == 2 && maleP[i][j] == 3)||(femaleT[i][j] == 2 && femaleP[i][j] == 3))sum6+=1;
                    else if ((maleT[i][j] == 3 && maleP[i][j] == 1)||(femaleT[i][j] == 3 && femaleP[i][j] == 1))sum7+=1;
                    else if ((maleT[i][j] == 3 && maleP[i][j] == 2)||(femaleT[i][j] == 3 && femaleP[i][j] == 2))sum8+=1;
                    else if ((maleT[i][j] == 3 && maleP[i][j] == 3)||(femaleT[i][j] == 3 && femaleP[i][j] == 3))sum9+=1;
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
        // fclose(data4);
        // fclose(data7);
        }
    }
    data1 = fopen(data_file1, "w");
    for (int n =0; n < buf_count; n++) {
        fprintf(data1, "%d\t%f\t%f\t%f\n",
         buft3p3[n].t, buft3p3[n].mt2, buft3p3[n].mp2, buft3p3[n].initP2);
    }
    fclose(data1);

    data4 = fopen(data_file4, "w");
    for (int n =0; n < buf_count; n++) {
        fprintf(data4, "%d\t%f\t%f\t%f\n",
        buffer[n].t, buffer[n].mt2, buffer[n].mp2, buffer[n].initP2);
    }
    fclose(data4);  

    data7 = fopen(data_file7, "w");
    for (int n = 0; n < geno_count; n++) {
        fprintf(data7, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
        genorepo[n].t, genorepo[n].sum1, genorepo[n].sum2, genorepo[n].sum3, genorepo[n].sum4,genorepo[n].sum5,\
         genorepo[n].sum6, genorepo[n].sum7, genorepo[n].sum8, genorepo[n].sum9);
    }
    fclose(data7);
    
    printf("Ok\n");
    //T3P3図
    data_file3 = malloc(100);
    sprintf(data_file3, "Three_env_final_T3P3_K_%f_V_%f_l_%f_a1_%f_a2_%f.dat", K,V,l,a1,a2);
    // data_file3 = "final_env.dat";
    gp = fopen(data_file1, "r");
    while (fscanf(gp, "%d %lf %lf %lf", &x1, &y1, &z1,&init) == 4)
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
    if (fscanf(data1, "%d %lf %lf %lf", &x1, &y1, &z1,&init) != 4)
        return 1;
    while (fscanf(data1, "%d %lf %lf %lf",&x2, &y2, &z2, &init1) == 4)
    {   if(fabs(init-init1)<1e-12){
            fprintf(data2, "%lf\t%lf\t%lf\t%lf\n", y1, z1, y2, z2);
            
        }
        x1 = x2;
        y1 = y2;
        z1 = z2;
        init=init1;
    }
    fclose(data1);
    fclose(data2);

    //T2P2図
    data_file6 = malloc(100);
    sprintf(data_file6, "Three_env_final_T2P2_K_%f_V_%f_l_%f_a1_%f_a2_%f.dat", K,V,l,a1,a2);
    // data_file3 = "final_env.dat";
    gp = fopen(data_file4, "r");
    while (fscanf(gp, "%d %lf %lf %lf", &x1, &y1, &z1,&init) == 4)
    {
        if (x1 == (tend - 10))
        {
            data6 = fopen(data_file6, "a");
            fprintf(data6, "%d\t%f\t%f\n",x1, y1, z1);
            fclose(data6);
        }
    }
    fclose(gp);

    data5 = fopen(data_file5, "w");
    data4 = fopen(data_file4, "r");
    if (fscanf(data4, "%d %lf %lf %lf", &x1, &y1, &z1,&init) != 4)
        return 1;
    while (fscanf(data4, "%d %lf %lf %lf", &x2, &y2, &z2,&init1) == 4)
    {
        if(fabs(init-init1)<1e-12){
            fprintf(data5, "%lf\t%lf\t%lf\t%lf\n", y1, z1, y2, z2);
        }
        x1 = x2;
        y1 = y2;
        z1 = z2;
        init=init1;
    }
    fclose(data4);
    fclose(data5);


    gp = popen("gnuplot -persist", "w");
    fprintf(gp, "set terminal png\n");
    fprintf(gp, "set output 'Threealleles_hyper_new/Three_env_T3P3_K_%f_V_%f_l_%f_a1_%f_a2_%f.png'\n", K,V,l,a1,a2);
    fprintf(gp, "set xrange [0:%f]\n", 1.0);
    fprintf(gp, "set xlabel 'T3'\n");
    fprintf(gp, "set yrange [0:%f]\n", 1.0);
    fprintf(gp, "set ylabel 'P3'\n");
    
    fprintf(gp, "plot \'%s\' using 2:3 with points pointtype 7 lc rgb 'blue' title \"survivalrateK=%f\",\'%s\' using 1:2:($3-$1):($4-$2) with vectors head filled lc rgb 'blue',\'%s\' using 2:3 with points pointtype 7 lc rgb 'red' title \"finalarrival\"\n", data_file1, K, data_file2, data_file3);
    
    pclose(gp);

    gp = popen("gnuplot -persist", "w");
    fprintf(gp, "set terminal png\n");
    fprintf(gp, "set output 'Threealleles_hyper_new/Three_env_T2P2_K_%f_V_%f_l_%f_a1_%f_a2_%f.png'\n", K,V,l,a1,a2);
    fprintf(gp, "set xrange [0:%f]\n", 1.0);
    fprintf(gp, "set xlabel 'T2'\n");
    fprintf(gp, "set yrange [0:%f]\n", 1.0);
    fprintf(gp, "set ylabel 'P2'\n");
    
    fprintf(gp, "plot \'%s\' using 2:3 with points pointtype 7 lc rgb 'blue' title \"survivalrateV=%f\",\'%s\' using 1:2:($3-$1):($4-$2) with vectors head filled lc rgb 'blue',\'%s\' using 2:3 with points pointtype 7 lc rgb 'red' title \"finalarrival\"\n", data_file4, V, data_file5, data_file6);
    
    pclose(gp);

    for(i=1;i<=9;i++){
        initP2=(double)0.1*i;
        gp = popen("gnuplot -persist", "w");
        fprintf(gp, "set terminal png\n");
        fprintf(gp, "set output 'Threealleles_hyper_genoportion_new/Three_env_genoport_K_%f_V_%f_l_%f_a1_%f_a2_%f_initP2_%f.png'\n", K,V,l,a1,a2,initP2);
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
        
        //fprintf(gp, "plot \'%s\' using 2:3 with points pointtype 7 lc rgb 'blue' title \"survivalrateV=%f\",\'%s\' using 1:2:($3-$1):($4-$2) with vectors head filled lc rgb 'blue',\'%s\' using 2:3 with points pointtype 7 lc rgb 'red' title \"finalarrival\"\n", data_file4, V, data_file5, data_file6);
        fprintf(gp, "plot for [j=2:10] \'%s\' every ::%d::%d using 1:j with lines ls (j-1) title word(titles, j-1)\n",data_file7,(i - 1) * (tend + 1),i * (tend + 1) - 1);
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
    free(data_file4);
    free(data_file5);
    free(data_file6);
    free(data_file7);
    free(buffer);
    free(buft3p3);
    free(genorepo);
    }}}}}
    free(rng_states);
}
