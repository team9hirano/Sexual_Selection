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
    int x,x1,x2,iK,iV,il,ia1,ia2,i,j,n,i2,j2,t,p;
    double K,y, z,V, l,a1,a2,sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9,init;
    double gsum1, gsum2, gsum3, gsum4, gsum5, gsum6, gsum7, gsum8, gsum9,init1;
    K=0.0;V=K/2;a2=3.0;
    FILE *gp, *data1, *data2, *data3;
    char *data_file2, *data_file3;// *data_file1,
    for(iK=0;iK<=1;iK++){//iK=1;iK<=3;iK++
        K=0.115+(double)0.01*(double)iK;
        // K=0.11+(double)iK*0.001;
        for(iV=2;iV<=2;iV++){
            V=(double)(iV*2-1)*K/6.0;
            // V=0.20+(double)iV*0.01;
            for(il=2;il<=2;il++){
                l=(double)(il*2-1)*u/6.0;
                
                for(ia1=3;ia1<=3;ia1++){
                    // if(ia1==3)continue;
                    // else a1=(double)ia1;
                    a1=(double)ia1;
                    for(ia2=3;ia2<=3;ia2++){
                        a2=(double)ia2;
    printf("K V l a2:%f %f %f %f\n",K,V,l,a2);

   char data_file1[256];
    snprintf(data_file1, sizeof(data_file1),
        "Three_env_T3P3_K_%f_V_%f_l_0.150000_a1_3.000000_a2_3.000000.dat",
        K, V
    );

    FILE *fp = fopen(data_file1, "r");
    if (!fp) {
        printf("ファイル読めん: %s\n", data_file1);
        continue;
    }

    // data_file1= f("Three_env_T3P3_K_%f_V_%f_l_%f_a1_%f_a2_%f.dat", K,V,l,a1,a2);
    

    data_file2 = malloc(100);
    sprintf(data_file2, "Three_env_T3P3_flow_K_%f_V_%f_l_%f_a1_%f_a2_%f.dat", K,V,l,a1,a2);

    data_file3 = malloc(100);
    sprintf(data_file3, "Three_env_final_T2P2_T3P3_K_%f_V_%f_l_%f_a1_%f_a2_%f.dat", K,V,l,a1,a2);
    // data_file3 = "final_env.dat";

    data2 = fopen(data_file2, "w");
    data1 = fopen(data_file1, "r");
    if (fscanf(data1, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &x1, &sum1, &sum2,&sum3,&sum4,&sum5,\
        &sum6,&sum7,&sum8,&sum9,&init) != 11)
        return 1;
    while (fscanf(data1, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &x2, &gsum1, &gsum2,&gsum3,&gsum4,&gsum5,\
        &gsum6,&gsum7,&gsum8,&gsum9,&init1) == 11)
    {   if(fabs(init-init1)<1e-12){
            fprintf(data2, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n"\
                , sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9,\
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
        init=init1;
    }
    fclose(data2);
    fclose(data1);
    gp = fopen(data_file1, "r");
    data3 = fopen(data_file3, "w");
    while (fscanf(gp, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &x1, &sum1, &sum2,&sum3,&sum4,&sum5,\
        &sum6,&sum7,&sum8,&sum9,&init) == 11)
    {
        if (x1 == (tend - 10))
        {
            fprintf(data3, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", x1, sum1, \
                sum2, sum3, sum4, sum5, sum6, sum7, sum8,sum9);
        }
    }
    fclose(data3);
    fclose(gp);

    //T3P3-T2P2図示
    gp = popen("gnuplot -persist", "w");
    fprintf(gp, "set terminal png\n");
    fprintf(gp, "set term pngcairo size 1000,700\n");
    fprintf(gp, "set output 'Genotype_Threealleles/Three_env_T3P3_T2P2_K_%f_V_%f_l_%f_a1_%f_a2_%f.png'\n", K,V,l,a1,a2);
    fprintf(gp, "set xrange [0:%f]\n", 1.0);
    fprintf(gp, "set xlabel 'T3P3'\n");
    fprintf(gp, "set yrange [0:%f]\n", 1.0);
    fprintf(gp, "set ylabel 'T2P2'\n");
    
    fprintf(gp, "plot \'%s\' using 10:6 with points pointtype 7 lc rgb 'blue'\
         title \"survivalrateK=%f\",\'%s\' using 9:5:($18-$9):($14-$5) with vectors head filled lc rgb 'blue'\
         ,\'%s\' using 10:6 with points pointtype 7 lc rgb 'red' title \"finalarrival\"\n"\
         , data_file1, K, data_file2, data_file3);
    
    pclose(gp);

    //T2P1-T2P2図示
    gp = popen("gnuplot -persist", "w");
    fprintf(gp, "set terminal png\n");
    fprintf(gp, "set term pngcairo size 1000,700\n");
    fprintf(gp, "set output 'Genotype_Threealleles/Three_env_T2P1_T2P2_K_%f_V_%f_l_%f_a1_%f_a2_%f.png'\n", K,V,l,a1,a2);
    fprintf(gp, "set xrange [0:%f]\n", 1.0);
    fprintf(gp, "set xlabel 'T2P1'\n");
    fprintf(gp, "set yrange [0:%f]\n", 1.0);
    fprintf(gp, "set ylabel 'T2P2'\n");
    
    fprintf(gp, "plot \'%s\' using 5:6 with points pointtype 7 lc rgb 'blue'\
         title \"survivalrateK=%f\",\'%s\' using 4:5:($13-$4):($14-$5) with vectors head filled lc rgb 'blue'\
         ,\'%s\' using 5:6 with points pointtype 7 lc rgb 'red' title \"finalarrival\"\n"\
         , data_file1, K, data_file2, data_file3);
    
    pclose(gp);
}}}}}
}
