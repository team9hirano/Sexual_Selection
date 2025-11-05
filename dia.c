/* nearest neighbor interaction */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "MT.h"

#define LH	1000	//10000
#define LV	1000	//10000
#define K  0.00
#define u	0.3
#define a	3.0
#define tend	20000 //2000

int main(void)
{
    int x;
    double y,z;
    FILE *gp,*data1,*data2,*data3;
    char *data_file1,*data_file2,*data_file3;
    data_file1="env_female_cost_self_0.000000.dat";
    data_file3 = "finalenv_cost_self_0.000000.dat";
    // gp=fopen(data_file1,"r");
    // while(fscanf(gp,"%d %lf %lf",&x,&y,&z)==3){
    //       if(x==79990){
    //           data3=fopen(data_file3,"a");
    //           fprintf(data3,"%d\t%f\t%f\n",x,y,z);
    //           fclose(data3);
    //       }
    //   }
    //   fclose(gp);

	gp=popen("gnuplot -persist","w");
    fprintf(gp,"set terminal png\n");
    fprintf(gp,"set output 'Fig_env_kirkpatric_cost_self_%f.png'\n",K);
    fprintf(gp,"set xrange [0:%f]\n",1.0);
    fprintf(gp,"set xlabel 'T2'\n");
    fprintf(gp,"set yrange [0:%f]\n",1.0);
    fprintf(gp,"set ylabel 'P2'\n");
    fprintf(gp,"set key top left\n");
    // fprintf(gp,"u=%f\n",u);
    // fprintf(gp,"a=%f\n",a);
    // fprintf(gp,"f(x)=u*(x*(a*(1-u)-1)+1)/((1-u)*(a-1))\n");
    //fprintf(gp,"plot \'%s\' using 1:2 with lines linetype 1 title \"a= %f  S \",\'%s\' using 1:3 with lines linetype 3 title \"I \",\'%s\' using 1:4 with lines linetype 4 title \"0\"\n",data_file3,A_list[h],data_file3,data_file3);
    // fprintf(gp,"plot \'%s\' using 1:2 with points pointtype 7 lc rgb 'blue' title \"point \"\n",data_file1);
    fprintf(gp,"plot \'%s\' using 2:3 with points pointtype 7 lc rgb 'blue' title \"survivalrateK=%f\",\'%s\' using 2:3 with points pointtype 7 lc rgb 'red' title \"finalarrival\",\'%s\' every ::302::302 using 2:3 with points pointtype 13 pointsize 4 lc rgb 'yellow' notitle,\'%s\' every ::312::312 using 2:3 with points pointtype 13 pointsize 4 lc rgb 'black' notitle,\'%s\' every ::351::351 using 2:3 with points pointtype 13 pointsize 4 lc rgb 'violet' notitle\n",data_file1,K,data_file3,data_file1,data_file1,data_file1);
    // fprintf(gp,"plot f(x) with lines linetype 1 lc rgb 'red' title \"equilibria line\"\n");
    pclose(gp);

    data_file1="env_female_cost_self_0.075000.dat";
    data_file3 = "finalenv_cost_self_0.075000.dat";
    // gp=fopen(data_file1,"r");
    // while(fscanf(gp,"%d %lf %lf",&x,&y,&z)==3){
    //       if(x==79990){
    //           data3=fopen(data_file3,"a");
    //           fprintf(data3,"%d\t%f\t%f\n",x,y,z);
    //           fclose(data3);
    //       }
    //   }
    //   fclose(gp);

	gp=popen("gnuplot -persist","w");
    fprintf(gp,"set terminal png\n");
    fprintf(gp,"set output 'Fig_env_kirkpatric_cost_self_%f.png'\n",0.075000);
    fprintf(gp,"set xrange [0:%f]\n",1.0);
    fprintf(gp,"set xlabel 'T2'\n");
    fprintf(gp,"set yrange [0:%f]\n",1.0);
    fprintf(gp,"set ylabel 'P2'\n");
    fprintf(gp,"set key top left\n");
    // fprintf(gp,"u=%f\n",u);
    // fprintf(gp,"a=%f\n",a);
    // fprintf(gp,"f(x)=u*(x*(a*(1-u)-1)+1)/((1-u)*(a-1))\n");
    //fprintf(gp,"plot \'%s\' using 1:2 with lines linetype 1 title \"a= %f  S \",\'%s\' using 1:3 with lines linetype 3 title \"I \",\'%s\' using 1:4 with lines linetype 4 title \"0\"\n",data_file3,A_list[h],data_file3,data_file3);
    // fprintf(gp,"plot \'%s\' using 1:2 with points pointtype 7 lc rgb 'blue' title \"point \"\n",data_file1);
    fprintf(gp,"plot \'%s\' using 2:3 with points pointtype 7 lc rgb 'blue' title \"survivalrateK=%f\",\'%s\' using 2:3 with points pointtype 7 lc rgb 'red' title \"finalarrival\",\'%s\' every ::8051::8051 using 2:3 with points pointtype 13 pointsize 4 lc rgb 'yellow' notitle,\'%s\' every ::8301::8301 using 2:3 with points pointtype 13 pointsize 4 lc rgb 'black' notitle,\'%s\' every ::8561::8561 using 2:3 with points pointtype 13 pointsize 4 lc rgb 'violet' notitle\n",data_file1,K,data_file3,data_file1,data_file1,data_file1);
    // fprintf(gp,"plot f(x) with lines linetype 1 lc rgb 'red' title \"equilibria line\"\n");
    pclose(gp);


    
	
}
