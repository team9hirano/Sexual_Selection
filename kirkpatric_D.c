/* if males' traits are 2, then these males cannot choose at probability u for each female */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "MT.h"

#define LH	1000	//10000
#define LV	1000	//10000

#define u	0.3
// #define K  0.1
#define a	3.0
#define tend	100//20000 	//100

// int main(int argc,char *argv[])
int main()
{
	int **maleT, *base_maleT;
	int **maleP, *base_maleP;
	int **femaleT, *base_femaleT;
	int **femaleP, *base_femaleP;
    int **D,*base_D;
	int **maleTdummy, *base_maleTdummy;
	int **malePdummy, *base_malePdummy;
	int **femaleTdummy, *base_femaleTdummy;
	int **femalePdummy, *base_femalePdummy;
    int **Ddummy,*base_Ddummy;
    double **dia,*base_dia;
	double initT2, initP2;
    // double K=atof(argv[1]);
	int k, k2, i, j, i2, j2, t, ok,x;
	int maleI, maleJ,femaleI,femaleJ;
	int numMT1, numMT2, numMP1, numMP2;
	int numFT1, numFT2, numFP1, numFP2;
	double rnd, rnd2, sum1, sum2, sum3, sum4,y,z,w;
	int maleT0, maleP0,femaleT0,femaleP0;
    double D_con;
    FILE *gp,*data1,*data2,*data3;
    char *data_file1,*data_file2,*data_file3;

    
	maleT = malloc(sizeof(int *) * LH);
	maleP = malloc(sizeof(int *) * LH);
	femaleT = malloc(sizeof(int *) * LH);
	femaleP = malloc(sizeof(int *) * LH);
    D=malloc(sizeof(int *) * LH);
	maleTdummy = malloc(sizeof(int *) * LH);
	malePdummy = malloc(sizeof(int *) * LH);
	femaleTdummy = malloc(sizeof(int *) * LH);
	femalePdummy = malloc(sizeof(int *) * LH);
    Ddummy=malloc(sizeof(int *) * LH);
	
    base_maleT = malloc(sizeof(int) * LH * LV);
	base_maleP = malloc(sizeof(int) * LH * LV);
	base_femaleT = malloc(sizeof(int) * LH * LV);
	base_femaleP = malloc(sizeof(int) * LH * LV);
    base_D = malloc(sizeof(int) * LH * LV);
	base_maleTdummy = malloc(sizeof(int) * LH * LV);
	base_malePdummy = malloc(sizeof(int) * LH * LV);
	base_femaleTdummy = malloc(sizeof(int) * LH * LV);
	base_femalePdummy = malloc(sizeof(int) * LH * LV);
    base_Ddummy = malloc(sizeof(int) * LH * LV);
	

    for (i=0;i<LH;i++) {
		maleT[i] = base_maleT + i * LV;
		maleP[i] = base_maleP + i * LV;
		femaleT[i] = base_femaleT + i * LV;
		femaleP[i] = base_femaleP + i * LV;
        D[i] = base_femaleP + i * LV;
		maleTdummy[i] = base_maleTdummy + i * LV;
		malePdummy[i] = base_malePdummy + i * LV;
		femaleTdummy[i] = base_femaleTdummy + i * LV;
		femalePdummy[i] = base_femalePdummy + i * LV;
        Ddummy[i] = base_femaleP + i * LV;
	}

	init_genrand(0);
    data_file1="kirkpatricD.dat";
	// data_file1 = malloc(100);
    // sprintf(data_file1, "kirkpatricD.dat");
    
	for (k=1; k<=1; k++)
	{
		initT2 = 0.1*k;
		for (k2=1; k2<=9; k2++)
		{
			initP2 = 0.1*k2;

			for (i=0; i<LH; i++) {
				for (j=0; j<LV; j++) {
					rnd = genrand_real2();
					if(rnd < initT2) maleT[i][j] = 2;
					else maleT[i][j] = 1;
					rnd = genrand_real2();
					if(rnd < initP2) maleP[i][j] = 2;
					else maleP[i][j] = 1;
					rnd = genrand_real2();
					if(rnd < initT2) femaleT[i][j] = 2;
					else femaleT[i][j] = 1;
					rnd = genrand_real2();
					if(rnd < initP2) femaleP[i][j] = 2;
					else femaleP[i][j] = 1;
				}
			}

			numMT1 = numMT2 = numFT1 = numFT2 = numMP1 = numMP2 = numFP1 = numFP2 = 0;
			for (i=0; i<LH; i++) {
				for (j=0; j<LV; j++) {
					if (maleT[i][j] == 1) numMT1++;
					else if (maleT[i][j] == 2) numMT2++;
					if (femaleT[i][j] == 1) numFT1++;
					else if (femaleT[i][j] == 2) numFT2++;
					if (maleP[i][j] == 1) numMP1++;
					else if (maleP[i][j] == 2) numMP2++;
					if (femaleP[i][j] == 1) numFP1++;
					else if (femaleP[i][j] == 2) numFP2++;
				}
			}
			printf("0 %f %f \n", (double)numMT2/(double)(LH*LV), (double)numMP2/(double)(LH*LV));
//			printf("0 %f %f %f %f %f %f %f %f \n", (double)numMT1/(double)(LH*LV), (double)numMT2/(double)(LH*LV), (double)numFT1/(double)(LH*LV), (double)numFT2/(double)(LH*LV), (double)numMP1/(double)(LH*LV), (double)numMP2/(double)(LH*LV), (double)numFP1/(double)(LH*LV), (double)numFP2/(double)(LH*LV));

			t = 0;
			while (t<tend)
			{   //D is inserted
                D_con=0.0;
                sum1 = sum2 = sum3 = sum4 = 0.0;
				for (i=0; i<LH; i++)
					for (j=0; j<LV; j++)
					{
						if (maleT[i][j]==1 && maleP[i][j]==1) sum1 += 1.0;
						else if (maleT[i][j]==1 && maleP[i][j]==2) sum2 += 1.0;
						else if (maleT[i][j]==2 && maleP[i][j]==1) sum3 += 1.0;
						else if (maleT[i][j]==2 && maleP[i][j]==2) sum4 += 1.0;
					}
                D_con=(sum1/(sum1+sum2+sum3+sum4))*(sum4/(sum1+sum2+sum3+sum4))-(sum2/(sum1+sum2+sum3+sum4))*(sum3/(sum1+sum2+sum3+sum4));

				data1=fopen(data_file1,"a");
				// if(t%10==0)
                fprintf(data1,"%d\t%f\t%f\t%f\n",t,(double)numMT2/(double)(LH*LV), (double)numMP2/(double)(LH*LV),D_con);
                fclose(data1);
                sum1 = sum2 = sum3 = sum4 = 0.0;
				for (i=0; i<LH; i++)
					for (j=0; j<LV; j++)
					{
						if (maleT[i][j]==1 && maleP[i][j]==1) sum1 += 1.0;
						else if (maleT[i][j]==1 && maleP[i][j]==2) sum2 += 1.0;
						else if (maleT[i][j]==2 && maleP[i][j]==1) sum3 += a;
						else if (maleT[i][j]==2 && maleP[i][j]==2) sum4 += a;
					}
				for (i=0; i<LH; i++)
					for (j=0; j<LV; j++) {
						if (femaleP[i][j]==1)
						{
							do
							{
								maleI = (int)(LH*genrand_real2());
								maleJ = (int)(LH*genrand_real2());
								if (maleT[maleI][maleJ]==1) rnd = 1.0;
								else rnd = genrand_real2();
							} while (rnd<u);

                            maleT0 = maleT[maleI][maleJ];
							maleP0 = maleP[maleI][maleJ];
                            //メスのコスト
                            // do
							// {
							// 	femaleI = (int)(LV*genrand_real2());
							// 	femaleJ = (int)(LV*genrand_real2());
							// 	if (femaleP[femaleI][femaleJ]==1) rnd = 1.0;
							// 	else rnd = genrand_real2();
							// } while (rnd<K);

                            // femaleT0 = femaleT[femaleI][femaleJ];
							// femaleP0 = femaleP[femaleI][femaleJ];
						}
						else
						{
							do
							{
								rnd = genrand_real2();
								if (rnd < sum1/(sum1+sum2+sum3+sum4))
								{
									maleT0 = 1;
									maleP0 = 1;
								}
								else if (rnd < (sum1+sum2)/(sum1+sum2+sum3+sum4))
								{
									maleT0 = 1;
									maleP0 = 2;
								}
								else if (rnd < (sum1+sum2+sum3)/(sum1+sum2+sum3+sum4))
								{
									maleT0 = 2;
									maleP0 = 1;
								}
								else
								{
									maleT0 = 2;
									maleP0 = 2;
								}
								if (maleT0==1) rnd2 = 1.0;
								else rnd2 = genrand_real2();
							} while (rnd2<u);
                            //メスのコスト
                            // do
							// {
							// 	femaleI = (int)(LV*genrand_real2());
							// 	femaleJ = (int)(LV*genrand_real2());
							// 	if (femaleP[femaleI][femaleJ]==1) rnd = 1.0;
							// 	else rnd = genrand_real2();
							// } while (rnd<K);

                            // femaleT0 = femaleT[femaleI][femaleJ];
							// femaleP0 = femaleP[femaleI][femaleJ];

						}
						rnd = genrand_real2();
						if (rnd < 0.5)
						{
							maleTdummy[i][j] = maleT0;
							malePdummy[i][j] = maleP0;
						}
						else
						{
							maleTdummy[i][j] = femaleT[i][j];//femaleT0
							malePdummy[i][j] = femaleP[i][j];//femaleP0
						}
						rnd = genrand_real2();
						if (rnd < 0.5)
						{
							femaleTdummy[i][j] = maleT0;
							femalePdummy[i][j] = maleP0;
						}
						else
						{
							femaleTdummy[i][j] = femaleT[i][j];//femaleT0
							femalePdummy[i][j] = femaleP[i][j];//femaleP0
						}
					}
				for (i=0; i<LH; i++)
					for (j=0; j<LV; j++) {
						maleT[i][j] = maleTdummy[i][j];
						maleP[i][j] = malePdummy[i][j];
						femaleT[i][j] = femaleTdummy[i][j];
						femaleP[i][j] = femalePdummy[i][j];
					}
				numMT1 = numMT2 = numFT1 = numFT2 = numMP1 = numMP2 = numFP1 = numFP2 = 0;
				for (i=0; i<LH; i++) {
					for (j=0; j<LV; j++) {
						if (maleT[i][j] == 1) numMT1++;
						else if (maleT[i][j] == 2) numMT2++;
						if (femaleT[i][j] == 1) numFT1++;
						else if (femaleT[i][j] == 2) numFT2++;
						if (maleP[i][j] == 1) numMP1++;
						else if (maleP[i][j] == 2) numMP2++;
						if (femaleP[i][j] == 1) numFP1++;
						else if (femaleP[i][j] == 2) numFP2++;
					}
				}
				t++;
				// printf("%d %f %f \n", t, (double)numMT2/(double)(LH*LV), (double)numMP2/(double)(LH*LV));
//				printf("%d %f %f %f %f %f %f %f %f \n", t, (double)numMT1/(double)(LH*LV), (double)numMT2/(double)(LH*LV), (double)numFT1/(double)(LH*LV), (double)numFT2/(double)(LH*LV), (double)numMP1/(double)(LH*LV), (double)numMP2/(double)(LH*LV), (double)numFP1/(double)(LH*LV), (double)numFP2/(double)(LH*LV));
			}
		}
	}
    // data_file2="female_cost.dat";
    data_file3 = "final_kirkpatricD.dat";
	// data_file3=malloc(100);
    // sprintf(data_file3, "finalnoenv_cost_%f.dat",K);
	gp=fopen(data_file1,"r");
    while(fscanf(gp,"%d %lf %lf %lf",&x,&y,&z,&w)==4){
         if(x==(tend-10)){
             data3=fopen(data_file3,"a");
             fprintf(data3,"%d\t%f\t%f\t%f\n",x,y,z,w);
             fclose(data3);
         }
     }
     fclose(gp);
	
	gp=popen("gnuplot -persist","w");
    fprintf(gp,"set terminal png\n");
    fprintf(gp,"set output 'Fig_kirkpatricD.png'\n");
    fprintf(gp,"set xrange [0:%f]\n",1.0);
    fprintf(gp,"set xlabel 'T2'\n");
    fprintf(gp,"set yrange [0:%f]\n",1.0);
    fprintf(gp,"set ylabel 'P2'\n");
    fprintf(gp,"set zlabel 'D'\n");
    // fprintf(gp,"u=%f\n",u);
    // fprintf(gp,"a=%f\n",a);
    // fprintf(gp,"f(x)=u*(x*(a*(1-u)-1)+1)/((1-u)*(a-1))\n");
    //fprintf(gp,"plot \'%s\' using 1:2 with lines linetype 1 title \"a= %f  S \",\'%s\' using 1:3 with lines linetype 3 title \"I \",\'%s\' using 1:4 with lines linetype 4 title \"0\"\n",data_file3,A_list[h],data_file3,data_file3);
    // fprintf(gp,"plot \'%s\' using 1:2 with points pointtype 7 lc rgb 'blue' title \"point \"\n",data_file1);
    // fprintf(gp,"plot \'%s\' using 2:3 with points pointtype 7 lc rgb 'blue' title \"survivalrateK=%f\",\'%s\' using 2:3 with points pointtype 7 lc rgb 'red' title \"finalarrival\",f(x) with lines linetype 1 lc rgb 'red' title \"equilibria line\"\n",data_file1,K,data_file3);
    // fprintf(gp,"plot f(x) with lines linetype 1 lc rgb 'red' title \"equilibria line\"\n");
    fprintf(gp,"splot \'%s\' using 2:3:4 with points pointtype 7 lc rgb 'blue' title \"frequency evolution\",\'%s\' using 2:3:4 with lines linetype 1 lc rgb 'red' title \"finalarrival\"\n",data_file1,data_file3);
    pclose(gp);
    
	free(base_femalePdummy);
	free(base_femaleTdummy);
	free(base_malePdummy);
	free(base_maleTdummy);
    free(base_Ddummy);
	free(base_femaleP);
	free(base_femaleT);
	free(base_maleP);
	free(base_maleT);
    free(base_D);
	free(femalePdummy);
	free(femaleTdummy);
	free(malePdummy);
	free(maleTdummy);
    free(Ddummy);
	free(femaleP);
	free(femaleT);
	free(maleP);
	free(maleT);
    free(D);
	// free(data_file1);
	// free(data_file3);
}
