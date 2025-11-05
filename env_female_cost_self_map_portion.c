/* nearest neighbor interaction */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MT.h"

#define LH	1000	//10000
#define LV	1000	//10000
#define K  0.075
#define u	0.3
#define a	3.0
#define tend	3000//4000 80000 10000
#define mapinitP 0.2
#define initialP 2


void Map(const char *sex,const char *filename,double initP,int t,int h,int v){
    FILE *gp;
    gp=popen("gnuplot -persist","w");
    fprintf(gp,"set term pngcairo size 1000,1000\n");
    //  fprintf(gp,"set terminal png\n");
    if(strcmp(sex, "male") == 0){
        if(strstr(filename, "stmap")) fprintf(gp,"set output 'MapSt_%s_env_cost_self_%g_initP_%g.png'\n",sex,K,initP);
        else if(strstr(filename, "intmap")) fprintf(gp,"set output 'MapInt_%s_env_cost_self_%g_initP_%g_t_%d_portion_%d.png'\n",sex,K,initP,t,h);
        else if(strstr(filename, "finmap")) fprintf(gp,"set output 'MapFin_%s_env_cost_self_%g_initP_%g.png'\n",sex,K,initP);
        fprintf(gp,"unset key\n");
        fprintf(gp,"set size ratio -1\n");
        fprintf(gp,"set xrange [%d:%d]\n",h,h+100);
        // fprintf(gp,"set xlabel 'T2'\n");
        fprintf(gp,"set yrange [%d:%d]\n",v,v+100);
        fprintf(gp,"set palette defined(1 'blue',2 'green',3 'orange',4 'red')\n");
        // fprintf(gp,"set multiplot layout 1,2 title 'Genotype map (T×P: 0=T1P1, 1=T1P2, 2=T2P1, 3=T2P2)'\n");
        fprintf(gp,"set title 'Male map'\n");
        fprintf(gp,"unset xtics;unset ytics\n");
        // fprintf(gp,"unset yticks\n");
        fprintf(gp,"plot \'%s\' using 1:2:3 with image\n",filename);
    }
    else if(strcmp(sex, "female") == 0){
        if(strstr(filename, "stmap")) fprintf(gp,"set output 'MapSt_%s_env_cost_self_%g_initP_%g.png'\n",sex,K,initP);
        else if(strstr(filename, "intmap")) fprintf(gp,"set output 'MapInt_%s_env_cost_self_%g_initP_%g_t_%d_portion_%d.png'\n",sex,K,initP,t,h);
        else if(strstr(filename, "finmap")) fprintf(gp,"set output 'MapFin_%s_env_cost_self_%g_initP_%g.png'\n",sex,K,initP);
        fprintf(gp,"unset key\n");
        fprintf(gp,"set size ratio -1\n");
        fprintf(gp,"set xrange [%d:%d]\n",h,h+100);
        // fprintf(gp,"set xlabel 'T2'\n");
        fprintf(gp,"set yrange [%d:%d]\n",v,v+100);
        fprintf(gp,"set palette defined(1 'blue',2 'green',3 'orange',4 'red')\n");
        // fprintf(gp,"set multiplot layout 1,2 title 'Genotype map (T×P: 0=T1P1, 1=T1P2, 2=T2P1, 3=T2P2)'\n");
        fprintf(gp,"set title 'Female map'\n");
        fprintf(gp,"unset xtics;unset ytics\n");
        // fprintf(gp,"unset yticks\n");
        fprintf(gp,"plot \'%s\' using 1:2:4 with image\n",filename);
    }else printf("それはだめよ");
    
    pclose(gp);
}



int main(void)
{
	//新しく挿入
	// int tend;
	int **maleT, *base_maleT;
	int **maleP, *base_maleP;
	int **femaleT, *base_femaleT;
	int **femaleP, *base_femaleP;
	int **maleTdummy, *base_maleTdummy;
	int **malePdummy, *base_malePdummy;
	int **femaleTdummy, *base_femaleTdummy;
	int **femalePdummy, *base_femalePdummy;
	double initT2, initP2;
	int k, k2, i, j, i2, j2, t, ok,x;
	int maleI, maleJ,femaleI,femaleJ;
	int numMT1, numMT2, numMP1, numMP2;
	int numFT1, numFT2, numFP1, numFP2;
	double rnd, rnd2, sum1, sum2, sum3,sum4,sum5,y,z;
	int maleT0, maleP0,femaleT0,femaleP0,mgenotype,fgenotype;
    FILE *gp,*data1,*data2,*data3,*snapshot1,*snapshot2,*snapshot3,*snapshot4,*snapshot5,*snapshot6;
    char *data_file1,*data_file2,*data_file3,*snapshot_file1,*snapshot_file2,*snapshot_file3,*snapshot_file4,*snapshot_file5,*snapshot_file6;

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
	for (i=0;i<LH;i++) {
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
	

    snapshot_file1=malloc(100);
    sprintf(snapshot_file1, "efcs_stmap_%f_initP_%g.dat", K,mapinitP);

    snapshot_file2=malloc(100);
    // sprintf(snapshot2, "efcs_intmap_%f.dat", K);

    snapshot_file3=malloc(100);
    sprintf(snapshot_file3, "efcs_finmap_%f_initP_%g.dat", K,mapinitP);

    
    

	for (k=1; k<=1; k++)
	{
		initT2 = 0.1*k;
		for (k2=initialP; k2<=initialP; k2++)//k2=1; k2<=9; k2++
		{
			initP2 = 0.1*k2;
			//新しく挿入
			// if(initP2==0.1){
			// 	tend=80000;
			// }else{
			// 	tend=40000;
			// }
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
            //初期の図
            if(initP2==mapinitP){
                snapshot1=fopen(snapshot_file1,"w");
                mgenotype=fgenotype= 0;
                for (i=0; i<LH; i++) {
                    for (j=0; j<LV; j++) {
                        if (maleT[i][j] == 1 && maleP[i][j]==1 ) mgenotype=1;
                        else if (maleT[i][j] == 1 && maleP[i][j]==2 ) mgenotype=2;
                        else if (maleT[i][j] == 2 && maleP[i][j]==1 ) mgenotype=3;
                        else if (maleT[i][j] == 2 && maleP[i][j]==2 ) mgenotype=4;
                        if (femaleT[i][j] == 1 && femaleP[i][j]==1 ) fgenotype=1;
                        else if (femaleT[i][j] == 1 && femaleP[i][j]==2 ) fgenotype=2;
                        else if (femaleT[i][j] == 2 && femaleP[i][j]==1 ) fgenotype=3;
                        else if (femaleT[i][j] == 2 && femaleP[i][j]==2 ) fgenotype=4;
                        fprintf(snapshot1,"%d\t%d\t%d\t%d\n",i,j,mgenotype,fgenotype);
                        
                    }
                }
                fclose(snapshot1);
            }

			t = 0;
			while (t<tend)
			{
				
				for (i=0; i<LH; i++)
					for (j=0; j<LV; j++) {
						if (femaleP[i][j]==1)
						{
							sum1 = sum2 = sum3 = sum4 = 0.0;
							i2 = (i-1+LH)%LH; j2 = j;
							if (maleT[i2][j2]==1 && maleP[i2][j2]==1) sum1 += 1.0;
							else if (maleT[i2][j2]==1 && maleP[i2][j2]==2) sum2 += 1.0;
							else if (maleT[i2][j2]==2 && maleP[i2][j2]==1) sum3 += 1.0;
							else if (maleT[i2][j2]==2 && maleP[i2][j2]==2) sum4 += 1.0;
							i2 = (i+1+LH)%LH; j2 = j;
							if (maleT[i2][j2]==1 && maleP[i2][j2]==1) sum1 += 1.0;
							else if (maleT[i2][j2]==1 && maleP[i2][j2]==2) sum2 += 1.0;
							else if (maleT[i2][j2]==2 && maleP[i2][j2]==1) sum3 += 1.0;
							else if (maleT[i2][j2]==2 && maleP[i2][j2]==2) sum4 += 1.0;
							i2 = i; j2 = (j-1+LV)%LV;
							if (maleT[i2][j2]==1 && maleP[i2][j2]==1) sum1 += 1.0;
							else if (maleT[i2][j2]==1 && maleP[i2][j2]==2) sum2 += 1.0;
							else if (maleT[i2][j2]==2 && maleP[i2][j2]==1) sum3 += 1.0;
							else if (maleT[i2][j2]==2 && maleP[i2][j2]==2) sum4 += 1.0;
							i2 = i; j2 = (j+1+LV)%LV;
							if (maleT[i2][j2]==1 && maleP[i2][j2]==1) sum1 += 1.0;
							else if (maleT[i2][j2]==1 && maleP[i2][j2]==2) sum2 += 1.0;
							else if (maleT[i2][j2]==2 && maleP[i2][j2]==1) sum3 += 1.0;
							else if (maleT[i2][j2]==2 && maleP[i2][j2]==2) sum4 += 1.0;

                            //新たに挿入
                            i2 = i; j2 = j;
							if (maleT[i2][j2]==1 && maleP[i2][j2]==1) sum1 += 1.0;
							else if (maleT[i2][j2]==1 && maleP[i2][j2]==2) sum2 += 1.0;
							else if (maleT[i2][j2]==2 && maleP[i2][j2]==1) sum3 += 1.0;
							else if (maleT[i2][j2]==2 && maleP[i2][j2]==2) sum4 += 1.0;
						}
						else
						{
							sum1 = sum2 = sum3 = sum4 = 0.0;
							i2 = (i-1+LH)%LH; j2 = j;
							if (maleT[i2][j2]==1 && maleP[i2][j2]==1) sum1 += 1.0;
							else if (maleT[i2][j2]==1 && maleP[i2][j2]==2) sum2 += 1.0;
							else if (maleT[i2][j2]==2 && maleP[i2][j2]==1) sum3 += a;
							else if (maleT[i2][j2]==2 && maleP[i2][j2]==2) sum4 += a;
							i2 = (i+1+LH)%LH; j2 = j;
							if (maleT[i2][j2]==1 && maleP[i2][j2]==1) sum1 += 1.0;
							else if (maleT[i2][j2]==1 && maleP[i2][j2]==2) sum2 += 1.0;
							else if (maleT[i2][j2]==2 && maleP[i2][j2]==1) sum3 += a;
							else if (maleT[i2][j2]==2 && maleP[i2][j2]==2) sum4 += a;
							i2 = i; j2 = (j-1+LV)%LV;
							if (maleT[i2][j2]==1 && maleP[i2][j2]==1) sum1 += 1.0;
							else if (maleT[i2][j2]==1 && maleP[i2][j2]==2) sum2 += 1.0;
							else if (maleT[i2][j2]==2 && maleP[i2][j2]==1) sum3 += a;
							else if (maleT[i2][j2]==2 && maleP[i2][j2]==2) sum4 += a;
							i2 = i; j2 = (j+1+LV)%LV;
							if (maleT[i2][j2]==1 && maleP[i2][j2]==1) sum1 += 1.0;
							else if (maleT[i2][j2]==1 && maleP[i2][j2]==2) sum2 += 1.0;
							else if (maleT[i2][j2]==2 && maleP[i2][j2]==1) sum3 += a;
							else if (maleT[i2][j2]==2 && maleP[i2][j2]==2) sum4 += a;
                            //新たに挿入
                            i2 = i; j2 = j;
							if (maleT[i2][j2]==1 && maleP[i2][j2]==1) sum1 += 1.0;
							else if (maleT[i2][j2]==1 && maleP[i2][j2]==2) sum2 += 1.0;
							else if (maleT[i2][j2]==2 && maleP[i2][j2]==1) sum3 += a;
							else if (maleT[i2][j2]==2 && maleP[i2][j2]==2) sum4 += a;
						}

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
                        sum1 = sum2 = sum3 = sum4 = 0.0;
						i2 = (i-1+LH)%LH; j2 = j;
						if (femaleT[i2][j2]==1 && femaleP[i2][j2]==1) sum1 += 1.0;
						else if (femaleT[i2][j2]==1 && femaleP[i2][j2]==2) sum2 += 1.0;
						else if (femaleT[i2][j2]==2 && femaleP[i2][j2]==1) sum3 += 1.0;
						else if (femaleT[i2][j2]==2 && femaleP[i2][j2]==2) sum4 += 1.0;
						i2 = (i+1+LH)%LH; j2 = j;
						if (femaleT[i2][j2]==1 && femaleP[i2][j2]==1) sum1 += 1.0;
						else if (femaleT[i2][j2]==1 && femaleP[i2][j2]==2) sum2 += 1.0;
						else if (femaleT[i2][j2]==2 && femaleP[i2][j2]==1) sum3 += 1.0;
						else if (femaleT[i2][j2]==2 && femaleP[i2][j2]==2) sum4 += 1.0;
						i2 = i; j2 = (j-1+LV)%LV;
						if (femaleT[i2][j2]==1 && femaleP[i2][j2]==1) sum1 += 1.0;
						else if (femaleT[i2][j2]==1 && femaleP[i2][j2]==2) sum2 += 1.0;
						else if (femaleT[i2][j2]==2 && femaleP[i2][j2]==1) sum3 += 1.0;
						else if (femaleT[i2][j2]==2 && femaleP[i2][j2]==2) sum4 += 1.0;
						i2 = i; j2 = (j+1+LV)%LV;
						if (femaleT[i2][j2]==1 && femaleP[i2][j2]==1) sum1 += 1.0;
						else if (femaleT[i2][j2]==1 && femaleP[i2][j2]==2) sum2 += 1.0;
						else if (femaleT[i2][j2]==2 && femaleP[i2][j2]==1) sum3 += 1.0;
						else if (femaleT[i2][j2]==2 && femaleP[i2][j2]==2) sum4 += 1.0;
                        //新たに挿入
                            i2 = i; j2 = j;
							if (femaleT[i2][j2]==1 && femaleP[i2][j2]==1) sum1 += 1.0;
							else if (femaleT[i2][j2]==1 && femaleP[i2][j2]==2) sum2 += 1.0;
							else if (femaleT[i2][j2]==2 && femaleP[i2][j2]==1) sum3 += 1.0;
							else if (femaleT[i2][j2]==2 && femaleP[i2][j2]==2) sum4 += 1.0;    
						
						do
						{
							rnd = genrand_real2();
							if (rnd < sum1/(sum1+sum2+sum3+sum4))
							{
								femaleT0 = 1;
								femaleP0 = 1;
							}
							else if (rnd < (sum1+sum2)/(sum1+sum2+sum3+sum4))
							{
								femaleT0 = 1;
								femaleP0 = 2;
							}
							else if (rnd < (sum1+sum2+sum3)/(sum1+sum2+sum3+sum4))
							{
								femaleT0 = 2;
								femaleP0 = 1;
							}
							else
							{
								femaleT0 = 2;
								femaleP0 = 2;
							}
							if (femaleP0==1) rnd2 = 1.0;
							else rnd2 = genrand_real2();	
							} while (rnd2<K);

                            

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

                //途中の図
                if(t>tend-101 && t< tend && initP2==mapinitP){//t<tend,t%10==0
                    sprintf(snapshot_file2, "efcs_intmap_t_%d_cost_%f_initP_%g_portion.dat", t,K,mapinitP);
                    snapshot2=fopen(snapshot_file2,"w");
                    mgenotype=fgenotype= 0;
                    for (i=0; i<LH; i++) {
                        for (j=0; j<LV; j++) {
                            if (maleT[i][j] == 1 && maleP[i][j]==1 ) mgenotype=1;
                            else if (maleT[i][j] == 1 && maleP[i][j]==2 ) mgenotype=2;
                            else if (maleT[i][j] == 2 && maleP[i][j]==1 ) mgenotype=3;
                            else if (maleT[i][j] == 2 && maleP[i][j]==2 ) mgenotype=4;
                            if (femaleT[i][j] == 1 && femaleP[i][j]==1 ) fgenotype=1;
                            else if (femaleT[i][j] == 1 && femaleP[i][j]==2 ) fgenotype=2;
                            else if (femaleT[i][j] == 2 && femaleP[i][j]==1 ) fgenotype=3;
                            else if (femaleT[i][j] == 2 && femaleP[i][j]==2 ) fgenotype=4;
                            fprintf(snapshot2,"%d\t%d\t%d\t%d\n",i,j,mgenotype,fgenotype);
                        }
                    }
                    fclose(snapshot2);
					
					Map("male",snapshot_file2,mapinitP,t,100,0);
    				Map("female",snapshot_file2,mapinitP,t,100,0);
					
                }
				
			}
			//最後の図
			if(initP2==mapinitP){
				snapshot3=fopen(snapshot_file3,"w");
				mgenotype=fgenotype= 0;
				for (i=0; i<LH; i++) {
					for (j=0; j<LV; j++) {
						if (maleT[i][j] == 1 && maleP[i][j]==1 ) mgenotype=1;
						else if (maleT[i][j] == 1 && maleP[i][j]==2 ) mgenotype=2;
						else if (maleT[i][j] == 2 && maleP[i][j]==1 ) mgenotype=3;
						else if (maleT[i][j] == 2 && maleP[i][j]==2 ) mgenotype=4;
						if (femaleT[i][j] == 1 && femaleP[i][j]==1 ) fgenotype=1;
						else if (femaleT[i][j] == 1 && femaleP[i][j]==2 ) fgenotype=2;
						else if (femaleT[i][j] == 2 && femaleP[i][j]==1 ) fgenotype=3;
						else if (femaleT[i][j] == 2 && femaleP[i][j]==2 ) fgenotype=4;
						fprintf(snapshot3,"%d\t%d\t%d\t%d\n",i,j,mgenotype,fgenotype);
					}
				}
				fclose(snapshot3);
			}
		}
	}
    



	
	// Map("male",snapshot_file1,mapinitP,0,LH,LV);
	// Map("female",snapshot_file1,mapinitP,0,LH,LV);
	// Map("male",snapshot_file3,mapinitP,0,LH,LV);
	// Map("female",snapshot_file3,mapinitP,0,LH,LV);
	

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
	// free(data_file1);
    free(snapshot_file1);
    free(snapshot_file2);
    free(snapshot_file3);
	// free(data_file3);

}
