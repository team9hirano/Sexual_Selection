/* nearest neighbor interaction */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MT.h"

#define LH 100	// 10000
#define LV 1000 // 10000
#define K 0.08
#define u 0.3
#define a 3.0
#define tend 300 // 80000
#define mapinitP 0.5
#define initialP 5
#define initialT 5

void Map(const char *sex, const char *filename, double initP) // int t, int h, int v
{
	FILE *gp;
	gp = popen("gnuplot -persist", "w");
	fprintf(gp, "set term pngcairo size 1000,1000\n");
	//  fprintf(gp,"set terminal png\n");
	if (strcmp(sex, "male") == 0)
	{
		if (strstr(filename, "stmap"))
			fprintf(gp, "set output 'MapSt_%s_env_cost_self_%g_initP_%g.png'\n", sex, K, initP);
		else if (strstr(filename, "intmap"))
			fprintf(gp, "set output 'OneInt_T1P1vsT2P2_%s_env_cost_self_%g_initP_%g.png'\n", sex, K, initP);
		else if (strstr(filename, "finmap"))
			fprintf(gp, "set output 'MapFin_%s_env_cost_self_%g_initP_%g.png'\n", sex, K, initP);
		fprintf(gp, "unset key\n");
		fprintf(gp, "set size ratio -1\n");
		fprintf(gp, "set xrange [%d:%d]\n", 0, LH);
		// fprintf(gp,"set xlabel 'T2'\n");
		fprintf(gp, "set yrange [%d:%d]\n", 1, tend);
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
			fprintf(gp, "set output 'MapSt_%s_env_cost_self_%g_initP_%g.png'\n", sex, K, initP);
		else if (strstr(filename, "intmap"))
			fprintf(gp, "set output 'OneInt_T1P1vsT2P2_%s_env_cost_self_%g_initP_%g.png'\n", sex, K, initP);
		else if (strstr(filename, "finmap"))
			fprintf(gp, "set output 'MapFin_%s_env_cost_self_%g_initP_%g.png'\n", sex, K, initP);
		fprintf(gp, "unset key\n");
		fprintf(gp, "set size ratio -1\n");
		fprintf(gp, "set xrange [%d:%d]\n", 0, LH);
		// fprintf(gp,"set xlabel 'T2'\n");
		fprintf(gp, "set yrange [%d:%d]\n", 1, tend);
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
	int *maleT, *base_maleT;
	int *maleP, *base_maleP;
	int *femaleT, *base_femaleT;
	int *femaleP, *base_femaleP;
	int *maleTdummy, *base_maleTdummy;
	int *malePdummy, *base_malePdummy;
	int *femaleTdummy, *base_femaleTdummy;
	int *femalePdummy, *base_femalePdummy;
	double initT2, initP2;
	int k, k2, i, j, i2, j2, t, ok, x1, x2;
	int maleI, maleJ, femaleI, femaleJ;
	int numMT1, numMT2, numMP1, numMP2;
	int numFT1, numFT2, numFP1, numFP2;
	double rnd, rnd2, sum1, sum2, sum3, sum4, sum5, y1, z1, y2, z2;
	int maleT0, maleP0, femaleT0, femaleP0, mgenotype, fgenotype;
	FILE *gp, *data1, *data2, *data3, *snapshot1, *snapshot2, *snapshot3, *snapshot4, *snapshot5, *snapshot6;
	char *data_file1, *data_file2, *data_file3, *snapshot_file1, *snapshot_file2, *snapshot_file3, *snapshot_file4, *snapshot_file5, *snapshot_file6;

	maleT = malloc(sizeof(int) * LH);
	maleP = malloc(sizeof(int) * LH);
	femaleT = malloc(sizeof(int) * LH);
	femaleP = malloc(sizeof(int) * LH);
	maleTdummy = malloc(sizeof(int) * LH);
	malePdummy = malloc(sizeof(int) * LH);
	femaleTdummy = malloc(sizeof(int) * LH);
	femalePdummy = malloc(sizeof(int) * LH);
	// base_maleT = malloc(sizeof(int) * LH * LV);
	// base_maleP = malloc(sizeof(int) * LH * LV);
	// base_femaleT = malloc(sizeof(int) * LH * LV);
	// base_femaleP = malloc(sizeof(int) * LH * LV);
	// base_maleTdummy = malloc(sizeof(int) * LH * LV);
	// base_malePdummy = malloc(sizeof(int) * LH * LV);
	// base_femaleTdummy = malloc(sizeof(int) * LH * LV);
	// base_femalePdummy = malloc(sizeof(int) * LH * LV);
	// for (i=0;i<LH;i++) {
	// 	maleT[i] = base_maleT + i * LV;
	// 	maleP[i] = base_maleP + i * LV;
	// 	femaleT[i] = base_femaleT + i * LV;
	// 	femaleP[i] = base_femaleP + i * LV;
	// 	maleTdummy[i] = base_maleTdummy + i * LV;
	// 	malePdummy[i] = base_malePdummy + i * LV;
	// 	femaleTdummy[i] = base_femaleTdummy + i * LV;
	// 	femalePdummy[i] = base_femalePdummy + i * LV;
	// }

	init_genrand(0);
	data_file1 = malloc(100);
	sprintf(data_file1, "Onedimention_%f.dat", K);
	snapshot_file2 = malloc(100);
	data_file2 = "flow.dat";

	for (k = initialT; k <= initialT; k++)
	{
		initT2 = 0.1 * k;
		for (k2 = initialP; k2 <= initialP; k2++)
		{
			initP2 = 0.1 * k2;

			// 2遺伝子型で対決T2P1vsT2P2 or T2P1vsT1P1
			for (i = 0; i < LH; i++)
			{
				// 一般パターン(4タイプ)
				// rnd = genrand_real2();
				// if (rnd < initT2)
				// 	maleT[i] = 2;
				// else
				// 	maleT[i] = 1;
				// rnd = genrand_real2();
				// if (rnd < initP2)
				// 	maleP[i] = 2;
				// else
				// 	maleP[i] = 1;
				// rnd = genrand_real2();
				// if (rnd < initT2)
				// 	femaleT[i] = 2;
				// else
				// 	femaleT[i] = 1;
				// rnd = genrand_real2();
				// if (rnd < initP2)
				// 	femaleP[i] = 2;
				// else
				// 	femaleP[i] = 1;

				// 2つの遺伝子型での対決
				// T2P1vsT2P2
				// maleT[i] = 2;
				// rnd = genrand_real2();
				// if (rnd < initP2)
				// 	maleP[i] = 2;
				// else
				//  	maleP[i] = 1;
				// femaleT[i] = 2;
				// rnd = genrand_real2();
				// if (rnd < initP2)
				// 	femaleP[i] = 2;
				// else
				// 	femaleP[i] = 1;

				// T1P1vsT2P1
				//  rnd = genrand_real2();
				//  if (rnd < initT2)
				//  	maleT[i] = 2;
				//  else
				//  	maleT[i] = 1;
				//  maleP[i] = 2;
				//  rnd = genrand_real2();
				//  if (rnd < initT2)
				//  	femaleT[i] = 2;
				//  else
				//  	femaleT[i] = 1;
				//  femaleP[i] = 2;

				// T1P1vaT2P2
				rnd = genrand_real2();
				if (rnd < initT2)
				{
					maleT[i] = 1;
					maleP[i] = 1;
				}
				else
				{
					maleT[i] = 2;
					maleP[i] = 2;
				}
				rnd = genrand_real2();
				if (rnd < initP2)
				{
					femaleT[i] = 1;
					femaleP[i] = 1;
				}
				else
				{
					femaleT[i] = 2;
					femaleP[i] = 2;
				}
			}

			numMT1 = numMT2 = numFT1 = numFT2 = numMP1 = numMP2 = numFP1 = numFP2 = 0;
			for (i = 0; i < LH; i++)
			{

				if (maleT[i] == 1)
					numMT1++;
				else if (maleT[i] == 2)
					numMT2++;
				if (femaleT[i] == 1)
					numFT1++;
				else if (femaleT[i] == 2)
					numFT2++;
				if (maleP[i] == 1)
					numMP1++;
				else if (maleP[i] == 2)
					numMP2++;
				if (femaleP[i] == 1)
					numFP1++;
				else if (femaleP[i] == 2)
					numFP2++;
			}
			printf("0 %f %f \n", (double)numMT2 / (double)(LH), (double)numMP2 / (double)(LH));
			//			printf("0 %f %f %f %f %f %f %f %f \n", (double)numMT1/(double)(LH*LV), (double)numMT2/(double)(LH*LV), (double)numFT1/(double)(LH*LV), (double)numFT2/(double)(LH*LV), (double)numMP1/(double)(LH*LV), (double)numMP2/(double)(LH*LV), (double)numFP1/(double)(LH*LV), (double)numFP2/(double)(LH*LV));

			t = 0;
			while (t < tend)
			{

				data1 = fopen(data_file1, "a");
				if (t % 10 == 0)
					fprintf(data1, "%d\t%f\t%f\n", t, (double)numMT2 / (double)(LH), (double)numMP2 / (double)(LH));
				fclose(data1);
				for (j = 0; j < LH; j++)
				{
					// i = genrand_int32() % LH;
					i = j;

					if (femaleP[i] == 1)
					{
						sum1 = sum2 = sum3 = sum4 = 0.0;
						i2 = (i - 1 + LH) % LH;
						if (maleT[i2] == 1 && maleP[i2] == 1)
							sum1 += 1.0;
						else if (maleT[i2] == 1 && maleP[i2] == 2)
							sum2 += 1.0;
						else if (maleT[i2] == 2 && maleP[i2] == 1)
							sum3 += 1.0;
						else if (maleT[i2] == 2 && maleP[i2] == 2)
							sum4 += 1.0;
						i2 = (i + 1 + LH) % LH;
						if (maleT[i2] == 1 && maleP[i2] == 1)
							sum1 += 1.0;
						else if (maleT[i2] == 1 && maleP[i2] == 2)
							sum2 += 1.0;
						else if (maleT[i2] == 2 && maleP[i2] == 1)
							sum3 += 1.0;
						else if (maleT[i2] == 2 && maleP[i2] == 2)
							sum4 += 1.0;
						i2 = i;
						if (maleT[i2] == 1 && maleP[i2] == 1)
							sum1 += 1.0;
						else if (maleT[i2] == 1 && maleP[i2] == 2)
							sum2 += 1.0;
						else if (maleT[i2] == 2 && maleP[i2] == 1)
							sum3 += 1.0;
						else if (maleT[i2] == 2 && maleP[i2] == 2)
							sum4 += 1.0;
					}
					else
					{
						sum1 = sum2 = sum3 = sum4 = 0.0;
						i2 = (i - 1 + LH) % LH;
						if (maleT[i2] == 1 && maleP[i2] == 1)
							sum1 += 1.0;
						else if (maleT[i2] == 1 && maleP[i2] == 2)
							sum2 += 1.0;
						else if (maleT[i2] == 2 && maleP[i2] == 1)
							sum3 += a;
						else if (maleT[i2] == 2 && maleP[i2] == 2)
							sum4 += a;
						i2 = (i + 1 + LH) % LH;
						if (maleT[i2] == 1 && maleP[i2] == 1)
							sum1 += 1.0;
						else if (maleT[i2] == 1 && maleP[i2] == 2)
							sum2 += 1.0;
						else if (maleT[i2] == 2 && maleP[i2] == 1)
							sum3 += a;
						else if (maleT[i2] == 2 && maleP[i2] == 2)
							sum4 += a;
						i2 = i;
						if (maleT[i2] == 1 && maleP[i2] == 1)
							sum1 += 1.0;
						else if (maleT[i2] == 1 && maleP[i2] == 2)
							sum2 += 1.0;
						else if (maleT[i2] == 2 && maleP[i2] == 1)
							sum3 += a;
						else if (maleT[i2] == 2 && maleP[i2] == 2)
							sum4 += a;
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
					if (femaleT[i2] == 1 && femaleP[i2] == 1)
						sum1 += 1.0;
					else if (femaleT[i2] == 1 && femaleP[i2] == 2)
						sum2 += 1.0;
					else if (femaleT[i2] == 2 && femaleP[i2] == 1)
						sum3 += 1.0;
					else if (femaleT[i2] == 2 && femaleP[i2] == 2)
						sum4 += 1.0;
					i2 = (i + 1 + LH) % LH;
					if (femaleT[i2] == 1 && femaleP[i2] == 1)
						sum1 += 1.0;
					else if (femaleT[i2] == 1 && femaleP[i2] == 2)
						sum2 += 1.0;
					else if (femaleT[i2] == 2 && femaleP[i2] == 1)
						sum3 += 1.0;
					else if (femaleT[i2] == 2 && femaleP[i2] == 2)
						sum4 += 1.0;
					i2 = i;
					if (femaleT[i2] == 1 && femaleP[i2] == 1)
						sum1 += 1.0;
					else if (femaleT[i2] == 1 && femaleP[i2] == 2)
						sum2 += 1.0;
					else if (femaleT[i2] == 2 && femaleP[i2] == 1)
						sum3 += 1.0;
					else if (femaleT[i2] == 2 && femaleP[i2] == 2)
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
						maleTdummy[i] = maleT0;
						malePdummy[i] = maleP0;
					}
					else
					{
						maleTdummy[i] = femaleT0;
						malePdummy[i] = femaleP0;
					}
					rnd = genrand_real2();
					if (rnd < 0.5)
					{
						femaleTdummy[i] = maleT0;
						femalePdummy[i] = maleP0;
					}
					else
					{
						femaleTdummy[i] = femaleT0;
						femalePdummy[i] = femaleP0;
					}
				}

				for (i = 0; i < LH; i++)
				{

					maleT[i] = maleTdummy[i];
					maleP[i] = malePdummy[i];
					femaleT[i] = femaleTdummy[i];
					femaleP[i] = femalePdummy[i];
				}

				numMT1 = numMT2 = numFT1 = numFT2 = numMP1 = numMP2 = numFP1 = numFP2 = 0;
				for (i = 0; i < LH; i++)
				{

					if (maleT[i] == 1)
						numMT1++;
					else if (maleT[i] == 2)
						numMT2++;
					if (femaleT[i] == 1)
						numFT1++;
					else if (femaleT[i] == 2)
						numFT2++;
					if (maleP[i] == 1)
						numMP1++;
					else if (maleP[i] == 2)
						numMP2++;
					if (femaleP[i] == 1)
						numFP1++;
					else if (femaleP[i] == 2)
						numFP2++;
				}
				t++;

				// 途中の図
				if (t < tend && fabs(initP2 - mapinitP) < 1e-12)
				{ // t<tend
					sprintf(snapshot_file2, "Onedi_intmap_cost_%f_initP_%g.dat", K, mapinitP);
					snapshot2 = fopen(snapshot_file2, "a");
					mgenotype = fgenotype = 0;
					for (i = 0; i < LH; i++)
					{

						if (maleT[i] == 1 && maleP[i] == 1)
							mgenotype = 1;
						else if (maleT[i] == 1 && maleP[i] == 2)
							mgenotype = 2;
						else if (maleT[i] == 2 && maleP[i] == 1)
							mgenotype = 3;
						else if (maleT[i] == 2 && maleP[i] == 2)
							mgenotype = 4;
						if (femaleT[i] == 1 && femaleP[i] == 1)
							fgenotype = 1;
						else if (femaleT[i] == 1 && femaleP[i] == 2)
							fgenotype = 2;
						else if (femaleT[i] == 2 && femaleP[i] == 1)
							fgenotype = 3;
						else if (femaleT[i] == 2 && femaleP[i] == 2)
							fgenotype = 4;
						fprintf(snapshot2, "%d\t%d\t%d\t%d\n", i, t, mgenotype, fgenotype);
					}
					fclose(snapshot2);
					Map("male", snapshot_file2, mapinitP);
					Map("female", snapshot_file2, mapinitP);
				}
			}
		}
	}
	data_file3 = malloc(100);
	sprintf(data_file3, "Oned_final_%f.dat", K);
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
	// 矢印
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
	fprintf(gp, "set output 'OnedimentionT2P2_%f.png'\n", K);
	fprintf(gp, "set xrange [0:%f]\n", 1.0);
	fprintf(gp, "set xlabel 'T2'\n");
	fprintf(gp, "set yrange [0:%f]\n", 1.0);
	fprintf(gp, "set ylabel 'P2'\n");
	fprintf(gp, "u=%f\n", u);
	fprintf(gp, "a=%f\n", a);
	fprintf(gp, "f(x)=u*(x*(a*(1-u)-1)+1)/((1-u)*(a-1))\n");
	// fprintf(gp,"plot \'%s\' using 1:2 with lines linetype 1 title \"a= %f  S \",\'%s\' using 1:3 with lines linetype 3 title \"I \",\'%s\' using 1:4 with lines linetype 4 title \"0\"\n",data_file3,A_list[h],data_file3,data_file3);
	//  fprintf(gp,"plot \'%s\' using 1:2 with points pointtype 7 lc rgb 'blue' title \"point \"\n",data_file1);
	fprintf(gp, "plot \'%s\' using 2:3 with points pointtype 7 lc rgb 'blue' title \"survivalrateK=%f\",\'%s\' using 1:2:3:4 with vectors head filled lc rgb 'blue',\'%s\' using 2:3 with points pointtype 7 lc rgb 'red' title \"finalarrival\"\n", data_file1, K, data_file2, data_file3);
	// fprintf(gp,"plot f(x) with lines linetype 1 lc rgb 'red' title \"equilibria line\"\n");
	pclose(gp);

	free(femalePdummy);
	free(femaleTdummy);
	free(malePdummy);
	free(maleTdummy);
	free(femaleP);
	free(femaleT);
	free(maleP);
	free(maleT);
	free(data_file1);
	free(data_file3);
	free(snapshot_file2);
}
