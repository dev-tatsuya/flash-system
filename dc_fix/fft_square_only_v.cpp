//直流
//組織変化無し

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
// #include "wingxa.h"

using namespace std;

#define DRND(x) ((double)(x)/RAND_MAX*rand())
#define ND 128 //分割数
#define NDX2 ND*2
#define INXY 400
#define PI 3.1415956535
#define RR 8.3144598
#define FF 9.64853415e+04                //[C/mol] ファラデー定数
#define VAVA 5300

int i, j;
int nd=ND, ndm=ND-1;        //組織の分割
int vava=VAVA+1;
double ch[ND][ND], Kh[ND][ND], Vh[ND][ND];
double K1, K2;
double V1, V2;
double XX;
double temp;
double E[ND][ND], W[ND][ND];
double delT[ND][ND], T[ND][ND];
double cVa[ND][ND];
double VVa[VAVA];
int ndx2 = NDX2;
double s0h[NDX2][NDX2], VH[NDX2][NDX2];
int ndx2m = ndx2 - 1;
double Q01, Q1, Q2;
double xr[NDX2][NDX2], xi[NDX2][NDX2], xrf[NDX2], xif[NDX2];
double s[NDX2], c[NDX2];
int ik[NDX2];
double qs;
int ig = 8;


void datin();
void datin2();
void shokiha_V();
void shokiha_S();
void table();
void fft();
void rcfft();
double G_T_V(double temp);
void datsave_c();
void datsave_V();
void datsave_E();
void datsave_W();
void datsave_T();
void datsave_delT();
void datsave_yVa();

void graph_c();
void graph_V();
void graph_W();
void graph_delT();
void graph_T();
void graph_yVa();

int main(void){
	cout << "main start!" << endl;

	int im, ip, jm, jp;
	double K01, K02;
	double al, b1;
	double V, V_E, V_W, V_N, V_S; 	//差分ブロックにおいてVを中心に、その上下左右
	double K, K_e, K_w, K_n, K_s;
	double beta;
	double time1, time1max;
	double sum1, sum2;
	double Vhdx[ND][ND], Vhdy[ND][ND];
	double Cp, den;

	// FILE	*stream;

	// stream = fopen("test_c.bin", "w");
	// fclose(stream);
	// stream = fopen("test_delT.bin", "w");
	// fclose(stream);
	// stream = fopen("test_T.bin", "w");
	// fclose(stream);
	// stream = fopen("test_W.bin", "w");
	// fclose(stream);
	// stream = fopen("test_E.bin", "w");
	// fclose(stream);
	// stream = fopen("test_yVa.bin", "w");
	// fclose(stream);
	// stream = fopen("test_V.bin", "w");
	// fclose(stream);


	printf("timestep=  ");
	scanf(" %lf",&XX);

	K01=1.0;      //[S/m]=[1/(Ωm)]  ZrO2表面
	K02=1.0e-18;  //[S/m]=[1/(Ωm)]  真空の値を適当に設定
	K1=K01/K02;  K2=K02/K02;   //無次元化

	al=3.0;            //計算領域の１辺(μm)
	al=1.0e-06*al;     //計算領域の１辺(m)
	b1=al/(double)nd;  //差分１ブロックのサイズ

	beta=1.5;
	V1=0.01;       //[V]低電位
	V2=0.3;      //[V]高電位
	double Vav=0.5*(V1+V2);

	time1=0.0;
	time1max=1.0e7;

	sum1=0.0;

	den=6.05e3;
	Cp=0.66e3*den;

	//設定電荷[C/(s・m^3)]
	//Q01=1.0e+7;		//初期設定値
	Q01=1.43e+07;		//換算値(No.800組織)
	//Q01=2.92e+07;		//換算値(No.1000組織)

	Q1=-Q01/(V2*K02/b1/b1);
	Q2=Q01/(V2*K02/b1/b1);

	//****** 場の読み込み ********************
	datin();
	datin2();
	shokiha_V();	//初期電位場の設定
	shokiha_S();

	//*** 電位計算 *******************************************
	double CH[ndx2][ndx2];
	int ii, jj;

	//濃度を4倍に拡張する
	for(i=0;i<=ndx2m;i++){
		for(j=0;j<=ndx2m;j++){
			if(i<nd){ii=i;} else{ii=ndx2m-i;}
			if(j<nd){jj=j;} else{jj=ndx2m-j;}
			CH[i][j]=ch[ii][jj];
		}
	}

	//濃度平均
	double c_0;
	sum1=0.0; for(i=0;i<=ndx2m;i++){ for(j=0;j<=ndx2m;j++){ sum1+=CH[i][j]; } }
    c_0=sum1/ndx2/ndx2;

	//導電率
	double K0;
	double dKh[ndx2][ndx2];
	K0=K1*c_0+K2*(1.0-c_0);
	for(i=0;i<=ndx2m;i++){
		for(j=0;j<=ndx2m;j++){
			dKh[i][j]=(K1-K2)*(CH[i][j]-c_0);
		}
	}

	//**** 電荷の発生（消滅）量のフ－リエ変換 s0h ---> s0qh1 ********************************

	double s0qrh1[ndx2][ndx2], s0qih1[ndx2][ndx2];
	double qs;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			xr[i][j]=s0h[i][j];  xi[i][j]=0.0;
		}
	}
	qs=-1.; rcfft();
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			s0qrh1[i][j]=xr[i][j];  s0qih1[i][j]=xi[i][j];
		}
	}

	//収束計算
	int ief, loopief=300;
	double VH2[ndx2][ndx2];
	double J1[ndx2][ndx2], J2[ndx2][ndx2];
	double a1_qrh1[ndx2][ndx2], a1_qih1[ndx2][ndx2];
	double a2_qrh1[ndx2][ndx2], a2_qih1[ndx2][ndx2];

	for(ief=0;ief<=loopief;ief++){

		for(i=0;i<=ndx2m;i++){
			for(j=0;j<=ndx2m;j++){
				VH2[i][j]=VH[i][j];//補助配列にコピー
			}
		}

        //*** 電位勾配場の計算(符号に注意) **********
		for(i=0;i<=ndx2m;i++){
			for(j=0;j<=ndx2m;j++){
				ip=i+1; im=i-1; jp=j+1; jm=j-1;
				if(i==ndx2m){ip=0;} 	if(i==0){im=ndx2m;}
				if(j==ndx2m){jp=0;}   if(j==0){jm=ndx2m;}

				J1[i][j]=0.5*(VH[ip][j]-VH[im][j]);
				J2[i][j]=0.5*(VH[i][jp]-VH[i][jm]);
			}
		}

        //**** 電位勾配*導電率の差のフ－リエ変換（dKh[][]*J1[][] ---> a1_qrh1） ********************************
		for(i=0;i<=ndx2m;i++){
			for(j=0;j<=ndx2m;j++){
				xr[i][j]=dKh[i][j]*J1[i][j]; xi[i][j]=0.0;
			}
		}
		qs=-1.; rcfft();
		for(i=0;i<=ndx2m;i++){
			for(j=0;j<=ndx2m;j++){
				a1_qrh1[i][j]=xr[i][j]; a1_qih1[i][j]=xi[i][j];
			}
		}
		//a1_qrh1[0][0]=a1_qih1[0][0]=0.;

        //**** 電位勾配*導電率の差のフ－リエ変換 (dKh[][]*J2[][] ---> a2_qrh2） ********************************
		for(i=0;i<=ndx2m;i++){
			for(j=0;j<=ndx2m;j++){
				xr[i][j]=dKh[i][j]*J2[i][j]; xi[i][j]=0.0;
			}
		}
		qs=-1.; rcfft();
		for(i=0;i<=ndx2m;i++){
			for(j=0;j<=ndx2m;j++){
				a2_qrh1[i][j]=xr[i][j]; a2_qih1[i][j]=xi[i][j];
			}
		}
		//a2_qrh1[0][0]=a2_qih1[0][0]=0.;

        //***** 次のステップの電位場の計算 *************************************
		for(i=0;i<=ndx2m;i++){
			if(i<=nd-1){ii=i;}  if(i>=nd){ii=i-ndx2;}
			for(j=0;j<=ndx2m;j++){
				if(j<=nd-1){jj=j;}  if(j>=nd){jj=j-ndx2;}
					double kx=2.0*PI/(double)ndx2*(double)ii;
					double ky=2.0*PI/(double)ndx2*(double)jj;
					double alnn=sqrt(kx*kx+ky*ky);  if(alnn==0.){alnn=1.;}
					//xr[i][j]=( s0qrh1[i][j]	-( kx*(a1_qrh1[i][j]+a1_qih1[i][j])
					//												  +ky*(a2_qrh1[i][j]+a2_qih1[i][j]) ) )/(K0*alnn*alnn);
					//xi[i][j]=( s0qih1[i][j]	+( kx*(a1_qrh1[i][j]+a1_qih1[i][j])
					//												  +ky*(a2_qrh1[i][j]+a2_qih1[i][j]) ) )/(K0*alnn*alnn);
					xr[i][j]=( s0qrh1[i][j]-(kx*a1_qih1[i][j]+ky*a2_qih1[i][j]) )/(K0*alnn*alnn);
					xi[i][j]=( s0qih1[i][j]+(kx*a1_qrh1[i][j]+ky*a2_qrh1[i][j]) )/(K0*alnn*alnn);
			}
		}
		qs=1.; rcfft();
		for(i=0;i<=ndx2m;i++){
			for(j=0;j<=ndx2m;j++){
				VH[i][j]=xr[i][j];
			}
		}

		double w0 = 0.4;
		for(i=0;i<=ndx2m;i++){
			for(j=0;j<=ndx2m;j++){
				VH[i][j]=w0*VH[i][j]+(1.0-w0)*VH2[i][j];//重み付き平均
			}
		}


		sum1=sum2=0.0;  for(j=0;j<=ndm;j++){ sum1+=VH[1][j];  sum2+=VH[ndm-1][j]; }
		double Vx1=sum1/nd;  double Vx2=sum2/nd;
		for(j=0;j<=ndx2m;j++){
			VH[0][j]=VH[1][j]=VH[ndx2m][j]=VH[ndx2m-1][j]=Vx1;
			VH[ndm][j]=VH[ndm-1][j]=VH[nd][j]=VH[nd+1][j]=Vx2;
		}

		// graph_V();
		// printf("ief= %d \n", ief);

	}//iefのloop

	//電位勾配場の計算
	for(i=0;i<=ndx2m;i++){
		for(j=0;j<=ndx2m;j++){
			VH[i][j]=VH[i][j]*V2+Vav;//Vの次元を戻し、平均電位をVavにシフト
		}
	}

	for(i=0;i<=ndx2m;i++){
		for(j=0;j<=ndx2m;j++){
			ip=i+1; im=i-1; jp=j+1; jm=j-1;
			if(i==ndx2m){ip=0;}   if(i==0){im=ndx2m;}
			if(j==ndx2m){jp=0;}   if(j==0){jm=ndx2m;}

			double fact1=K01*CH[i][j]+K02*(1.0-CH[i][j]);
			J1[i][j]=-fact1*0.5*(VH[ip][j]-VH[im][j]);
			J2[i][j]=-fact1*0.5*(VH[i][jp]-VH[i][jm]);
		}
	}

	sum1=sum2=0.0;  for(j=0;j<=ndm;j++){ sum1+=VH[1][j];  sum2+=VH[ndm-1][j]; }
	double Vx1=sum1/nd;  double Vx2=sum2/nd;
	printf("Vx1, Vx2, ΔVx= %f  %f  %f \n", Vx1, Vx2, Vx2-Vx1);
	printf("V1, V2, ΔV= %f  %f  %f \n", V1, V2, V2-V1);

	double Q00=Q01*(V2-V1)/(Vx2-Vx1);
	printf("Q01, Q00= %e  %e \n", Q01, Q00);

	//*** [平均の導電率の計算] ***********************************************
	sum1=0.0; for(j=0;j<=ndm;j++){ sum1+=J1[nd/2][j]; }
	double Km=fabs( sum1/nd*(nd-3)/(Vx2-Vx1) );
	//Km=fabs( sum1/nd2*(nd2-3)/(V2-V1) );
    //平均の導電率の計算(3を引いているのは左右の１ブロックづつと、半ブロック２つ分)
	printf("Km= %e \n", Km);
	//***********************************************

	// VH->Vh
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			Vh[i][j] = VH[i][j];
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1;  im=i-1;   jp=j+1;  jm=j-1;
			if(i==ndm){ip=ndm;}  if(i==0){im=0;}
			if(j==ndm){jp=ndm;}  if(j==0) {jm=0;}
			Vhdx[i][j] = (Vh[ip][j] - Vh[im][j]) / 2.0 / b1;
			Vhdy[i][j] = (Vh[i][jp] - Vh[i][jm]) / 2.0 / b1;
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			E[i][j] = sqrt(Vhdx[i][j] * Vhdx[i][j] + Vhdy[i][j] * Vhdy[i][j]);
			W[i][j]=E[i][j]*E[i][j]*Kh[i][j]*K02;

			delT[i][j]=W[i][j]/Cp*0.1;
			T[i][j]=delT[i][j]+1700.0;         //炉温1500K
			//printf("%lf ",T[i][j]);
		}
	}

	//*****************************************

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			cVa[i][j]=G_T_V(T[i][j]);
		}
		//printf("e");
	}

	// gwinsize(INXY,INXY); ginit(2); gsetorg(0,0);      //描画Window表示
	//graph_c();
	//graph_V();
	//graph_delT();
	//graph_yVa();

	// datsave_c();
	datsave_V();
	// datsave_E();
	// datsave_W();
	// datsave_T();
	// datsave_delT();
	// datsave_yVa();
	// if(keypress()){return 0;}
	cout << "main end" << endl;
	return 0;
}//main


//*********** 組織形態情報の入力 **************************
void datin(){
	cout << "datin" << endl;

	FILE   *datin0;

 	double ceM, ceP;
	double cmax, cmin;
	double time2;

    ceM=0.0001;  ceP=0.9999;

	datin0 = fopen("dc_fix/data/test.dat", "r");

	start: ;
	fscanf(datin0, "%lf\n", &time2);

    for(i=0;i<=ndm;i++){
        for(j=0;j<=ndm;j++){
            fscanf(datin0, "%lf  ", &ch[i][j]);
            if(ch[i][j]>=ceP){ch[i][j]=ceP;}
            if(ch[i][j]<=ceM){ch[i][j]=ceM;}
        }
    }

	cmax=0.0; cmin=1.0;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			if(cmax<ch[i][j]){cmax=ch[i][j];}
			if(cmin>ch[i][j]){cmin=ch[i][j];}
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ch[i][j]=(ch[i][j]-cmin)/(cmax-cmin);//規格化
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			Kh[i][j]=K2+(K1-K2)*ch[i][j];//導電率マトリックスの設定
		}
	}

	if (time2 != XX) {goto start;}
	fclose(datin0);

	cout << "datin end" << endl;
}


//*********** 組織形態情報の入力 **************************
void datin2(){
	cout << "datin2" << endl;

	FILE		*datin1;

	datin1 = fopen("dc_fix/data/temp3.dat", "r");
	//fscanf(datin0, "%lf", &time1);

	for(i=0;i<vava;i++){
		fscanf(datin1, "%lf\n", &VVa[i]);
	}

	fclose(datin1);

	cout << "datin2 end" << endl;
}

//************[初期電位場]*************************
void shokiha_V(){
	cout << "shokiha_V" << endl;

	int i, j;
	int ii, jj;
    srand(time(NULL)); // 乱数初期化

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			VH[i][j]=V1;
			//Vh[i][j]=0.00001*DRND(1);
		}
	}

	for(j=0;j<=ndm;j++){  VH[ndm][j]=VH[ndm-1][j]=V2;  VH[0][j]=VH[1][j]=V1; }

	for(i=0;i<=ndx2m;i++){
		for(j=0;j<=ndx2m;j++){
			if(i<nd){ii=i;} else{ii=ndx2m-i;}
			if(j<nd){jj=j;} else{jj=ndx2m-j;}
			VH[i][j]=VH[ii][jj];
		}
	}

	for(i=0;i<=ndx2m;i++){
		for(j=0;j<=ndx2m;j++){
			VH[i][j]=VH[i][j]/V2;
		}
	}

	cout << "shokiha_V end" << endl;
}

//**************************************************:
void shokiha_S(){
	int i, j;
	int ii, jj;
    srand(time(NULL)); // 乱数初期化

	for(i=0;i<=ndx2m;i++){
		for(j=0;j<=ndx2m;j++){
			s0h[i][j]=0.0;
		}
	}

	for(j=0;j<=ndm;j++){ s0h[ndm][j]=s0h[ndm-1][j]=Q2;  s0h[0][j]=s0h[1][j]=Q1; }

	for(i=0;i<=ndx2m;i++){
		for(j=0;j<=ndx2m;j++){
			if(i<nd){ii=i;} else{ii=ndx2m-i;}
			if(j<nd){jj=j;} else{jj=ndx2m-j;}
			s0h[i][j]=s0h[ii][jj];
		}
	}
}

//******************************************************************
double G_T_V(double temp){

	int T_left, T_right;
	double xT, alpha, y;

		xT=temp-700.0;
		T_left = (int)xT;
		T_right =(int)xT + 1;          //xの含まれる範囲の両端を表すインデックス
		if (T_left == 5300){ T_right = 5300; }    //データ点領域右端の補正
		alpha = xT - (double)T_left;              //α、範囲左側との距離
		y = VVa[T_left] * (1.0 - alpha) + VVa[T_right] * alpha;    //yの補間

	return(y);
}

//************ フェーズフィールド場のデータ保存サブルーチン *******************************
void datsave_c()
{
	FILE		*stream;		//ストリームのポインタ設定
	int 		i, j, k, p;		//整数

	stream = fopen("test_c.bin", "ab");	//書き込む先のファイルを追記方式でオープン
	//printf("time %lf, 炉温%lf, T0%lf\n", time1, T_ro, T0);						//計算カウント数の表示

		for (i = 0; i < nd; i++) {
			for (j = 0; j < nd; j++) {
				fwrite(&ch[i][j], sizeof(double), 1, stream);
			}
		}
		//fprintf(stream, "\n");	//改行の書き込み

	fclose(stream);					//ファイルをクローズ
}

//************ フェーズフィールド場のデータ保存サブルーチン *******************************
void datsave_V()
{
	FILE		*stream;		//ストリームのポインタ設定
	int 		i, j, k, p;		//整数

	stream = fopen("dc_fix/bin/test_V.bin", "ab");	//書き込む先のファイルを追記方式でオープン
	//printf("time %f\n", time1);						//計算カウント数の表示

		for (i = 0; i < nd; i++) {
			for (j = 0; j < nd; j++) {
				fwrite(&Vh[i][j], sizeof(double), 1, stream);
			}
		}
		//fprintf(stream, "\n");	//改行の書き込み

	fclose(stream);					//ファイルをクローズ
}

//************ フェーズフィールド場のデータ保存サブルーチン *******************************
void datsave_E()
{
	FILE		*stream;		//ストリームのポインタ設定
	int 		i, j, k, p;		//整数

	stream = fopen("test_E.bin", "ab");	//書き込む先のファイルを追記方式でオープン
	//printf("time %f\n", time1);						//計算カウント数の表示

		for (i = 0; i < nd; i++) {
			for (j = 0; j < nd; j++) {
				fwrite(&E[i][j], sizeof(double), 1, stream);
			}
		}
		//fprintf(stream, "\n");	//改行の書き込み

	fclose(stream);					//ファイルをクローズ
}
//************ フェーズフィールド場のデータ保存サブルーチン *******************************
void datsave_W()
{
	FILE		*stream;		//ストリームのポインタ設定
	int 		i, j, k, p;		//整数

	stream = fopen("test_W.bin", "ab");	//書き込む先のファイルを追記方式でオープン
	//printf("time %f\n", time1);						//計算カウント数の表示

		for (i = 0; i < nd; i++) {
			for (j = 0; j < nd; j++) {
				fwrite(&W[i][j], sizeof(double), 1, stream);
			}
		}
		//fprintf(stream, "\n");	//改行の書き込み

	fclose(stream);					//ファイルをクローズ
}
//************ フェーズフィールド場のデータ保存サブルーチン *******************************
void datsave_T()
{
	FILE		*stream;		//ストリームのポインタ設定
	int 		i, j, k, p;		//整数

	stream = fopen("test_T.bin", "ab");	//書き込む先のファイルを追記方式でオープン
	//printf("time %f\n", time1);						//計算カウント数の表示

		for (i = 0; i < nd; i++) {
			for (j = 0; j < nd; j++) {
				fwrite(&T[i][j], sizeof(double), 1, stream);
			}
		}
		//fprintf(stream, "\n");	//改行の書き込み

	fclose(stream);					//ファイルをクローズ
}
//************ フェーズフィールド場のデータ保存サブルーチン *******************************
void datsave_delT()
{
	FILE		*stream;		//ストリームのポインタ設定
	int 		i, j, k, p;		//整数

	stream = fopen("test_delT.bin", "ab");	//書き込む先のファイルを追記方式でオープン
	//printf("time %f\n", time1);						//計算カウント数の表示

		for (i = 0; i < nd; i++) {
			for (j = 0; j < nd; j++) {
				fwrite(&delT[i][j], sizeof(double), 1, stream);
			}
		}
		//fprintf(stream, "\n");	//改行の書き込み

	fclose(stream);					//ファイルをクローズ
}
//************ フェーズフィールド場のデータ保存サブルーチン *******************************
void datsave_yVa()
{
	FILE		*stream;		//ストリームのポインタ設定
	int 		i, j, k, p;		//整数

	stream = fopen("test_yVa.bin", "ab");	//書き込む先のファイルを追記方式でオープン
	//printf("time %f\n", time1);						//計算カウント数の表示

		for (i = 0; i < nd; i++) {
			for (j = 0; j < nd; j++) {
				fwrite(&cVa[i][j], sizeof(double), 1, stream);
			}
		}
		//fprintf(stream, "\n");	//改行の書き込み

	fclose(stream);					//ファイルをクローズ
}

//*******[組織の描画]**************************************************
// void graph_c()
// {
// 	int i, j, ii, jj, A,B;
// 	double col,col_R,col_G,col_B;
// 	double c, x, xmax, xmin, y, ymax, ymin, rad0;
// 	int ixmin=0, iymin=0, igx, igy, irad0;
// 	int ixmax=INXY;
// 	int iymax=INXY;

// 	gcls(); //画面クリア
// 	xmin=0.; xmax=1.;
// 	ymin=0.; ymax=1.;

// 	rad0=1./nd/2.;
// 	irad0=(ixmax-ixmin)/(xmax-xmin)*rad0+1;

// 	for(i=0;i<=nd;i++){
// 		for(j=0;j<=nd;j++){
// 			x=1./nd*i+rad0;  igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
// 			y=1./nd*j+rad0;  igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
// 			ii=i; jj=j;
// 			if(i==nd){ii=0;}  if(j==nd){jj=0;}

// 			col=1.0-ch[i][j];

// 			col_R=col_G=col_B=col;
// 			if(ch[i][j]<0.3){col_B=col_G=col_R=1.0;}
// 			gcolor((int)(255*col_R),(int)(255*col_G),(int)(255*col_B));
// 			grect(igx-irad0,igy-irad0,igx+irad0,igy+irad0);
// 		}
// 	}

// 	save_screen("c_1.bmp");
// 	//printf("hello");
// 	swapbuffers();     //画面スワップ
// }

//*******[組織の描画]**************************************************
// void graph_V()
// {
// 	int i, j, ii, jj;
// 	double col,col_R,col_G,col_B;
// 	double c, x, xmax, xmin, y, ymax, ymin, rad0;
// 	int ixmin=0, iymin=0, igx, igy, irad0;
// 	int ixmax=INXY;
// 	int iymax=INXY;
// 	double A,B,C,D,E,F;

// 	gcls(); //画面クリア

// 	A=178.0; B=102.0; C=255.0;   //低
// 	//A=255.0; B=255.0; C=255.0;   //低
// 	D=255.0; E=102.0; F=178.0;   //高

// 	xmin=0.; xmax=1.;
// 	ymin=0.; ymax=1.;

// 	rad0=1./nd/2.;
// 	irad0=(ixmax-ixmin)/(xmax-xmin)*rad0+1;

// 	for(i=0;i<=nd;i++){
// 		for(j=0;j<=nd;j++){
// 			x=1./nd*i+rad0;  igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
// 			y=1./nd*j+rad0;  igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
// 			ii=i; jj=j;
// 			if(i==nd){ii=0;}  if(j==nd){jj=0;}

// 			col=Vh[i][j]/0.3*255.0;

// 			col_R=A+col*(D-A)/255.0;
// 			col_G=B+col*(E-B)/255.0;
// 			col_B=C+col*(F-C)/255.0;

// 			if(ch[i][j]<0.3){col_B=col_G=col_R=255.0;}

// 			gcolor((int)(col_R),(int)(col_G),(int)(col_B));
// 			grect(igx-irad0,igy-irad0,igx+irad0,igy+irad0);
// 		}
// 	}
// 	//save_screen("V_500000.bmp");
// 	swapbuffers();     //画面スワップ
// }


//*******[組織の描画]**************************************************
// void graph_T()
// {
// 	int i, j, ii, jj;
// 	double col,col_R,col_G,col_B;
// 	double c, x, xmax, xmin, y, ymax, ymin, rad0;
// 	int ixmin=0, iymin=0, igx, igy, irad0;
// 	int ixmax=INXY;
// 	int iymax=INXY;
// 	double A,B,C,D,E,F;

// 	gcls(); //画面クリア

// 	D=255.0; E=0.0; F=0.0;   //高
// 	A=0.0; B=0.0; C=255.0;   //低
// 	//A=255.0; B=255.0; C=255.0;   //低

// 	xmin=0.; xmax=1.;
// 	ymin=0.; ymax=1.;

// 	rad0=1./nd/2.;
// 	irad0=(ixmax-ixmin)/(xmax-xmin)*rad0+1;

// 	for(i=0;i<=nd;i++){
// 		for(j=0;j<=nd;j++){
// 			x=1./nd*i+rad0;  igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
// 			y=1./nd*j+rad0;  igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
// 			ii=i; jj=j;
// 			if(i==nd){ii=0;}  if(j==nd){jj=0;}

// 			col=T[i][j]/1000.0*255.0;

// 			col_R=A+col*(D-A)/255.0;
// 			col_G=B+col*(E-B)/255.0;
// 			col_B=C+col*(F-C)/255.0;

// 			if(ch[i][j]<0.3){col_B=col_G=col_R=255.0;}

// 			gcolor((int)(col_R),(int)(col_G),(int)(col_B));
// 			grect(igx-irad0,igy-irad0,igx+irad0,igy+irad0);
// 		}
// 	}
// 	//save_screen("3.bmp");
// 	swapbuffers();     //画面スワップ
// }

//*******[組織の描画]**************************************************
// void graph_delT()
// {
// 	int i, j, ii, jj;
// 	double col,col_R,col_G,col_B;
// 	double c, x, xmax, xmin, y, ymax, ymin, rad0;
// 	int ixmin=0, iymin=0, igx, igy, irad0;
// 	int ixmax=INXY;
// 	int iymax=INXY;
// 	double A,B,C,D,E,F;

// 	gcls(); //画面クリア

// 	D=255.0; E=0.0; F=0.0;   //高
// 	//A=193.0; B=224.0; C=255.0;   //低
// 	A=255.0; B=255.0; C=255.0;   //低

// 	xmin=0.; xmax=1.;
// 	ymin=0.; ymax=1.;

// 	rad0=1./nd/2.;
// 	irad0=(ixmax-ixmin)/(xmax-xmin)*rad0+1;

// 	for(i=0;i<=nd;i++){
// 		for(j=0;j<=nd;j++){
// 			x=1./nd*i+rad0;  igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
// 			y=1./nd*j+rad0;  igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
// 			ii=i; jj=j;
// 			if(i==nd){ii=0;}  if(j==nd){jj=0;}

// 			col=delT[i][j]/1000.0*255.0;

// 			col_R=A+col*(D-A)/255.0;
// 			col_G=B+col*(E-B)/255.0;
// 			col_B=C+col*(F-C)/255.0;

// 			if(ch[i][j]<0.3){col_B=col_G=col_R=255.0;}

// 			gcolor((int)(col_R),(int)(col_G),(int)(col_B));
// 			grect(igx-irad0,igy-irad0,igx+irad0,igy+irad0);
// 		}
// 	}
// 	//save_screen("delT2_500000.bmp");
// 	swapbuffers();     //画面スワップ
// }

//*******[組織の描画]**************************************************
// void graph_yVa()
// {
// 	int i, j, ii, jj;
// 	double col,col_R,col_G,col_B;
// 	double c, x, xmax, xmin, y, ymax, ymin, rad0;
// 	int ixmin=0, iymin=0, igx, igy, irad0;
// 	int ixmax=INXY;
// 	int iymax=INXY;
// 	double A,B,C,D,E,F;

// 	gcls(); //画面クリア

// 	A=220.0; B=220.0; C=220.0;   //低
// 	//A=255.0; B=255.0; C=255.0;   //低
// 	// D=242.0; E=135.0; F=18.0;   //高
// 	D=255.0; E=0.0; F=127.0;   //高

// 	xmin=0.; xmax=1.;
// 	ymin=0.; ymax=1.;

// 	rad0=1./nd/2.;
// 	irad0=(ixmax-ixmin)/(xmax-xmin)*rad0+1;

// 	for(i=0;i<=nd;i++){
// 		for(j=0;j<=nd;j++){
// 			x=1./nd*i+rad0;  igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
// 			y=1./nd*j+rad0;  igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
// 			ii=i; jj=j;
// 			if(i==nd){ii=0;}  if(j==nd){jj=0;}

// 			col=(cVa[i][j]-0.015)/0.000001*255.0;
// 			if(col>255.0){col=255.0;}
// 			//printf("%lf ",col);
// 			col_R=A+col*(D-A)/255.0;
// 			col_G=B+col*(E-B)/255.0;
// 			col_B=C+col*(F-C)/255.0;

// 			if(ch[i][j]<0.3){col_B=col_G=col_R=255.0;}

// 			gcolor((int)(col_R),(int)(col_G),(int)(col_B));
// 			grect(igx-irad0,igy-irad0,igx+irad0,igy+irad0);
// 		}
// 	}
// 	save_screen("yVa_500000.bmp");
// 	swapbuffers();     //画面スワップ
// }

//******* Sin, Cos のテーブルおよびビット反転テーブルの設定 ***************
void table()
{
	int it, it1, it2, mc, mn;
	double q;

	q=2.0*PI/ndx2;
	for(it=0;it<=nd-1;it++){ c[it]=cos(q*it); s[it]=sin(q*it); }//Sin, Cos のテーブル

	ik[0]=0; mn=nd; mc=1;
	for(it1=1;it1<=ig;it1++){
		for(it2=0;it2<=mc-1;it2++){
			ik[it2+mc]=ik[it2]+mn;				//ビット反転テーブル
		}
		mn=mn/2; mc=2*mc;
	}
}

//********** １次元高速フーリエ変換 **************************************
void fft()
{
	int ix, ka, kb, l2, lf, mf, n2, nf;
	double tj, tr;

	l2=1;
	for(lf=1;lf<=ig;lf++){
		n2=nd/l2;
		for(mf=1;mf<=l2;mf++){
			for(nf=0;nf<=n2-1;nf++){
				ix=nf*l2;
				ka=nf+2*n2*(mf-1);
				kb=ka+n2;
				tr=xrf[ka]-xrf[kb];  					tj=xif[ka]-xif[kb];
				xrf[ka]=xrf[ka]+xrf[kb]; 			xif[ka]=xif[ka]+xif[kb];
				xrf[kb]=tr*c[ix]-tj*qs*s[ix];	xif[kb]=tj*c[ix]+tr*qs*s[ix];
			}
		}
		l2=l2*2;
	}
}

//************ ２次元高速フーリエ変換 ***********************************
void rcfft()
{
	int i, ic, ir, j;

	for(ir=0;ir<=ndx2m;ir++){
		for(ic=0;ic<=ndx2m;ic++){
			xrf[ic]=xr[ir][ic];	xif[ic]=xi[ir][ic];
		}
	fft();
		for(ic=0;ic<=ndx2m;ic++){
			xr[ir][ic]=xrf[ik[ic]];	xi[ir][ic]=xif[ik[ic]];
		}
	}
	for(ic=0;ic<=ndx2m;ic++){
		for(ir=0;ir<=ndx2m;ir++){
			xrf[ir]=xr[ir][ic];	xif[ir]=xi[ir][ic];
		}
	fft();
		for(ir=0;ir<=ndx2m;ir++){
			xr[ir][ic]=xrf[ik[ir]];	xi[ir][ic]=xif[ik[ir]];
		}
	}
	if(qs>0.){return;}
	for(i=0;i<=ndx2m;i++){
		for(j=0;j<=ndx2m;j++){
			xr[i][j]=xr[i][j]/ndx2/ndx2;	xi[i][j]=xi[i][j]/ndx2/ndx2;
		}
	}
}
