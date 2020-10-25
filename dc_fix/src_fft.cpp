//直流
//組織変化無し

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <chrono>

using namespace std::chrono;

#define DRND(x) ((double)(x)/RAND_MAX*rand())
#define ND 128 //分割数
#define NDX 256
#define NDY 128
#define IGX 8
#define IGY 7
#define PI 3.1415956535
#define RR 8.3144598
#define FF 9.64853415e+04                //[C/mol] ファラデー定数
#define VAVA 5300

int i, j;
int nd=ND, ndm=ND-1;        //組織の分割
int ndx=NDX, ndxm=NDX-1;
int ndy=NDY, ndym=NDY-1;
int vava=VAVA+1;
double ch[ND][ND], Kh[ND][ND], Vh[NDX][NDY];
double K1, K2;
double V1, V2;
double timestep;
double temp;
double E[ND][ND], W[ND][ND];
double delT[ND][ND], T[ND][ND];
double cVa[ND][ND];
double VVa[VAVA];

void datin();
void datin2();
void shokiha_V();
double G_T_V(double temp);
void datsave_c();
void datsave_V();
void datsave_E();
void datsave_W();
void datsave_T();
void datsave_delT();
void datsave_yVa();

int ig = 7;
int nd2=ND/2, nd2m=ND/2-1;	//(組織の分割数)／2：フ－リエ変換内で使用
int nx2=NDX/2, nx2m=NDX/2-1;	//(組織の分割数)／2：フ－リエ変換内で使用
int ny2=NDY/2, ny2m=NDY/2-1;	//(組織の分割数)／2：フ－リエ変換内で使用v
double Q01, Q1, Q02, Q2, Q00; //電荷
double s0h[NDX][NDY];			//電荷場（発生と消滅）
double qs;							//フ－リエ変換(qs:-1)と逆フ－リエ変換(qs:1)の区別
double xi[NDX][NDY], xr[NDX][NDY];		//フ－リエ変換の実部と虚部に使用する配列
double ffrx[NDX], ffix[NDX];
double ffry[NDY], ffiy[NDY];
double s[ND],c[ND];			//sinとcosのテーブル
int ik[ND];							//ビット反転配列
double Vh2[NDX][NDY];			//組織内の電位デ－タ補助配列
double J1[NDX][NDY], J2[NDX][NDY];	//流束
double Vx1, Vx2;    //電位[V]
double Vav;

void shokiha_S();				//初期電荷場設定サブル－チン
void table();				//sinとcosのテーブル作成サブル－チン
void fft_x();		  			//１次元ＦＦＴ
void fft_y();
void fft_1Dx();		  			//１次元ＦＦＴ
void fft_1Dy();		  			//１次元ＦＦＴ
void rcfft();				//２次元ＦＦＴ

int main(void){
	printf("main start!\n");

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

	int loopief, ief;
	double c_0;
	double K0, Km;  //伝導率(平均)
	double s0qrh1[NDX][NDY],	s0qih1[NDX][NDY];			//組織の振幅配列
	double dKh[NDX][NDY];                             //伝導率(変動量)
	double a1_qrh1[NDX][NDY],	a1_qih1[NDX][NDY];		//dummy配列
	double a2_qrh1[NDX][NDY],	a2_qih1[NDX][NDY];		//dummy配列
	int ii, jj;
	double kx, ky, alnn;
	double w0 = 0.4;

	printf("timestep = ");
	scanf("%lf",&timestep);

	printf("loop(6)=  "); scanf(" %d",&loopief); //収束計算ループ回数

	K01=1.0;      //[S/m]=[1/(Ωm)]  ZrO2表面
	K02=1.0e-18;  //[S/m]=[1/(Ωm)]  真空の値を適当に設定
	K1=K01/K02;  K2=K02/K02;   //無次元化

	al=3.0;            //計算領域の１辺(μm)
	al=1.0e-06*al;     //計算領域の１辺(m)
	b1=al/(double)nd;  //差分１ブロックのサイズ

	beta=1.5;
	V1=0.01;      //[V]低電位
	V2=0.3;       //[V]高電位

	time1=0.0;
	time1max=1.0e7;

	sum1=0.0;

	den=6.05e3;
	Cp=0.66e3*den;

	Q01=0;
	Q1=-Q01/(V2*K02/b1/b1);
	Q2=Q01/(V2*K02/b1/b1);
	Vav = 0.5*(V1+V2);

	//****** 場の読み込み ********************
	datin();
	// datin2();
	shokiha_V();	//初期電位場の設定
	shokiha_S();
	table();

	datsave_c();
	datsave_V();

	//*** 繰り返し計算スタ－ト *******************************************

	printf("potential calc start\n");
	auto start = system_clock::now();

	sum1=0.0; for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ sum1+=ch[i][j]; } }
    c_0=sum1/nd/nd;

	K0=K1*c_0+K2*(1.0-c_0);

	//TODO 拡張する
	double dKh2[ND][ND];
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			dKh2[i][j]=(K1-K2)*(ch[i][j]-c_0);
		}
	}

	for(i=0;i<=ndxm;i++){
		for(j=0;j<=ndym;j++){
			if(i<nx2){ii=i;} else{ii=ndxm-i;}
			dKh[i][j]=dKh2[ii][jj];
		}
	}

	for(i=0;i<=ndxm;i++){
		for(j=0;j<=ndym;j++){
			xr[i][j]=s0h[i][j];  xi[i][j]=0.0;
		}
	}
	qs=-1.; rcfft();
	for(i=0;i<=ndxm;i++){
		for(j=0;j<=ndym;j++){
			s0qrh1[i][j]=xr[i][j];  s0qih1[i][j]=xi[i][j];
		}
	}

	//***** 収束計算 *******************************************************************************************
	for(ief=0;ief<=loopief;ief++){

		for(i=0;i<=ndxm;i++){
			for(j=0;j<=ndym;j++){
				Vh2[i][j]=Vh[i][j];//補助配列にコピー
			}
		}

        //*** 電位勾配場の計算(符号に注意) **********
		for(i=0;i<=ndxm;i++){
			for(j=0;j<=ndym;j++){
				ip=i+1; im=i-1; jp=j+1; jm=j-1;
				if(i==ndxm){ip=0;} 	if(i==0){im=ndxm;}
				if(j==ndym){jp=0;}   if(j==0){jm=ndym;}

				J1[i][j]=0.5*(Vh[ip][j]-Vh[im][j]);
				J2[i][j]=0.5*(Vh[i][jp]-Vh[i][jm]);
			}
		}

		//**** 電位勾配*導電率の差のフ－リエ変換（dKh[][]*J1[][] ---> a1_qrh1） ********************************
		for(i=0;i<=ndxm;i++){
			for(j=0;j<=ndym;j++){
				xr[i][j]=dKh[i][j]*J1[i][j];
				xi[i][j]=0.0;
			}
		}

		qs=-1.; rcfft();
		for(i=0;i<=ndxm;i++){
			for(j=0;j<=ndym;j++){
				a1_qrh1[i][j]=xr[i][j];
				a1_qih1[i][j]=xi[i][j];
			}
		}
		//a1_qrh1[0][0]=a1_qih1[0][0]=0.;

		//**** 電位勾配*導電率の差のフ－リエ変換 (dKh[][]*J2[][] ---> a2_qrh1） ********************************
		for(i=0;i<=ndxm;i++){
			for(j=0;j<=ndym;j++){
				xr[i][j]=dKh[i][j]*J2[i][j];
				xi[i][j]=0.0;
			}
		}
		qs=-1.; rcfft();
		for(i=0;i<=ndxm;i++){
			for(j=0;j<=ndym;j++){
				a2_qrh1[i][j]=xr[i][j];
				a2_qih1[i][j]=xi[i][j];
			}
		}
		//a2_qrh1[0][0]=a2_qih1[0][0]=0.;

		//***** 次のステップの電位場の計算 *************************************
		for(i=0;i<=ndxm;i++){
			if(i<=nx2-1){ii=i;} else{ii=i-ndx;}
			for(j=0;j<=ndym;j++){
				if(j<=ny2-1){jj=j;} else{jj=j-ndy;}
				kx=2.0*PI/(double)ndx*(double)ii;
				ky=2.0*PI/(double)ndy*(double)jj;
				alnn=sqrt(kx*kx+ky*ky);  if(alnn==0.){alnn=1.;}
				//xr[i][j]=( s0qrh1[i][j]	-( kx*(a1_qrh1[i][j]+a1_qih1[i][j])
				//												  +ky*(a2_qrh1[i][j]+a2_qih1[i][j]) ) )/(K0*alnn*alnn);
				//xi[i][j]=( s0qih1[i][j]	+( kx*(a1_qrh1[i][j]+a1_qih1[i][j])
				//												  +ky*(a2_qrh1[i][j]+a2_qih1[i][j]) ) )/(K0*alnn*alnn);
				xr[i][j]=( s0qrh1[i][j]-(kx*a1_qih1[i][j]+ky*a2_qih1[i][j]) )/(K0*alnn*alnn);
				xi[i][j]=( s0qih1[i][j]+(kx*a1_qrh1[i][j]+ky*a2_qrh1[i][j]) )/(K0*alnn*alnn);
			}
		}
		qs=1.; rcfft();
		for(i=0;i<=ndxm;i++){
			for(j=0;j<=ndym;j++){
				Vh[i][j]=xr[i][j];
			}
		}

		for(i=0;i<=ndxm;i++){
			for(j=0;j<=ndym;j++){
				Vh[i][j]=w0*Vh[i][j]+(1.0-w0)*Vh2[i][j];//重み付き平均
			}
		}

		sum1=sum2=0.0;  for(j=0;j<=ny2m;j++){ sum1+=Vh[1][j];  sum2+=Vh[nx2m-1][j]; }
		Vx1=sum1/nx2;  Vx2=sum2/ny2;
		for(j=0;j<=ndym;j++){
			Vh[0][j]=Vh[1][j]=Vh[ndxm][j]=Vh[ndxm-1][j]=Vx1;
			Vh[nx2m][j]=Vh[nx2m-1][j]=Vh[nx2][j]=Vh[nx2+1][j]=Vx2;
		}

		// printf("ief= %d \n", ief);s
	}

	for(i=0;i<=ndxm;i++){
		for(j=0;j<=ndym;j++){
			Vh[i][j]=Vh[i][j]*V2+Vav;//Vの次元を戻し、平均電位をVavにシフト
		}
	}

	for(i=0;i<=ndxm;i++){
		for(j=0;j<=ndym;j++){
			ip=i+1; im=i-1; jp=j+1; jm=j-1;
			if(i==ndxm){ip=0;} 	if(i==0){im=ndxm;}
			if(j==ndym){jp=0;}   if(j==0){jm=ndym;}

			double fact1=K01*ch[i][j]+K02*(1.0-ch[i][j]);
			J1[i][j]=-fact1*0.5*(Vh[ip][j]-Vh[im][j]);
			J2[i][j]=-fact1*0.5*(Vh[i][jp]-Vh[i][jm]);
		}
	}

	sum1=sum2=0.0;  for(j=0;j<=ny2m;j++){ sum1+=Vh[1][j];  sum2+=Vh[nx2m-1][j]; }
	Vx1=sum1/nx2/ny2;  Vx2=sum2/ny2;
	printf("Vx1, Vx2, dVx= %f  %f  %f \n", Vx1, Vx2, Vx2-Vx1);
	printf("V1, V2, dV= %f  %f  %f \n", V1, V2, V2-V1);

	Q00=Q01*(V2-V1)/(Vx2-Vx1);
	printf("Q01, Q00= %e  %e \n", Q01, Q00);

	int nx4 = NDX/4;
	sum1=0.0; for(j=0;j<=ny2m;j++){ sum1+=J1[nx4][j]; }
	Km=fabs( sum1/ny2*(nx2-3)/(Vx2-Vx1) );
	//Km=fabs( sum1/nd2*(nd2-3)/(V2-V1) );
    //平均の導電率の計算(3を引いているのは左右の１ブロックづつと、半ブロック２つ分)
	printf("Km= %e \n", Km);

	printf("potential calc end\n");
	auto end = system_clock::now();
	auto msec = duration_cast<milliseconds>(end - start).count();
	printf("duration = %f sec.\n", msec/1000.0);

	//***********************************************

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

	// for(i=0;i<=ndm;i++){
	// 	for(j=0;j<=ndm;j++){
	// 		cVa[i][j]=G_T_V(T[i][j]);
	// 	}
	// 	//printf("e");
	// }

	datsave_c();
	datsave_V();
	// datsave_E();
	// datsave_W();
	// datsave_T();
	// datsave_delT();
	// datsave_yVa();
	// if(keypress()){return 0;}
	printf("main end\n");
	return 0;
}//main


//*********** 組織形態情報の入力 **************************
void datin(){

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

	if (time2 != timestep) {goto start;}
	fclose(datin0);
}


//*********** 組織形態情報の入力 **************************
void datin2(){

	FILE		*datin1;

	datin1 = fopen("dc_fix/data/temp3.dat", "r");
	//fscanf(datin0, "%lf", &time1);

	for(i=0;i<vava;i++){
		fscanf(datin1, "%lf\n", &VVa[i]);
	}

	fclose(datin1);
}

//************[初期電位場]*************************
void shokiha_V(){

	int ii, jj;
  	srand(time(NULL)); // 乱数初期化

	for(i=0;i<=nx2m;i++){
		for(j=0;j<=ny2m;j++){
			Vh[i][j]=V2;
		}
	}

	for(j=0;j<=ny2m;j++){  Vh[nx2m][j]=Vh[nx2m-1][j]=V2;  Vh[0][j]=Vh[1][j]=V1; }

	for(i=0;i<=ndxm;i++){
		for(j=0;j<=ndym;j++){
			if(i<nx2){ii=i;} else{ii=ndxm-i;}
			if(j<ny2){jj=j;} else{jj=ndym-j;}
			Vh[i][j]=Vh[ii][jj];
		}
	}

	for(i=0;i<=ndxm;i++){
		for(j=0;j<=ndym;j++){
			Vh[i][j]=Vh[i][j]/V2;
		}
	}
}

//************[初期電荷場]*************************
void shokiha_S()
{
	int i, j;
	int ii, jj;
	srand(time(NULL)); // 乱数初期化

	for(i=0;i<=ndxm;i++){
		for(j=0;j<=ndym;j++){
			s0h[i][j]=0.0;
		}
	}

	for(j=0;j<=ny2m;j++){
		s0h[nx2m][j]=s0h[nx2m-1][j]=Q2;
		s0h[0][j]=s0h[1][j]=Q1;
	}

	for(i=0;i<=ndxm;i++){
		for(j=0;j<=ndym;j++){
			if(i<nx2){ii=i;} else{ii=ndxm-i;}
			if(j<ny2){jj=j;} else{jj=ndym-j;}
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

	stream = fopen("dc_fix/bin/test_c.bin", "ab");	//書き込む先のファイルを追記方式でオープン
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

//************ ２次元高速フーリエ変換 ***********************************
void rcfft()
{
	int i, ic, ir, j;

	//x方向FFT
	for(ir=0;ir<=ndym;ir++){ //yを走査
		for(ic=0;ic<=ndxm;ic++){ //xを走査
			ffrx[ic]=xr[ir][ic];	ffix[ic]=xi[ir][ic];
		}
		fft_1Dx();
		for(ic=0;ic<=ndym;ic++){ //xを走査
			xr[ir][ic]=ffrx[ik[ic]];	xi[ir][ic]=ffix[ik[ic]];
		}
	}
	//y方向FFT
	for(ic=0;ic<=ndxm;ic++){ //xを走査
		for(ir=0;ir<=ndym;ir++){ //yを走査
			ffry[ir]=xr[ir][ic];	ffiy[ir]=xi[ir][ic];
		}
		fft_1Dy();
		for(ir=0;ir<=ndym;ir++){ //yを走査
			xr[ir][ic]=ffry[ik[ir]];	xi[ir][ic]=ffiy[ik[ir]];
		}
	}
	if(qs>0.){return;}
	for(i=0;i<=ndxm;i++){
		for(j=0;j<=ndym;j++){
			xr[i][j]=xr[i][j]/ndx/ndy;	xi[i][j]=xi[i][j]/ndx/ndy;
		}
	}
}

//********** x方向１次元高速フーリエ変換 **************************************
void fft_x()
{
	int ix, ka, kb, l2, lf, mf, n2, nf;
	double tj, tr;

	l2=1;
	for(lf=1;lf<=ig;lf++){
		n2=nx2/l2;
		for(mf=1;mf<=l2;mf++){
			for(nf=0;nf<=n2-1;nf++){
				ix=nf*l2;
				ka=nf+2*n2*(mf-1);
				kb=ka+n2;
				tr=ffrx[ka]-ffrx[kb];  					tj=ffix[ka]-ffix[kb];
				ffrx[ka]=ffrx[ka]+ffrx[kb]; 	  ffix[ka]=ffix[ka]+ffix[kb];
				ffrx[kb]=tr*c[ix]-tj*qs*s[ix]; ffix[kb]=tj*c[ix]+tr*qs*s[ix];
			}
		}
		l2=l2*2;
	}
}

void fft_1Dx(){
	int i, j, k;
	int ka, kb;
	int th;
	double tr, ti;
	double ffrx2[NDX], ffix2[NDX];

	//初期設定用
	static int ini_flag;
	static double s[NDX], c[NDX];
	static int ik[NDX];
	double q;
	int mc, mn;

	if(ini_flag==0){
		if((double)NDX!=pow(2.0,(double)IGX)){printf("ERROR \n"); return;}
		q=2.0*PI/(double)(NDX);
		for(i=0;i<=NDX-1;i++){ s[i]=sin(q*(double)i); c[i]=cos(q*(double)i);}

		ik[0]=0; mn=NDX/2; mc=1;
		for(i=1;i<=IGX;i++){
			for(j=0;j<=mc-1;j++){
				ik[j+mc]=ik[j]+mn;
			}
			mn=mn/2; mc=2*mc;
		}
		ini_flag=1;

		for(i=0;i<=NDX-1;i++){printf("ik[%d]=%d \n",i,ik[i]);}
		//getkey();
	}

	//フーリエ変換
	mn=NDX/2; mc=1;
	for(i=1;i<=IGX;i++){
		for(j=0;j<=mc-1;j++){
			for(k=0;k<=mn-1;k++){
				th=k*mc;
				ka=j*2*mn+k; kb=ka+mn;
				tr=ffrx[ka]-ffrx[kb]; ti=ffix[ka]-ffix[kb];
				ffrx[ka]=ffrx[ka]+ffrx[kb];
				ffix[ka]=ffix[ka]+ffix[kb];
				ffrx[kb]=tr*c[th]-qs*ti*s[th];
				ffix[kb]=ti*c[th]+qs*tr*s[th];
				//printf("th=%d ka%d kb=%d mn=%d mc=%d \n",th, ka, kb, mn, mc);
				//printf("fra=%f fia=%f frb=%f fib=%f \n",ffrx[ka],ffix[ka],ffrx[kb],ffix[kb]);
				//printf("tr=%f ti=%f \n",tr, ti);
			}
		}
		mn=mn/2; mc=2*mc;
	}

	//フーリエ変換後の順序変換
	for(i=0;i<=NDX-1;i++){ffrx2[i]=ffrx[i]; ffix2[i]=ffix[i];}
	for(i=0;i<=NDX-1;i++){
		ffrx[i]=ffrx2[ik[i]]; ffix[i]=ffix2[ik[i]];
		if(qs==1.0){ffrx[i]=ffrx[i]/(double)NDX; ffix[i]=ffix[i]/(double)NDX;}
	}
}

//********** y方向１次元高速フーリエ変換 **************************************
void fft_y()
{
	int ix, ka, kb, l2, lf, mf, n2, nf;
	double tj, tr;

	l2=1;
	for(lf=1;lf<=ig;lf++){
		n2=ny2/l2;
		for(mf=1;mf<=l2;mf++){
			for(nf=0;nf<=n2-1;nf++){
				ix=nf*l2;
				ka=nf+2*n2*(mf-1);
				kb=ka+n2;
				tr=ffry[ka]-ffry[kb];  					tj=ffiy[ka]-ffiy[kb];
				ffry[ka]=ffry[ka]+ffry[kb]; 	  ffiy[ka]=ffiy[ka]+ffiy[kb];
				ffry[kb]=tr*c[ix]-tj*qs*s[ix]; ffiy[kb]=tj*c[ix]+tr*qs*s[ix];
			}
		}
		l2=l2*2;
	}
}

void fft_1Dy(){
	int i, j, k;
	int ka, kb;
	int th;
	double tr, ti;
	double ffry2[NDY], ffiy2[NDY];

	//初期設定用
	static int ini_flag;
	static double s[NDY], c[NDY];
	static int ik[NDY];
	double q;
	int mc, mn;

	if(ini_flag==0){
		if((double)NDY!=pow(2.0,(double)IGY)){printf("ERROR \n"); return;}
		q=2.0*PI/(double)(NDY);
		for(i=0;i<=NDY-1;i++){ s[i]=sin(q*(double)i); c[i]=cos(q*(double)i);}

		ik[0]=0; mn=NDY/2; mc=1;
		for(i=1;i<=IGY;i++){
			for(j=0;j<=mc-1;j++){
				ik[j+mc]=ik[j]+mn;
			}
			mn=mn/2; mc=2*mc;
		}
		ini_flag=1;

		for(i=0;i<=NDY-1;i++){printf("ik[%d]=%d \n",i,ik[i]);}
		//getkey();
	}

	//フーリエ変換
	mn=NDY/2; mc=1;
	for(i=1;i<=IGY;i++){
		for(j=0;j<=mc-1;j++){
			for(k=0;k<=mn-1;k++){
				th=k*mc;
				ka=j*2*mn+k; kb=ka+mn;
				tr=ffry[ka]-ffry[kb]; ti=ffiy[ka]-ffiy[kb];
				ffry[ka]=ffry[ka]+ffry[kb];
				ffiy[ka]=ffiy[ka]+ffiy[kb];
				ffry[kb]=tr*c[th]-qs*ti*s[th];
				ffiy[kb]=ti*c[th]+qs*tr*s[th];
				//printf("th=%d ka%d kb=%d mn=%d mc=%d \n",th, ka, kb, mn, mc);
				//printf("fra=%f fia=%f frb=%f fib=%f \n",ffry[ka],ffiy[ka],ffry[kb],ffiy[kb]);
				//printf("tr=%f ti=%f \n",tr, ti);
			}
		}
		mn=mn/2; mc=2*mc;
	}

	//フーリエ変換後の順序変換
	for(i=0;i<=NDY-1;i++){ffry2[i]=ffry[i]; ffiy2[i]=ffiy[i];}
	for(i=0;i<=NDY-1;i++){
		ffry[i]=ffry2[ik[i]]; ffiy[i]=ffiy2[ik[i]];
		if(qs==1.0){ffry[i]=ffry[i]/(double)NDY; ffiy[i]=ffiy[i]/(double)NDY;}
	}
}

//******* Sin, Cos のテーブルおよびビット反転テーブルの設定 ***************
void table()
{
	int it, it1, it2, mc, mn;
	double q;

	q=2.0*PI/nd;
	for(it=0;it<=nd2-1;it++){ c[it]=cos(q*it); s[it]=sin(q*it); }//Sin, Cos のテーブル

	ik[0]=0; mn=nd2; mc=1;
	for(it1=1;it1<=ig;it1++){
		for(it2=0;it2<=mc-1;it2++){
			ik[it2+mc]=ik[it2]+mn;				//ビット反転テーブル
		}
		mn=mn/2; mc=2*mc;
	}
}
