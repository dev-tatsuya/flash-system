//直流
//全部無次元化
//熱拡散

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define DRND(x) ((double)(x)/RAND_MAX*rand())
#define ND 128 //分割数
#define RR 8.3144598      //気体定数 [J/Kmol]
#define FF 9.64853415e+04 //ファラデー定数 [C/mol]
#define KB 1.38064852e-23 //ボルツマン定数 [J/K]
#define NA 6.02214086e23  //アボガドロ数 [/mol]
#define VAVA 1300 //これはなんだ？

int nd=ND, ndm=ND-1; //組織の分割
int vava=VAVA+1;
int i, j; //index

double D[ND][ND]; //拡散係数
double delt; //時間刻み

double c2a;	        //平均組成(1:Cu,2:Co)
double time1;       //時間(ル－プ回数)
double c2h[ND][ND];	//組織内の濃度デ－タ配列
double XX;

double Vh[ND][ND]; //電位場
double V1,V2;
double Q[ND][ND]; //ジュール熱

double temp;
double Dm;

double E[ND][ND]; //電場
double T_ro; //初期温度
double T0;

double yVa, yVa_min;
double delT[ND][ND], T[ND][ND];
double cO2[ND][ND];
double cVa[ND][ND]; //平衡酸素空孔濃度
double cVay[ND][ND];

double sigma[ND][ND];
double VVa[VAVA]; //構造空孔濃度 なぜ1300個？
double jI[ND][ND]; //電流？

double ramuda[ND][ND]; //λ/(ρ×Cp) = α

double ceM, ceP;
double cmax, cmin;

double al; //計算領域
double b1; //差分ブロックサイズ

double sigma_0; //真空の電気伝導率
double e0; //電気伝導率
double Kh[ND][ND], sigsum, Wsum;
double Isum, sigmasum;

void set_init_conc(); //初期濃度波設定サブル－チン
void set_init_electric_potential(); //初期電位場の設定

void datsave_c();
void datsave_V();
void datsave_W();
void datsave_T();
void datsave_delT();
void datsave_yVa();
void datsave_si();
void datsave_D();
void datsave_I();

void datin();
void datin2();
void datin3();

double G_T_V(double temp);

int main(void){

	int k;
	int ip, im, jp, jm;

	double c2h2[ND][ND];
	double c2k[ND][ND];
	double kai2_chem[ND][ND], kai2_str, kai2_surf[ND][ND];

	double max_times; 	             //最大時間
	double amob_c22, amob_c[ND][ND]; //原子の易動度定数
	double c1, c2;                   //濃度
	double c2ddtt[ND][ND];			 //濃度の増分
	double kapa_c2;                  //濃度勾配エネルギ－定数
	double L12;				         //原子間相互作用パラメ－タ

	double D0, Q0; //拡散係数
	double yy;

	double sumc2, dc2a;
	double atom_n;

	double beta;
	double time2, time2max;
	double sumv1, sumv2; //電位の合計はどのように使用しているのか
	double V, V_E, V_W, V_N, V_S; //差分ブロックにおいてVを中心に、その上下左右
	double K, K_e, K_w, K_n, K_s;
	double Si_e, Si_w, Si_n, Si_s;

	double M[ND][ND]; //拡散の易動度
	double M_ion[ND][ND]; //酸素空孔の移動度

	double Vhdx[ND][ND], Vhdy[ND][ND]; //dV/dx, dV/dy
	double den[ND][ND]; //平均密度
	double Cp[ND][ND]; //密度×熱容量
	double n_el[ND][ND]; //単位体積あたりの電子数
	double n_ion[ND][ND]; //単位体積あたりの酸素イオン濃度
	double T1[ND][ND], T2[ND][ND]; //違いは？
	double ZZ; //YSZ分子量？にしては少し高め

	double Tddtt[ND][ND], Trddtt[ND][ND];
	double M0; //電子の移動度に関する定数
	double M1; //酸素空孔の移動度に関する定数
	double M_el[ND][ND]; //電子の移動度

	double Gc0;
	double cVa0; //熱平衡により生じる酸素空孔濃度
	double cVa1; //構造空孔濃度

	int interval;

	FILE *stream;

	// printf("炉温= ");
	// scanf("%lf",&T_ro);

	// printf("yy= ");
	// scanf("%lf",&yy);

	printf("save interval= ");
	scanf("%d",&interval);

	printf("max times= ");
	scanf("%lf",&max_times);

	T_ro=1650.0; //初期温度

	e0=1.6021766208e-19; //[c]電気素量
	T0=1200.0; //なに？どこからでてきた値？

	al=3.0;           //計算領域の１辺(μm)
	al=1.0e-06*al;    //計算領域の１辺(m)
	b1=al/(double)nd; //差分１ブロックのサイズ

	delt=0.0005;
	//delt=0.00008;

	amob_c22=1.0;
	time1=0.0;
	// max_times=500000.0;
	// max_times=80000.;

	L12=1.85e4/RR/T0;             //無次元化
	kapa_c2=10.0e-15/b1/b1/RR/T0; //無次元化

	//D0=0.35e-6/b1/b1;
	D0=2.5e-6/b1/b1;   //[/s]無次元化
	// Q0=64.5e3/RR/T0;   //[J/mol]→無次元化
	Q0=0.67*e0/KB;
	Dm=D0*exp(-Q0/T0); //基準の拡散係数

	V1=0.3*e0/KB/T0;  //[V]=[J/C] 高電位→無次元化
	V2=0.01*e0/KB/T0; //[V] 低電位
	beta=1.5;

	time2=0.;        //スタ－ト時間(ル－プの回数)
	time2max=1.0e+7; //ループ数の上限打ち切り

	sumv1=0.0;

	// den=6.05e3;                           //YSZ密度[kg/m3]
	// Cp=0.66e3*den*b1*b1*b1/KB;         //[J/kgK]×[kg/m3]=[J/Km3]→無次元化
	ZZ=126.2977e-3;         //[kg/mol]

	M0=8.02e-2/b1/b1*e0/KB/T0/Dm; //[m2/Vs]→無次元化
	M1=1.0e-3/b1/b1*e0/KB/T0/Dm;

	sigma_0=1.0e-25*b1/KB/T0/Dm; //[1/Ωm]  真空の値を適当に設定 卒論ではe+25になっている<-これが間違いっぽい

	//ramuda=3.0*b1/Dm/KB/Cp; //[W/mK] 熱伝導率→無次元化

	set_init_conc();
	set_init_electric_potential();
	datin();
	//datin3();
	//printf("%e\n",b1*KB*T0*Dm);

	cVa0=G_T_V(T_ro)-0.015;
	cVa1=0.015;

	// ρ(r)*Cp(r)
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			den[i][j]=6.05e3*c2h[i][j]+0.232*(1.0-c2h[i][j]);
			//ρ(r)*Cp(r)のこと なぜ式(2.20)に従ってないのか？
			Cp[i][j]=0.66e3*den[i][j]*b1*b1*b1/KB*c2h[i][j]+1227.0*den[i][j]*b1*b1*b1/KB*(1.0-c2h[i][j]);         //[J/kgK]×[kg/m3]=[J/Km3]→無次元化
			//printf("%d %d %e\n",i,j,Cp[i][j]);
		}
	}

	// σ(r)
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			M_el[i][j]=M0*exp(-1.89*FF/RR/T_ro);    //[m2/Vs] 電子移動度
			n_el[i][j]=2.0*2.0*cVa0*NA*den[i][j]/ZZ*b1*b1*b1; //なぜ 2*2 なのか(式(2.41))

			M_ion[i][j]=M1*exp(-1.0*FF/RR/T_ro);
			n_ion[i][j]=cVa1*NA*den[i][j]/ZZ*b1*b1*b1;

			sigma[i][j]=(n_el[i][j]*M_el[i][j] + M_ion[i][j]*n_ion[i][j])*c2h[i][j] + sigma_0*(1.0-c2h[i][j]);     // 無次元化→e0/e0
		}
	}

	//SOR1回目****************************

	// V(r)
	for(time2=0.;time2<=time2max;time2+=1.){
		//[非線形ラプラス方程式の解法]
		for(i=2;i<=ndm-2;i++){//左端と右端を省く
			for(j=0;j<=ndm;j++){
				ip=i+1;  im=i-1;   jp=j+1;  jm=j-1;
				//printf("%d %d %lf\n",i, j, Vh[i][j]);

				if(i==ndm){ip=ndm-1;} 	if(i==0){im=1;}
				if(j==ndm){jp=0;}       if(j==0){jm=ndm;}

				V=Vh[i][j]; V_E=Vh[ip][j];  V_W=Vh[im][j];  V_N=Vh[i][jp];  V_S=Vh[i][jm];

				K=sigma[i][j];
				K_e=(sigma[ip][j]+K)/2.0;  K_w=(sigma[im][j]+K)/2.0;
				K_n=(sigma[i][jp]+K)/2.0;  K_s=(sigma[i][jm]+K)/2.0;
				Vh[i][j]=beta*(K_e*V_E+K_w*V_W+K_n*V_N+K_s*V_S)/(K_e+K_w+K_n+K_s)+(1.0-beta)*V;
			}
		}

		//境界条件再設定
		for(j=0;j<=ndm;j++){ Vh[ndm][j]=Vh[ndm-1][j]=V1;  Vh[0][j]=Vh[1][j]=V2; }

		sumv2=0.0;
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				sumv2+=Vh[i][j];
			}
		}
		if(fabs(sumv1-sumv2)<=5.0e-5){break;}
		else {sumv1=sumv2;}
	}

	//********************************************************************

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1;  im=i-1;   jp=j+1;  jm=j-1;
			if(i==ndm){ip=ndm;}  if(i==0){im=0;}
			if(j==ndm){jp=ndm;}  if(j==0) {jm=0;}
			Vhdx[i][j] = (Vh[ip][j] - Vh[im][j]) / 2.0;
			Vhdy[i][j] = (Vh[i][jp] - Vh[i][jm]) / 2.0;
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			E[i][j] = sqrt(Vhdx[i][j] * Vhdx[i][j] + Vhdy[i][j] * Vhdy[i][j]);
			Q[i][j]=sigma[i][j]*E[i][j]*E[i][j]*delt;

			delT[i][j]=Q[i][j]/Cp[i][j];
			T[i][j]=delT[i][j]+T_ro/T0;         //炉温 なにしてる？
			T1[i][j]=T[i][j];
			if(T[i][j]>3.0){T[i][j]=3.0;}
		}
	}

	//******************

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			cVa[i][j]=G_T_V(T[i][j]*T0)-0.015;
		}
	}
	//******************

	//再度、σ(r)計算
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			M_el[i][j]=M0*exp(-1.89*FF/RR/T[i][j]/T0);    //[m2/Vs] 電子移動度
			n_el[i][j]=2.0*2.0*cVa[i][j]*NA*den[i][j]/ZZ*b1*b1*b1;

			M_ion[i][j]=M1*exp(-1.0*FF/RR/T[i][j]/T0);
			n_ion[i][j]=cVa1*NA*den[i][j]/ZZ*b1*b1*b1;

			sigma[i][j]=(n_el[i][j]*M_el[i][j] + M_ion[i][j]*n_ion[i][j])*c2h[i][j] + sigma_0*(1.0-c2h[i][j]);
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			jI[i][j]=sigma[i][j]*E[i][j];
		}
	}
	//***********************************

	start: ;

	//拡散ポテンシャルと拡散の易動度の計算
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1;  im=i-1;   jp=j+1;  jm=j-1;
			if(i==ndm){ip=ndm;}  if(i==0){im=0;}
			if(j==ndm){jp=ndm;}  if(j==0) {jm=0;}

			c2=c2h[i][j];  c1=1.-c2; if(c1<=0.){c1=1.0e-06;}
			// D[i][j]=D0*exp(-Q0/T[i][j]/T0)/Dm*cVa[i][j]/0.015;           //無次元化
			D[i][j]=D0*exp(-Q0/T[i][j]/T0)/Dm;

			//--- 勾配拡散ポテンシャルの計算 --------------------------------------------
			kai2_surf[i][j]=-2.0*kapa_c2*(c2h[ip][j]+c2h[im][j]+c2h[i][jp]+c2h[i][jm]-4.0*c2);    //無次元化

			//--- 化学拡散ポテンシャルの計算 --------------------------------------------
			kai2_chem[i][j]=(log(c2)-log(c1))+L12*(c1-c2);                           //無次元化
			//kai2_chem[i][j]=L12*(c1-c2);

			//--- 拡散ポテンシャルの計算 --------------------------------------------
			c2k[i][j]=kai2_chem[i][j]+kai2_surf[i][j];
			M[i][j]=D[i][j]*c2h[i][j]*(1.0-c2h[i][j]);
		}
	}

	//濃度場の時間変化[vi, vii]
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1;  im=i-1;   jp=j+1;  jm=j-1;

			if(i==ndm){ip=ndm;}  if(i==0){im=0;}
			if(j==ndm){jp=ndm;}  if(j==0){jm=0;}

			c2ddtt[i][j]= 0.5*(M[i][j]+M[ip][j])*(c2k[ip][j]-c2k[i][j])
						- 0.5*(M[i][j]+M[im][j])*(c2k[i][j]-c2k[im][j])
						+ 0.5*(M[i][j]+M[i][jp])*(c2k[i][jp]-c2k[i][j])
						- 0.5*(M[i][j]+M[i][jm])*(c2k[i][j]-c2k[i][jm]);
			c2h2[i][j]=c2h[i][j]+c2ddtt[i][j]*delt;
			// c2h2[i][j]=c2h[i][j]+c2ddtt[i][j]*delt+4.0e-3*(2.*DRND(1.)-1.);
		}
	}

	//濃度場の収支の補正
	if((((int)(time1) % 200)==0)){
		sumc2=0.;  for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ sumc2+=c2h2[i][j]; } }
  		dc2a=sumc2/nd/nd-c2a;
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				c2h[i][j]=c2h2[i][j]-dc2a;
				if(c2h[i][j]<=0.){c2h[i][j]=0.00001;}
				if(c2h[i][j]>=1.){c2h[i][j]=0.99999;}
	    	}
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			c2h[i][j]=c2h2[i][j];
			if(c2h[i][j]<=0.){c2h[i][j]=0.00001;}
			if(c2h[i][j]>=1.){c2h[i][j]=0.99999;}
		}
	}

	//*****電位(SOR)***************************************
	for(time2=0.;time2<=time2max;time2+=1.){
		//非線形ラプラス方程式の解法
		for(i=2;i<=ndm-2;i++){//左端と右端を省く
			for(j=0;j<=ndm;j++){
				ip=i+1;  im=i-1;   jp=j+1;  jm=j-1;
				//printf("%d %d %lf\n",i, j, Vh[i][j]);

				if(i==ndm){ip=ndm-1;} 	if(i==0){im=1;}
				if(j==ndm){jp=0;}       if(j==0){jm=ndm;}

				V=Vh[i][j]; V_E=Vh[ip][j];  V_W=Vh[im][j];  V_N=Vh[i][jp];  V_S=Vh[i][jm];

				K=sigma[i][j];
				Si_e=(sigma[ip][j]+K)/2.0;  Si_w=(sigma[im][j]+K)/2.0;
				Si_n=(sigma[i][jp]+K)/2.0;  Si_s=(sigma[i][jm]+K)/2.0;
				Vh[i][j]=beta*(Si_e*V_E+Si_w*V_W+Si_n*V_N+Si_s*V_S)/(Si_e+Si_w+Si_n+Si_s)+(1.0-beta)*V;
			}
		}
		//境界条件再設定
		for(j=0;j<=ndm;j++){ Vh[ndm][j]=Vh[ndm-1][j]=V1;  Vh[0][j]=Vh[1][j]=V2; }

		sumv2=0.0;
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				sumv2+=Vh[i][j];
			}
		}
		//printf("%lf ",fabs(sumv1-sumv2));
		if(fabs(sumv1-sumv2)<=5.0e-5){break;}
		else {sumv1=sumv2;}
	}

	//*****電場,ジュール熱,温度**************************************************

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1;  im=i-1;   jp=j+1;  jm=j-1;
			if(i==ndm){ip=ndm;}  if(i==0){im=0;}
			if(j==ndm){jp=ndm;}  if(j==0) {jm=0;}
			Vhdx[i][j] = (Vh[ip][j] - Vh[im][j]) / 2.0;
			Vhdy[i][j] = (Vh[i][jp] - Vh[i][jm]) / 2.0;
		}
	}

	//熱拡散[viii]
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			den[i][j]=6.05e3*c2h[i][j]+0.232*(1.0-c2h[i][j]);
			Cp[i][j]=0.66e3*den[i][j]*b1*b1*b1/KB*c2h[i][j]+1227.0*den[i][j]*b1*b1*b1/KB*(1.0-c2h[i][j]); //[J/kgK]×[kg/m3]=[J/Km3]→無次元化
			ramuda[i][j]=3.0*b1/Dm/KB/Cp[i][j]*c2h[i][j] + 0.0891*b1/Dm/KB/Cp[i][j]*(1.0-c2h[i][j]);          //[Q/mK] 熱伝導率→無次元化
			//printf("%d %d %e %e\n",i,j,Cp[i][j], ramuda[i][j]);
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			E[i][j] = sqrt(Vhdx[i][j] * Vhdx[i][j] + Vhdy[i][j] * Vhdy[i][j]);
			Q[i][j]=E[i][j]*E[i][j]*sigma[i][j]*delt;
			delT[i][j]=Q[i][j]/Cp[i][j];
			T[i][j]=delT[i][j]+T1[i][j];
			T1[i][j]=T[i][j];
			T2[i][j]=T[i][j];
		}
	}
	//for(k=0;k<=9;k++){
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1;  im=i-1;   jp=j+1;  jm=j-1;
			if(i==ndm) {ip=ndm;}  if(i==0) {im=0;}
			if(j==ndm) {jp=ndm;}  if(j==0) {jm=0;}
			Tddtt[i][j]=T[ip][j]+T[im][j]+T[i][jp]+T[i][jm]-4.0*T[i][j];    //無次元
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1;  im=i-1;   jp=j+1;  jm=j-1;
			if(i==ndm){ip=ndm;}  if(i==0){im=0;}
			if(j==ndm){jp=ndm;}  if(j==0){jm=0;}

			Trddtt[i][j]= 0.5*(ramuda[i][j]+ramuda[ip][j])*(Tddtt[ip][j]-Tddtt[i][j])
						- 0.5*(ramuda[i][j]+ramuda[im][j])*(Tddtt[i][j]-Tddtt[im][j])
						+ 0.5*(ramuda[i][j]+ramuda[i][jp])*(Tddtt[i][jp]-Tddtt[i][j])
						- 0.5*(ramuda[i][j]+ramuda[i][jm])*(Tddtt[i][j]-Tddtt[i][jm]);

			T[i][j]=T2[i][j]+Trddtt[i][j]*delt;
			T2[i][j]=T[i][j];
			if(T[i][j]>3.0){T[i][j]=3.0;}
		}
	}
	//}//k

	//*****空孔濃度**************************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			cVa[i][j]=G_T_V(T[i][j]*T0)-0.015;
			//printf("%lf %.10e\n",T[i][j]*T0, cVa[i][j]);
		}
	}

	//printf("sigma");
	//*****電気伝導率******************************

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			M_el[i][j]=M0*exp(-1.89*FF/RR/T[i][j]/T0);    //[m2/Vs] 電子移動度
			n_el[i][j]=2.0*2.0*cVa[i][j]*NA*den[i][j]/ZZ*b1*b1*b1;

			M_ion[i][j]=M1*exp(-1.0*FF/RR/T[i][j]/T0);
			n_ion[i][j]=cVa1*NA*den[i][j]/ZZ*b1*b1*b1;
			sigma[i][j]=(n_el[i][j]*M_el[i][j] + M_ion[i][j]*n_ion[i][j])*c2h[i][j] + sigma_0*(1.0-c2h[i][j]);
		}
	}
	//printf("sigma");
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			jI[i][j]=sigma[i][j]*E[i][j];
		}
	}
	//printf("sigma");
		// sigmasum=0.0;
		// for(j=21;j<=107;j++){
		// 	sigmasum=sigmasum+sigma[64][j];
		// 	// sigsum=sigsum+sigma[i][32];
		// 	// Wsum=Wsum+Q[i][32];
		// 	//printf("%e %e\n",sigma[i][32],sigsum);
		// }
	//printf("sigma");

	if((((int)(time1) % interval)==0)) {
		datsave_c();
		datsave_V();
		datsave_W();
		datsave_T();
		datsave_delT();
		datsave_yVa();
		datsave_si();
		datsave_D();
		datsave_I();
	}

	time1=time1+1.0;
	if(time1<max_times){goto start;}//最大カウント数に到達したかどうかの判断
	end: ;
	return 0;
} //main

//********** [初期濃度場の設定] **********
void set_init_conc(){

	int k, l;
	int ir0, ir1, ix, iy, iddd;
	double r1, dy;
	double ceM, ceP;
	double rnd0, sum;
	double xmin, xmax, ymin, ymax;
 	srand(time(NULL)); // 乱数初期化

	//ceM=0.00001;  ceP=0.99999;
	ceM=0.07;  ceP=0.93; //なぜ0と1でないのか

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			c2h[i][j]=ceM;
		}
	}

	ir0=ND/6.+0.5;
	ir1=ND/6+1.0;

	ix=0; iy=2*ir0;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			iddd=(i-ix)*(i-ix)+(j-iy)*(j-iy);
			if(iddd<=(ir1*ir1)){c2h[i][j]=ceP;}
		}
	}

	ix=0; iy=4*ir0;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			iddd=(i-ix)*(i-ix)+(j-iy)*(j-iy);
			if(iddd<=(ir1*ir1)){c2h[i][j]=ceP;}
		}
	}

	ix=2*ir0; iy=2*ir0;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			iddd=(i-ix)*(i-ix)+(j-iy)*(j-iy);
			if(iddd<=(ir1*ir1)){c2h[i][j]=ceP;}
		}
	}

	ix=2*ir0; iy=4*ir0;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			iddd=(i-ix)*(i-ix)+(j-iy)*(j-iy);
			if(iddd<=(ir1*ir1)){c2h[i][j]=ceP;}
		}
	}

	ix=4*ir0; iy=2*ir0;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			iddd=(i-ix)*(i-ix)+(j-iy)*(j-iy);
			if(iddd<=(ir1*ir1)){c2h[i][j]=ceP;}
		}
	}

	ix=4*ir0; iy=4*ir0;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			iddd=(i-ix)*(i-ix)+(j-iy)*(j-iy);
			if(iddd<=(ir1*ir1)){c2h[i][j]=ceP;}
		}
	}

	ix=6*ir0; iy=2*ir0;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			iddd=(i-ix)*(i-ix)+(j-iy)*(j-iy);
			if(iddd<=(ir1*ir1)){c2h[i][j]=ceP;}
		}
	}

	ix=6*ir0; iy=4*ir0;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			iddd=(i-ix)*(i-ix)+(j-iy)*(j-iy);
			if(iddd<=(ir1*ir1)){c2h[i][j]=ceP;}
		}
	}

	sum=0.; for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ sum+=c2h[i][j]; } }
	c2a=sum/nd/nd;
	//printf("c2a = %f  \n", c2a);
}

//********** [初期電位場の設定] **********
void set_init_electric_potential(){

	double rnd0;

	ceM=0.00000001;  ceP=0.99999999;
	//cmax=1.0; cmin=0.0;

  	srand(time(NULL)); // 乱数初期化

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			Vh[i][j]=V2;
		}
	}

	for(j=0;j<=ndm;j++){
		Vh[ndm][j]=Vh[ndm-1][j]=V1;
		Vh[0][j]=Vh[1][j]=V2;
	}

	for(i=0;i<=ndm;i++){
    	for(j=0;j<=ndm;j++){
			if(c2h[i][j]>=ceP){c2h[i][j]=ceP;}
            if(c2h[i][j]<=ceM){c2h[i][j]=ceM;}
		}
	}

}

//********** [なにしてる？] **********
double G_T_V(double temp){

	int T_left, T_right;
	double xT, alpha, y;

	xT=temp-700.0;
	T_left = (int)xT;
	T_right =(int)xT + 1;//xの含まれる範囲の両端を表すインデックス
	if (T_left == 1300){ T_right = 1300; }//データ点領域右端の補正
	alpha = xT - (double)T_left;//α、範囲左側との距離
	y = VVa[T_left] * (1.0 - alpha) + VVa[T_right] * alpha;//yの補間
	//printf("%lf %d %d %lf %.10e %.10e %.10e\n",xT, T_left, T_right, alpha, VVa[T_left], VVa[T_right], y);

	return(y);
}

//*********** 組織形態情報の入力 **************************
void datin()
{
	FILE		*stream;
	int 		i;

	stream = fopen("dc_dev/data/temp10.dat", "r");
	//fscanf(stream, "%lf", &time1);

	for(i=0;i<vava;i++){
		fscanf(stream, "%lf\n", &VVa[i]);
	}

	//printf("hello");
	fclose(stream);
}

//*********** 組織形態情報の入力 **************************
// void datin()
// {
// 	FILE		*stream;
// 	int 		k;

// 	datin0 = fopen("test_c.bin", "rb");
// 	//fscanf(datin0, "%lf", &time1);
// 	for(k=0;k<1000;k++){
// 	fread(&c2h,nd*nd*sizeof(double), 1, stream);
//  	graph();
// 	}
// 	fclose(stream);
// }

//************ フェーズフィールド場のデータ保存サブルーチン *******************************
void datsave_c()
{
	FILE		*stream;		//ストリームのポインタ設定
	int 		i, j, k, p;		//整数

	stream = fopen("dc_dev/bin/test_c.bin", "ab");	//書き込む先のファイルを追記方式でオープン
	//printf("time %lf, 炉温%lf, T0%lf\n", time1, T_ro, T0);						//計算カウント数の表示
	printf("time %lf\n", time1);

	for (i = 0; i < nd; i++) {
		for (j = 0; j < nd; j++) {
			fwrite(&c2h[i][j], sizeof(double), 1, stream);
		}
	}
	//fprintf(stream, "\n");	//改行の書き込み

	fclose(stream);					//ファイルをクローズ
}

//************ フェーズフィールド場のデータ保存サブルーチン *******************************
void datsave_V()
{
	FILE		*stream;		//ストリームのポインタ設定
	double V00[ND][ND];

	stream = fopen("dc_dev/bin/test_V.bin", "ab");	//書き込む先のファイルを追記方式でオープン
	//printf("time %f\n", time1);						//計算カウント数の表示

		for (i = 0; i < nd; i++) {
			for (j = 0; j < nd; j++) {
				V00[i][j]=Vh[i][j]/e0*KB*T0;
				fwrite(&V00[i][j], sizeof(double), 1, stream);
			}
		}
		//fprintf(stream, "\n");	//改行の書き込み

	fclose(stream);					//ファイルをクローズ
}

//************ フェーズフィールド場のデータ保存サブルーチン *******************************
void datsave_W()
{
	FILE		*stream;		//ストリームのポインタ設定
	double W00[ND][ND];

	stream = fopen("dc_dev/bin/test_W.bin", "ab");	//書き込む先のファイルを追記方式でオープン
	//printf("time %f\n", time1);						//計算カウント数の表示

		for (i = 0; i < nd; i++) {
			for (j = 0; j < nd; j++) {
				W00[i][j]=Q[i][j]*T0/b1/b1/b1*KB;
				fwrite(&W00[i][j], sizeof(double), 1, stream);
			}
		}
		//fprintf(stream, "\n");	//改行の書き込み

	fclose(stream);					//ファイルをクローズ
}
//************ フェーズフィールド場のデータ保存サブルーチン *******************************
void datsave_T()
{
	FILE		*stream;		//ストリームのポインタ設定
	double T00[ND][ND];

	stream = fopen("dc_dev/bin/test_T.bin", "ab");	//書き込む先のファイルを追記方式でオープン
	//printf("time %f\n", time1);						//計算カウント数の表示

		for (i = 0; i < nd; i++) {
			for (j = 0; j < nd; j++) {
				T00[i][j]=T0*T[i][j];
				fwrite(&T00[i][j], sizeof(double), 1, stream);
			}
		}
		//fprintf(stream, "\n");	//改行の書き込み

	fclose(stream);					//ファイルをクローズ
}
//************ フェーズフィールド場のデータ保存サブルーチン *******************************
void datsave_delT()
{
	FILE		*stream;		//ストリームのポインタ設定
	double delT00[ND][ND];

	stream = fopen("dc_dev/bin/test_delT.bin", "ab");	//書き込む先のファイルを追記方式でオープン
	//printf("time %f\n", time1);						//計算カウント数の表示

		for (i = 0; i < nd; i++) {
			for (j = 0; j < nd; j++) {
				delT00[i][j]=delT[i][j]*T0*Dm/delt;
				fwrite(&delT00[i][j], sizeof(double), 1, stream);
			}
		}
		//fprintf(stream, "\n");	//改行の書き込み

	fclose(stream);					//ファイルをクローズ
}
//************ フェーズフィールド場のデータ保存サブルーチン *******************************
void datsave_yVa()
{
	FILE		*stream;		//ストリームのポインタ設定

	stream = fopen("dc_dev/bin/test_yVa.bin", "ab");	//書き込む先のファイルを追記方式でオープン
	//printf("time %f\n", time1);						//計算カウント数の表示

		for (i = 0; i < nd; i++) {
			for (j = 0; j < nd; j++) {
				cVay[i][j]=cVa[i][j]+0.015;
				if(cVay[i][j]<0.015){cVay[i][j]=0.015;}
				fwrite(&cVay[i][j], sizeof(double), 1, stream);
			}
		}
		//fprintf(stream, "\n");	//改行の書き込み

	fclose(stream);					//ファイルをクローズ
}
//************ フェーズフィールド場のデータ保存サブルーチン *******************************
void datsave_si()
{
	FILE		*stream;		//ストリームのポインタ設定
	double sig00[ND][ND];

	stream = fopen("dc_dev/bin/test_sig.bin", "ab");	//書き込む先のファイルを追記方式でオープン
	//printf("time %f\n", time1);						//計算カウント数の表示

		for (i = 0; i < nd; i++) {
			for (j = 0; j < nd; j++) {
				sig00[i][j]=sigma[i][j]/b1*KB*T0*Dm;
				fwrite(&sig00[i][j], sizeof(double), 1, stream);
			}
		}
		//fprintf(stream, "\n");	//改行の書き込み

	fclose(stream);					//ファイルをクローズ
}
//************ フェーズフィールド場のデータ保存サブルーチン *******************************
void datsave_D()
{
	FILE		*stream;		//ストリームのポインタ設定
	double D00[ND][ND];

	stream = fopen("dc_dev/bin/test_D.bin", "ab");	//書き込む先のファイルを追記方式でオープン
	//printf("time %f\n", time1);						//計算カウント数の表示

		for (i = 0; i < nd; i++) {
			for (j = 0; j < nd; j++) {
				D00[i][j]=D[i][j]*b1*b1*Dm;
				fwrite(&D00[i][j], sizeof(double), 1, stream);
			}
		}
		//fprintf(stream, "\n");	//改行の書き込み

	fclose(stream);					//ファイルをクローズ
}

//************ フェーズフィールド場のデータ保存サブルーチン *******************************
void datsave_I()
{
	FILE		*stream;		//ストリームのポインタ設定
	double jI00[ND][ND];

	stream = fopen("dc_dev/bin/test_I.bin", "ab");	//書き込む先のファイルを追記方式でオープン

		for (i = 0; i < nd; i++) {
			for (j = 0; j < nd; j++) {
				jI00[i][j]=jI[i][j]*e0*Dm;
				fwrite(&jI00[i][j], sizeof(double), 1, stream);
			}
		}
		//fprintf(stream, "\n");	//改行の書き込み

	fclose(stream);					//ファイルをクローズ
}


//*********** 組織形態情報の入力 **************************
void datin2()
{
	FILE   *datin1;

 	double ceM, ceP;
	double cmax, cmin;
	double time3=-1.0;

    ceM=0.0001;  ceP=0.9999;

	datin1 = fopen("test_2.dat", "r");

	start: ;
	fscanf(datin1, "%lf\n", &time3);
	printf("%lf\n",time3);

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			fscanf(datin1, "%lf  ", &c2h[i][j]);
			if(c2h[i][j]>=ceP){c2h[i][j]=ceP;}
			if(c2h[i][j]<=ceM){c2h[i][j]=ceM;}
		}
	}

	cmax=0.0; cmin=1.0;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			if(cmax<c2h[i][j]){cmax=c2h[i][j];}
			if(cmin>c2h[i][j]){cmin=c2h[i][j];}
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			c2h[i][j]=(c2h[i][j]-cmin)/(cmax-cmin);//規格化
		}
	}

	printf("%lf\n",time3);
	if (time3 != XX) {printf("%lf\n",time3);goto start;}
	printf("%lf\n",time3);
	fclose(datin1);
	printf("%lf\n",time3);
}


//*********** 組織形態情報の入力 **************************
// void datin3()
// {
// 	FILE		*datin4;

// 	datin4 = fopen("test_4.dat", "r");
// 	//fscanf(datin0, "%lf", &time1);

// 	for(i=0;i<=ndm;i++){
// 		for(j=0;j<=ndm;j++){

// 			fscanf(datin4, "%lf ", &c2h[i][j]);
// 			//printf("%e ", Q[i][j]);

// 		}
// 	}
// 	//printf("hello");
// 	fclose(datin4);
// }
