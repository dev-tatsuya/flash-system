//繰り返し計算に重みを付ける場合

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "wingxa.h"

#define DRND(x) ((double)(x)/RAND_MAX*rand())
#define ND 128
#define IG 7

#define INX 400				//描画window１辺(x方向)のピクセルサイズ
#define INY 400				//描画window１辺(y方向)のピクセルサイズ

	int nd=ND, ndm=ND-1;		//計算領域の一辺の差分分割数，ND-1を定義
	int nd2=ND/2, nd2m=ND/2-1;	//(組織の分割数)／2：フ－リエ変換内で使用
	int nd4=ND/4;
	int ig=IG;
	double PI=3.141592, time1;	//π，計算カウント数
	double RR=8.3145;			//ガス定数
	double V1, V2, Vav;
	double Vx1, Vx2;
	double Q01, Q1, Q02, Q2, Q00;
		
	double ch[ND][ND];			//組織内の濃度デ－タ配列
	double Vh[ND][ND];			//組織内の電位デ－タ配列
	double Vh2[ND][ND];			//組織内の電位デ－タ補助配列
	double J1[ND][ND], J2[ND][ND];	//流束
	double s0h[ND][ND];			//電荷場（発生と消滅）
	
	double qs;							//フ－リエ変換(qs:-1)と逆フ－リエ変換(qs:1)の区別
	double xi[ND][ND], xr[ND][ND], xif[ND], xrf[ND];		//フ－リエ変換の実部と虚部に使用する配列
	double s[ND],c[ND];			//sinとcosのテーブル
	int ik[ND];							//ビット反転配列

	void shokiha_V();				//初期電位場設定サブル－チン
	void shokiha_S();				//初期電荷場設定サブル－チン
	void graph_V();					//電位場表示サブル－チン
	void graph_c();
	void graph_J();			//ベクトル場描画サブル－チン
	void datsave();			//濃度デ－タ保存サブル－チン
	void datload0();			//初期波読み込み用サブル－チン
	void datload();			//初期波読み込み用サブル－チン
	void table();				//sinとcosのテーブル作成サブル－チン
	void fft();					//１次元ＦＦＴ
	void rcfft();				//２次元ＦＦＴ

//********メインプログラム************************************
int main(void)
{
	int loopief, ief;
	
	double V;
	double c_0;
	double al;									//計算領域
	double b1;									//差分ブロックサイズ

	double s0qrh1[ND][ND],	s0qih1[ND][ND];			//組織の振幅配列
	double dKh[ND][ND];
	double a1_qrh1[ND][ND],	a1_qih1[ND][ND];		//dummy配列
	double a2_qrh1[ND][ND],	a2_qih1[ND][ND];		//dummy配列
	
	int   i, j, k, l, ii, jj, kk, iii, jjj;			//整数
	int   p, q, m, n;									//整数
	int   ip, im, jp, jm;							//整数
	
	double K01, K1;
	double K02, K2;
	double K0, Km;

	double nx, ny, nz, alnn;
	double kx, ky, kz;
	double sum1, sum2, fact1;
	double w0;

//********計算条件および物質定数の設定************************

	printf("loop(6)=  "); scanf(" %d",&loopief); //収束計算ループ回数
	//loopief=10;
	
	al=1.0e-02;						//計算領域一辺の長さ(m)
	b1=al/ND;							//差分ブロック一辺の長さ(m)

	V1=0.0;//低電位側の電位(V)
	V2=1.0;//高電位側の電位(V)
	Vav=0.5*(V1+V2);

	K01=1.0/1.0e-02; //[1/(Ωm)]  カーボンブラック (c)
	K02=1.0/1.0e+12; //[1/(Ωm)]  オレフィン系エラストマー  (1-c)

	//設定電荷[C/(s・m^3)]
	//Q01=1.0e+7;		//初期設定値
	Q01=1.43e+07;		//換算値(No.800組織)
	//Q01=2.92e+07;		//換算値(No.1000組織)

	w0=0.4;//新しい計算結果への重み、１つ前の計算結果の重みは1-w0
	//フィラー内での移動が速い場合には、w0の値を下げる必要がある。

//**** 無次元化 ********************
	K1=K01/K02;
	K2=K02/K02;

	Q1=-Q01/(V2*K02/b1/b1); 
	Q2=Q01/(V2*K02/b1/b1);

//****** 場の読み込みとsin,cosテーブルの設定 ********************
	datload();
	shokiha_S();
	shokiha_V();	//初期電位場の設定
	table();			//sin,cosテーブル

//******** 画像window ***********************************
	gwinsize(INX,INY); ginit(2); gsetorg(0,0);	//描画のwindow表示
	graph_c();//組織形態の描画
	//graph_V();//電位勾配

//******** スタート ***************************
//start: ;
	
//********* PFの平均値の算出 ******************************************************
	sum1=0.0; for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ sum1+=ch[i][j]; } }
  c_0=sum1/nd/nd;
	//printf("c_0= %f \n", c_0);

//******** 導電率 K (平均値と変動量) ******************************************************
	K0=K1*c_0+K2*(1.0-c_0);

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			dKh[i][j]=(K1-K2)*(ch[i][j]-c_0);
		}
	}

//**** 電荷の発生（消滅）量のフ－リエ変換 s0h ---> s0qh1 ********************************
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
	//s0qrh1[0][0]=s0qih1[0][0]=0.;

//***** 収束計算 *******************************************************************************************
	for(ief=0;ief<=loopief;ief++){

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			Vh2[i][j]=Vh[i][j];//補助配列にコピー
		}
	}

//*** 電位勾配場の計算(符号に注意) **********
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				ip=i+1; im=i-1; jp=j+1; jm=j-1;
				if(i==ndm){ip=0;} 	if(i==0){im=ndm;}
				if(j==ndm){jp=0;}   if(j==0){jm=ndm;}

				J1[i][j]=0.5*(Vh[ip][j]-Vh[im][j]);
				J2[i][j]=0.5*(Vh[i][jp]-Vh[i][jm]);
			}
		}

//**** 電位勾配*導電率の差のフ－リエ変換（dKh[][]*J1[][] ---> a1_qrh1） ********************************
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				xr[i][j]=dKh[i][j]*J1[i][j]; xi[i][j]=0.0;
			}
		}
		qs=-1.; rcfft();
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				a1_qrh1[i][j]=xr[i][j]; a1_qih1[i][j]=xi[i][j];
			}
		}
		//a1_qrh1[0][0]=a1_qih1[0][0]=0.;
	
//**** 電位勾配*導電率の差のフ－リエ変換 (dKh[][]*J2[][] ---> a2_qrh2） ********************************
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				xr[i][j]=dKh[i][j]*J2[i][j]; xi[i][j]=0.0;
			}
		}
		qs=-1.; rcfft();
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				a2_qrh1[i][j]=xr[i][j]; a2_qih1[i][j]=xi[i][j];
			}
		}
		//a2_qrh1[0][0]=a2_qih1[0][0]=0.;
	
//***** 次のステップの電位場の計算 *************************************
		for(i=0;i<=ndm;i++){
			if(i<=nd2-1){ii=i;}  if(i>=nd2){ii=i-nd;}
			for(j=0;j<=ndm;j++){
				if(j<=nd2-1){jj=j;}  if(j>=nd2){jj=j-nd;}
					kx=2.0*PI/(double)nd*(double)ii; 
					ky=2.0*PI/(double)nd*(double)jj;
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
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				Vh[i][j]=xr[i][j];
				
			}
		}

		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				Vh[i][j]=w0*Vh[i][j]+(1.0-w0)*Vh2[i][j];//重み付き平均
			}
		}

		sum1=sum2=0.0;  for(j=0;j<=nd2m;j++){ sum1+=Vh[1][j];  sum2+=Vh[nd2m-1][j]; }
		Vx1=sum1/nd2;  Vx2=sum2/nd2;
		for(j=0;j<=ndm;j++){
			Vh[0][j]=Vh[1][j]=Vh[ndm][j]=Vh[ndm-1][j]=Vx1;
			Vh[nd2m][j]=Vh[nd2m-1][j]=Vh[nd2][j]=Vh[nd2+1][j]=Vx2;
		}

		graph_V();
		printf("ief= %d \n", ief);

	}//iefのloop


//*** 電位勾配場の計算 *****************************************************
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				Vh[i][j]=Vh[i][j]*V2+Vav;//Vの次元を戻し、平均電位をVavにシフト
			}
		}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1; im=i-1; jp=j+1; jm=j-1;
			if(i==ndm){ip=0;} 	if(i==0){im=ndm;}
			if(j==ndm){jp=0;}   if(j==0){jm=ndm;}

			fact1=K01*ch[i][j]+K02*(1.0-ch[i][j]);
			J1[i][j]=-fact1*0.5*(Vh[ip][j]-Vh[im][j]);
			J2[i][j]=-fact1*0.5*(Vh[i][jp]-Vh[i][jm]);
		}
	}

	sum1=sum2=0.0;  for(j=0;j<=nd2m;j++){ sum1+=Vh[1][j];  sum2+=Vh[nd2m-1][j]; }
	Vx1=sum1/nd2;  Vx2=sum2/nd2;
	printf("Vx1, Vx2, ΔVx= %f  %f  %f \n", Vx1, Vx2, Vx2-Vx1);
	printf("V1, V2, ΔV= %f  %f  %f \n", V1, V2, V2-V1);

	Q00=Q01*(V2-V1)/(Vx2-Vx1);
	printf("Q01, Q00= %e  %e \n", Q01, Q00);

//*** 結果の表示 *****************************************************
	graph_J();//電流場
	//graph_V();//電位勾配場


//*** [平均の導電率の計算] ***********************************************
	sum1=0.0; for(j=0;j<=nd2m;j++){ sum1+=J1[nd4][j]; }
	Km=fabs( sum1/nd2*(nd2-3)/(Vx2-Vx1) );
	//Km=fabs( sum1/nd2*(nd2-3)/(V2-V1) );
  //平均の導電率の計算(3を引いているのは左右の１ブロックづつと、半ブロック２つ分)
	printf("Km= %e \n", Km);

//************************************************************
	if(keypress()){return 0;}//キー待ち状態
	//time1=time1+1.;  if(time1<time1max){goto start;}//最大カウント数に到達したかどうかの判断

	end:;
  return 0;
}
			
//************[初期電位場]*****************************************
void shokiha_V()
{
	int i, j;
	int ii, jj;
  srand(time(NULL)); // 乱数初期化

	//for(j=0;j<=nd2m;j++){  Vh[nd2m][j]=Vh[nd2m-1][j]=V2;  Vh[0][j]=Vh[1][j]=V1; }
	//for(i=0;i<=nd2m-1;i++){
	//	for(j=0;j<=nd2m;j++){
	//		Vh[i][j]=V1+(V2-V1)*(double)(i)/(nd2m-1);
	//	}
	//}

	for(i=0;i<=nd2m;i++){
		for(j=0;j<=nd2m;j++){
			Vh[i][j]=0.0;
			//Vh[i][j]=0.00001*DRND(1);
		}
	}
	for(j=0;j<=nd2m;j++){  Vh[nd2m][j]=Vh[nd2m-1][j]=V2;  Vh[0][j]=Vh[1][j]=V1; }

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			if(i<nd2){ii=i;} else{ii=ndm-i;}
			if(j<nd2){jj=j;} else{jj=ndm-j;}
			Vh[i][j]=Vh[ii][jj];
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			Vh[i][j]=Vh[i][j]/V2;
		}
	}

}

//************[電荷の発生と消滅]*****************************************
void shokiha_S()
{
	int i, j;
	int ii, jj;
  srand(time(NULL)); // 乱数初期化

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			s0h[i][j]=0.0;
		}
	}

	for(j=0;j<=nd2m;j++){  s0h[nd2m][j]=s0h[nd2m-1][j]=Q2;  s0h[0][j]=s0h[1][j]=Q1; }

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			if(i<nd2){ii=i;} else{ii=ndm-i;}
			if(j<nd2){jj=j;} else{jj=ndm-j;}
			s0h[i][j]=s0h[ii][jj];
		}
	}

}

//*******[電位場の描画]**************************************************
void graph_V()
{
	int i, j, ii, jj;
	int col_R, col_G, col_B;
	double col, fact1;
	double c, x, xmax, xmin, y, ymax, ymin, rad0;
	int ixmin=0, iymin=0, igx, igy, irad0;
	int ixmax=INX;
	int iymax=INY;

	gcls(); //画面クリア
	xmin=0.; xmax=1.;
	ymin=0.; ymax=1.;

	//printf("time %f\n",time1);
	rad0=1./nd/2.;
	irad0=(ixmax-ixmin)/(xmax-xmin)*rad0+1;

	for(i=0;i<=nd;i++){
		for(j=0;j<=nd;j++){
			x=1./nd*i+rad0;  igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
			y=1./nd*j+rad0;  igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
			ii=i; jj=j;  if(i==nd){ii=0;}  if(j==nd){jj=0;}
			col=(Vh[ii][jj]-Vh[0][1])/(Vh[nd2m][1]-Vh[0][1]);
			//col=(V2*Vh[ii][jj]-V2)/(V1-V2);
			if(col>=1.0){col=1.0;}  if(col<=0.0){col=0.0;}
			if(ch[i][j]>0.5){fact1=1.0;} else{fact1=0.0;}
			//col_R=(int)(255*col*fact1);
			col_R=(int)(255*col);
			col_G=(int)(255*col);
			col_B=(int)(255*col);
			gcolor(col_R,col_G,col_B);
			grect(igx-irad0,igy-irad0,igx+irad0,igy+irad0);
		}
	}

	swapbuffers();

} 

//*******[組織の描画]**************************************************
void graph_c()
{
	int i, j, ii, jj;
	double col;
	double c, x, xmax, xmin, y, ymax, ymin, rad0;
	int ixmin=0, iymin=0, igx, igy, irad0;
	int ixmax=INX;
	int iymax=INY;

	gcls(); //画面クリア
	xmin=0.; xmax=1.;
	ymin=0.; ymax=1.;

	printf("time %f\n",time1);
	rad0=1./nd/2.;
	irad0=(ixmax-ixmin)/(xmax-xmin)*rad0+1;

	for(i=0;i<=nd;i++){
		for(j=0;j<=nd;j++){
			x=1./nd*i+rad0;  igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
			y=1./nd*j+rad0;  igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
			ii=i; jj=j;  if(i==nd){ii=0;}  if(j==nd){jj=0;}
			col=1.-ch[ii][jj];
			if(col>=1.0){col=1.0;}  if(col<=0.0){col=0.0;}
			gcolor((int)(255*col),(int)(255*col),(int)(255*col));
			grect(igx-irad0,igy-irad0,igx+irad0,igy+irad0);
		}
	}
	swapbuffers();

} 

//******* ベクトル場の表示 ***************************************
void graph_J()
{
	int i, j, k, ii, jj, kk;//整数
	double col;//色
	double x, xmax, xmin, y, ymax, ymin, rad0, dia0;//規格化座標系の設定
	int ixmin=0, iymin=0, igx, igy, irad0;//スクリーン座標系の設定
	double x1, y1, x2, y2, x3, y3, x4, y4, fact1, fact2, th0;//矢印関連の座標設定
	int igx1, igy1, igx2, igy2, igx3, igy3, igx4, igy4;//矢印関連のスクリーン座標の設定
	int ixmax=INX, iymax=INY;//描画Window範囲

  gcls(); //画面クリア
	xmin=0.; xmax=1.; ymin=0.; ymax=1.;//描画領域（規格化されている）

	printf("time %f\n",time1);//計算カウント数の表示
	dia0=1./nd;  
	rad0=dia0/2.;               irad0=(ixmax-ixmin)/(xmax-xmin)*rad0+1;
	//差分ブロックの半分の長さ	//スクリーン座標系に変換（+1は整数化時の切捨て補正）

//**** 下地の場を薄くして描画 ***************************
	for(i=0;i<=nd;i++){
		for(j=0;j<=nd;j++){
			x=rad0+dia0*i;  igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
			y=rad0+dia0*j;  igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
			//座標計算			//スクリーン座標系に変換
			ii=i; jj=j;  if(i==nd){ii=0;}  if(j==nd){jj=0;}//周期的境界条件
			//col=ch[ii][jj];
			col=1.0-(1.0-ch[ii][jj])*0.5;//場の色をRGBにて設定（0.5をかけて薄くしている）
			if(col>=1.){col=1.;}  if(col<=0.){col=0.;}
			gcolor((int)(255*col),(int)(255*col),(int)(255*col));//色設定
			grect(igx-irad0,igy-irad0,igx+irad0,igy+irad0);//中塗り四角形の描画
		}
	}

//**** 流れ場を矢印にて描画 ***************************
	//fact1=0.01;//矢印の長さ
	fact1=0.1;//矢印の長さ
	fact2=0.2;//矢印先端位置の設定
	th0=0.5*3.14159/2.;//矢の部分の角度

	for(i=0;i<=nd;i+=4){
		for(j=0;j<=nd;j+=4){
			x=rad0+dia0*i;  igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
			y=rad0+dia0*j;  igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
			//座標計算			//スクリーン座標系に変換
			ii=i; jj=j; if(i==nd){ii=0;}  if(j==nd){jj=0;}//周期的境界条件

//--- 矢印の各部位の座標計算 ---
			x1=x-0.5*J1[ii][jj]*fact1; y1=y-0.5*J2[ii][jj]*fact1;
			igx1=(ixmax-ixmin)/(xmax-xmin)*(x1-xmin)+ixmin;
			igy1=(iymax-iymin)/(ymax-ymin)*(y1-ymin)+iymin;

			x2=x+0.5*J1[ii][jj]*fact1; y2=y+0.5*J2[ii][jj]*fact1;
			igx2=(ixmax-ixmin)/(xmax-xmin)*(x2-xmin)+ixmin;
			igy2=(iymax-iymin)/(ymax-ymin)*(y2-ymin)+iymin;

			x3=x+0.5*J1[ii][jj]*fact1*fact2*cos(th0)-0.5*J2[ii][jj]*fact1*fact2*sin(th0); 
			y3=y+0.5*J1[ii][jj]*fact1*fact2*sin(th0)+0.5*J2[ii][jj]*fact1*fact2*cos(th0);
			igx3=(ixmax-ixmin)/(xmax-xmin)*(x3-xmin)+ixmin;
			igy3=(iymax-iymin)/(ymax-ymin)*(y3-ymin)+iymin;

			x4=x+0.5*J1[ii][jj]*fact1*fact2*cos(-th0)-0.5*J2[ii][jj]*fact1*fact2*sin(-th0); 
			y4=y+0.5*J1[ii][jj]*fact1*fact2*sin(-th0)+0.5*J2[ii][jj]*fact1*fact2*cos(-th0);
			igx4=(ixmax-ixmin)/(xmax-xmin)*(x4-xmin)+ixmin;
			igy4=(iymax-iymin)/(ymax-ymin)*(y4-ymin)+iymin;
//----------------------

			//circle(igx2, igy2, 1, 0);  paint(igx2, igy2,0,0);
			gcolor(0,0,0);//色設定
			gline(igx1, igy1, igx2, igy2);//線を引く（矢印の部品）
			gline(igx2, igy2, igx3, igy3);//線を引く（矢印の部品）
			gline(igx2, igy2, igx4, igy4);//線を引く（矢印の部品）
		}
	}
	swapbuffers();//画面スワップ
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

//********** １次元高速フーリエ変換 **************************************
void fft()
{
	int ix, ka, kb, l2, lf, mf, n2, nf;
	double tj, tr;

	l2=1;
	for(lf=1;lf<=ig;lf++){
		n2=nd2/l2;
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

	for(ir=0;ir<=ndm;ir++){
		for(ic=0;ic<=ndm;ic++){
			xrf[ic]=xr[ir][ic];	xif[ic]=xi[ir][ic];
		}
	fft();
		for(ic=0;ic<=ndm;ic++){
			xr[ir][ic]=xrf[ik[ic]];	xi[ir][ic]=xif[ik[ic]];
		}
	}
	for(ic=0;ic<=ndm;ic++){
		for(ir=0;ir<=ndm;ir++){
			xrf[ir]=xr[ir][ic];	xif[ir]=xi[ir][ic];
		}
	fft();
		for(ir=0;ir<=ndm;ir++){
			xr[ir][ic]=xrf[ik[ir]];	xi[ir][ic]=xif[ik[ir]];
		}
	}
	if(qs>0.){return;}
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			xr[i][j]=xr[i][j]/nd/nd;	xi[i][j]=xi[i][j]/nd/nd;
		}
	}

}

//************[場デ－タの保存]************************************
void datsave()
{
	FILE		*stream;
	int 		i, j;

	stream = fopen("test.dat", "a");		//保存ファイル名をtest.datとしている。
	fprintf(stream, "%e\n", time1);		//時間の書き込み
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			fprintf(stream, "%e  %e  %e  ", Vh[i][j], J1[i][j], J2[i][j]);		//場の書き込み
		}
	}
	fprintf(stream, "\n");
	fclose(stream);
}
//*********** 組織形態情報の入力 **************************
void datload0()
{
	FILE		*datload0;
	int 		i, j, ip, im, jp, jm;
	int ii, jj;
	double cmax, cmin;

	cmax=0.0; cmin=1.0;
	//datload0 = fopen("ini_field.dat", "r");
	datload0 = fopen("test_sphere.dat", "r");
	//datload0 = fopen("test_double_sphere.dat", "r");

	fscanf(datload0, "%lf", &time1);
	for(i=0;i<=nd2m;i++){
		for(j=0;j<=nd2m;j++){
			fscanf(datload0, "%lf", &ch[i][j]);
			if(cmax<ch[i][j]){cmax=ch[i][j];}
			if(cmin>ch[i][j]){cmin=ch[i][j];}
		}
	}
	fclose(datload0);

	for(i=0;i<=nd2m;i++){
		for(j=0;j<=nd2m;j++){
			ch[i][j]=(ch[i][j]-cmin)/(cmax-cmin);//規格化
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			if(i<nd2){ii=i;} else{ii=ndm-i;}
			if(j<nd2){jj=j;} else{jj=ndm-j;}
			ch[i][j]=ch[ii][jj];
		}
	}

}

//*********** 組織形態情報の入力 **************************
void datload()
{
	FILE *datin0;
	int i, j, ip, im, jp, jm;
	int ii, jj;
	double cmax, cmin;

	cmax=0.0; cmin=1.0;
	//datin0 = fopen("test500.dat", "r");
	//datin0 = fopen("test600.dat", "r");
	//datin0 = fopen("test700.dat", "r");
	datin0 = fopen("test800.dat", "r");
	////datin0 = fopen("test900.dat", "r");
	//datin0 = fopen("test1000.dat", "r");
	////datin0 = fopen("test1100.dat", "r");
	////datin0 = fopen("test1200.dat", "r");
	////datin0 = fopen("test1300.dat", "r");
	//datin0 = fopen("test2000.dat", "r");
	//datin0 = fopen("test3000.dat", "r");
	//datin0 = fopen("test4000.dat", "r");

	fscanf(datin0, "%lf", &time1);
	for(i=0;i<=nd2m;i++){
		for(j=0;j<=nd2m;j++){
			fscanf(datin0, "%lf", &ch[i][j]);
			if(cmax<ch[i][j]){cmax=ch[i][j];}
			if(cmin>ch[i][j]){cmin=ch[i][j];}
		}
	}

	for(i=0;i<=nd2m;i++){
		for(j=0;j<=nd2m;j++){
			//ch[i][j]=0.0;
			ch[i][j]=(ch[i][j]-cmin)/(cmax-cmin);//規格化
			//ch[i][j]=1.-(ch[i][j]-cmin)/(cmax-cmin);//規格化
		}
	}

	fclose(datin0);

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			if(i<nd2){ii=i;} else{ii=ndm-i;}
			if(j<nd2){jj=j;} else{jj=ndm-j;}
			ch[i][j]=ch[ii][jj];
		}
	}

	//for(i=0;i<=ndm;i++){
	//	for(j=0;j<=ndm;j++){
	//		Kh[i][j]=K2+(K1-K2)*ch[i][j];//導電率マトリックスの設定
	//	}
	//}

}


