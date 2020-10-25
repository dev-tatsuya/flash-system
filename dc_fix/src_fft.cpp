//����
//�g�D�ω�����

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <chrono>

using namespace std::chrono;

#define DRND(x) ((double)(x)/RAND_MAX*rand())
#define ND 128 //������
#define NDX 256
#define NDY 128
#define IGX 8
#define IGY 7
#define PI 3.1415956535
#define RR 8.3144598
#define FF 9.64853415e+04                //[C/mol] �t�@���f�[�萔
#define VAVA 5300

int i, j;
int nd=ND, ndm=ND-1;        //�g�D�̕���
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
int nd2=ND/2, nd2m=ND/2-1;	//(�g�D�̕�����)�^2�F�t�|���G�ϊ����Ŏg�p
int nx2=NDX/2, nx2m=NDX/2-1;	//(�g�D�̕�����)�^2�F�t�|���G�ϊ����Ŏg�p
int ny2=NDY/2, ny2m=NDY/2-1;	//(�g�D�̕�����)�^2�F�t�|���G�ϊ����Ŏg�pv
double Q01, Q1, Q02, Q2, Q00; //�d��
double s0h[NDX][NDY];			//�d�׏�i�����Ə��Łj
double qs;							//�t�|���G�ϊ�(qs:-1)�Ƌt�t�|���G�ϊ�(qs:1)�̋��
double xi[NDX][NDY], xr[NDX][NDY];		//�t�|���G�ϊ��̎����Ƌ����Ɏg�p����z��
double ffrx[NDX], ffix[NDX];
double ffry[NDY], ffiy[NDY];
double s[ND],c[ND];			//sin��cos�̃e�[�u��
int ik[ND];							//�r�b�g���]�z��
double Vh2[NDX][NDY];			//�g�D���̓d�ʃf�|�^�⏕�z��
double J1[NDX][NDY], J2[NDX][NDY];	//����
double Vx1, Vx2;    //�d��[V]
double Vav;

void shokiha_S();				//�����d�׏�ݒ�T�u���|�`��
void table();				//sin��cos�̃e�[�u���쐬�T�u���|�`��
void fft_x();		  			//�P�����e�e�s
void fft_y();
void fft_1Dx();		  			//�P�����e�e�s
void fft_1Dy();		  			//�P�����e�e�s
void rcfft();				//�Q�����e�e�s

int main(void){
	printf("main start!\n");

	int im, ip, jm, jp;
	double K01, K02;
	double al, b1;
	double V, V_E, V_W, V_N, V_S; 	//�����u���b�N�ɂ�����V�𒆐S�ɁA���̏㉺���E
	double K, K_e, K_w, K_n, K_s;
	double beta;
	double time1, time1max;
	double sum1, sum2;
	double Vhdx[ND][ND], Vhdy[ND][ND];
	double Cp, den;

	int loopief, ief;
	double c_0;
	double K0, Km;  //�`����(����)
	double s0qrh1[NDX][NDY],	s0qih1[NDX][NDY];			//�g�D�̐U���z��
	double dKh[NDX][NDY];                             //�`����(�ϓ���)
	double a1_qrh1[NDX][NDY],	a1_qih1[NDX][NDY];		//dummy�z��
	double a2_qrh1[NDX][NDY],	a2_qih1[NDX][NDY];		//dummy�z��
	int ii, jj;
	double kx, ky, alnn;
	double w0 = 0.4;

	printf("timestep = ");
	scanf("%lf",&timestep);

	printf("loop(6)=  "); scanf(" %d",&loopief); //�����v�Z���[�v��

	K01=1.0;      //[S/m]=[1/(��m)]  ZrO2�\��
	K02=1.0e-18;  //[S/m]=[1/(��m)]  �^��̒l��K���ɐݒ�
	K1=K01/K02;  K2=K02/K02;   //��������

	al=3.0;            //�v�Z�̈�̂P��(��m)
	al=1.0e-06*al;     //�v�Z�̈�̂P��(m)
	b1=al/(double)nd;  //�����P�u���b�N�̃T�C�Y

	beta=1.5;
	V1=0.01;      //[V]��d��
	V2=0.3;       //[V]���d��

	time1=0.0;
	time1max=1.0e7;

	sum1=0.0;

	den=6.05e3;
	Cp=0.66e3*den;

	Q01=0;
	Q1=-Q01/(V2*K02/b1/b1);
	Q2=Q01/(V2*K02/b1/b1);
	Vav = 0.5*(V1+V2);

	//****** ��̓ǂݍ��� ********************
	datin();
	// datin2();
	shokiha_V();	//�����d�ʏ�̐ݒ�
	shokiha_S();
	table();

	datsave_c();
	datsave_V();

	//*** �J��Ԃ��v�Z�X�^�|�g *******************************************

	printf("potential calc start\n");
	auto start = system_clock::now();

	sum1=0.0; for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ sum1+=ch[i][j]; } }
    c_0=sum1/nd/nd;

	K0=K1*c_0+K2*(1.0-c_0);

	//TODO �g������
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

	//***** �����v�Z *******************************************************************************************
	for(ief=0;ief<=loopief;ief++){

		for(i=0;i<=ndxm;i++){
			for(j=0;j<=ndym;j++){
				Vh2[i][j]=Vh[i][j];//�⏕�z��ɃR�s�[
			}
		}

        //*** �d�ʌ��z��̌v�Z(�����ɒ���) **********
		for(i=0;i<=ndxm;i++){
			for(j=0;j<=ndym;j++){
				ip=i+1; im=i-1; jp=j+1; jm=j-1;
				if(i==ndxm){ip=0;} 	if(i==0){im=ndxm;}
				if(j==ndym){jp=0;}   if(j==0){jm=ndym;}

				J1[i][j]=0.5*(Vh[ip][j]-Vh[im][j]);
				J2[i][j]=0.5*(Vh[i][jp]-Vh[i][jm]);
			}
		}

		//**** �d�ʌ��z*���d���̍��̃t�|���G�ϊ��idKh[][]*J1[][] ---> a1_qrh1�j ********************************
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

		//**** �d�ʌ��z*���d���̍��̃t�|���G�ϊ� (dKh[][]*J2[][] ---> a2_qrh1�j ********************************
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

		//***** ���̃X�e�b�v�̓d�ʏ�̌v�Z *************************************
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
				Vh[i][j]=w0*Vh[i][j]+(1.0-w0)*Vh2[i][j];//�d�ݕt������
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
			Vh[i][j]=Vh[i][j]*V2+Vav;//V�̎�����߂��A���ϓd�ʂ�Vav�ɃV�t�g
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
    //���ς̓��d���̌v�Z(3�������Ă���͍̂��E�̂P�u���b�N�ÂƁA���u���b�N�Q��)
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
			T[i][j]=delT[i][j]+1700.0;         //�F��1500K
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


//*********** �g�D�`�ԏ��̓��� **************************
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
			ch[i][j]=(ch[i][j]-cmin)/(cmax-cmin);//�K�i��
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			Kh[i][j]=K2+(K1-K2)*ch[i][j];//���d���}�g���b�N�X�̐ݒ�
		}
	}

	if (time2 != timestep) {goto start;}
	fclose(datin0);
}


//*********** �g�D�`�ԏ��̓��� **************************
void datin2(){

	FILE		*datin1;

	datin1 = fopen("dc_fix/data/temp3.dat", "r");
	//fscanf(datin0, "%lf", &time1);

	for(i=0;i<vava;i++){
		fscanf(datin1, "%lf\n", &VVa[i]);
	}

	fclose(datin1);
}

//************[�����d�ʏ�]*************************
void shokiha_V(){

	int ii, jj;
  	srand(time(NULL)); // ����������

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

//************[�����d�׏�]*************************
void shokiha_S()
{
	int i, j;
	int ii, jj;
	srand(time(NULL)); // ����������

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
		T_right =(int)xT + 1;          //x�̊܂܂��͈̗͂��[��\���C���f�b�N�X
		if (T_left == 5300){ T_right = 5300; }    //�f�[�^�_�̈�E�[�̕␳
		alpha = xT - (double)T_left;              //���A�͈͍����Ƃ̋���
		y = VVa[T_left] * (1.0 - alpha) + VVa[T_right] * alpha;    //y�̕��

	return(y);
}

//************ �t�F�[�Y�t�B�[���h��̃f�[�^�ۑ��T�u���[�`�� *******************************
void datsave_c()
{
	FILE		*stream;		//�X�g���[���̃|�C���^�ݒ�
	int 		i, j, k, p;		//����

	stream = fopen("dc_fix/bin/test_c.bin", "ab");	//�������ސ�̃t�@�C����ǋL�����ŃI�[�v��
	//printf("time %lf, �F��%lf, T0%lf\n", time1, T_ro, T0);						//�v�Z�J�E���g���̕\��

		for (i = 0; i < nd; i++) {
			for (j = 0; j < nd; j++) {
				fwrite(&ch[i][j], sizeof(double), 1, stream);
			}
		}
		//fprintf(stream, "\n");	//���s�̏�������

	fclose(stream);					//�t�@�C�����N���[�Y
}

//************ �t�F�[�Y�t�B�[���h��̃f�[�^�ۑ��T�u���[�`�� *******************************
void datsave_V()
{
	FILE		*stream;		//�X�g���[���̃|�C���^�ݒ�
	int 		i, j, k, p;		//����

	stream = fopen("dc_fix/bin/test_V.bin", "ab");	//�������ސ�̃t�@�C����ǋL�����ŃI�[�v��
	//printf("time %f\n", time1);						//�v�Z�J�E���g���̕\��

		for (i = 0; i < nd; i++) {
			for (j = 0; j < nd; j++) {
				fwrite(&Vh[i][j], sizeof(double), 1, stream);
			}
		}
		//fprintf(stream, "\n");	//���s�̏�������

	fclose(stream);					//�t�@�C�����N���[�Y
}

//************ �t�F�[�Y�t�B�[���h��̃f�[�^�ۑ��T�u���[�`�� *******************************
void datsave_E()
{
	FILE		*stream;		//�X�g���[���̃|�C���^�ݒ�
	int 		i, j, k, p;		//����

	stream = fopen("test_E.bin", "ab");	//�������ސ�̃t�@�C����ǋL�����ŃI�[�v��
	//printf("time %f\n", time1);						//�v�Z�J�E���g���̕\��

		for (i = 0; i < nd; i++) {
			for (j = 0; j < nd; j++) {
				fwrite(&E[i][j], sizeof(double), 1, stream);
			}
		}
		//fprintf(stream, "\n");	//���s�̏�������

	fclose(stream);					//�t�@�C�����N���[�Y
}
//************ �t�F�[�Y�t�B�[���h��̃f�[�^�ۑ��T�u���[�`�� *******************************
void datsave_W()
{
	FILE		*stream;		//�X�g���[���̃|�C���^�ݒ�
	int 		i, j, k, p;		//����

	stream = fopen("test_W.bin", "ab");	//�������ސ�̃t�@�C����ǋL�����ŃI�[�v��
	//printf("time %f\n", time1);						//�v�Z�J�E���g���̕\��

		for (i = 0; i < nd; i++) {
			for (j = 0; j < nd; j++) {
				fwrite(&W[i][j], sizeof(double), 1, stream);
			}
		}
		//fprintf(stream, "\n");	//���s�̏�������

	fclose(stream);					//�t�@�C�����N���[�Y
}
//************ �t�F�[�Y�t�B�[���h��̃f�[�^�ۑ��T�u���[�`�� *******************************
void datsave_T()
{
	FILE		*stream;		//�X�g���[���̃|�C���^�ݒ�
	int 		i, j, k, p;		//����

	stream = fopen("test_T.bin", "ab");	//�������ސ�̃t�@�C����ǋL�����ŃI�[�v��
	//printf("time %f\n", time1);						//�v�Z�J�E���g���̕\��

		for (i = 0; i < nd; i++) {
			for (j = 0; j < nd; j++) {
				fwrite(&T[i][j], sizeof(double), 1, stream);
			}
		}
		//fprintf(stream, "\n");	//���s�̏�������

	fclose(stream);					//�t�@�C�����N���[�Y
}
//************ �t�F�[�Y�t�B�[���h��̃f�[�^�ۑ��T�u���[�`�� *******************************
void datsave_delT()
{
	FILE		*stream;		//�X�g���[���̃|�C���^�ݒ�
	int 		i, j, k, p;		//����

	stream = fopen("test_delT.bin", "ab");	//�������ސ�̃t�@�C����ǋL�����ŃI�[�v��
	//printf("time %f\n", time1);						//�v�Z�J�E���g���̕\��

		for (i = 0; i < nd; i++) {
			for (j = 0; j < nd; j++) {
				fwrite(&delT[i][j], sizeof(double), 1, stream);
			}
		}
		//fprintf(stream, "\n");	//���s�̏�������

	fclose(stream);					//�t�@�C�����N���[�Y
}
//************ �t�F�[�Y�t�B�[���h��̃f�[�^�ۑ��T�u���[�`�� *******************************
void datsave_yVa()
{
	FILE		*stream;		//�X�g���[���̃|�C���^�ݒ�
	int 		i, j, k, p;		//����

	stream = fopen("test_yVa.bin", "ab");	//�������ސ�̃t�@�C����ǋL�����ŃI�[�v��
	//printf("time %f\n", time1);						//�v�Z�J�E���g���̕\��

		for (i = 0; i < nd; i++) {
			for (j = 0; j < nd; j++) {
				fwrite(&cVa[i][j], sizeof(double), 1, stream);
			}
		}
		//fprintf(stream, "\n");	//���s�̏�������

	fclose(stream);					//�t�@�C�����N���[�Y
}

//************ �Q���������t�[���G�ϊ� ***********************************
void rcfft()
{
	int i, ic, ir, j;

	//x����FFT
	for(ir=0;ir<=ndym;ir++){ //y�𑖍�
		for(ic=0;ic<=ndxm;ic++){ //x�𑖍�
			ffrx[ic]=xr[ir][ic];	ffix[ic]=xi[ir][ic];
		}
		fft_1Dx();
		for(ic=0;ic<=ndym;ic++){ //x�𑖍�
			xr[ir][ic]=ffrx[ik[ic]];	xi[ir][ic]=ffix[ik[ic]];
		}
	}
	//y����FFT
	for(ic=0;ic<=ndxm;ic++){ //x�𑖍�
		for(ir=0;ir<=ndym;ir++){ //y�𑖍�
			ffry[ir]=xr[ir][ic];	ffiy[ir]=xi[ir][ic];
		}
		fft_1Dy();
		for(ir=0;ir<=ndym;ir++){ //y�𑖍�
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

//********** x�����P���������t�[���G�ϊ� **************************************
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

	//�����ݒ�p
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

	//�t�[���G�ϊ�
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

	//�t�[���G�ϊ���̏����ϊ�
	for(i=0;i<=NDX-1;i++){ffrx2[i]=ffrx[i]; ffix2[i]=ffix[i];}
	for(i=0;i<=NDX-1;i++){
		ffrx[i]=ffrx2[ik[i]]; ffix[i]=ffix2[ik[i]];
		if(qs==1.0){ffrx[i]=ffrx[i]/(double)NDX; ffix[i]=ffix[i]/(double)NDX;}
	}
}

//********** y�����P���������t�[���G�ϊ� **************************************
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

	//�����ݒ�p
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

	//�t�[���G�ϊ�
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

	//�t�[���G�ϊ���̏����ϊ�
	for(i=0;i<=NDY-1;i++){ffry2[i]=ffry[i]; ffiy2[i]=ffiy[i];}
	for(i=0;i<=NDY-1;i++){
		ffry[i]=ffry2[ik[i]]; ffiy[i]=ffiy2[ik[i]];
		if(qs==1.0){ffry[i]=ffry[i]/(double)NDY; ffiy[i]=ffiy[i]/(double)NDY;}
	}
}

//******* Sin, Cos �̃e�[�u������уr�b�g���]�e�[�u���̐ݒ� ***************
void table()
{
	int it, it1, it2, mc, mn;
	double q;

	q=2.0*PI/nd;
	for(it=0;it<=nd2-1;it++){ c[it]=cos(q*it); s[it]=sin(q*it); }//Sin, Cos �̃e�[�u��

	ik[0]=0; mn=nd2; mc=1;
	for(it1=1;it1<=ig;it1++){
		for(it2=0;it2<=mc-1;it2++){
			ik[it2+mc]=ik[it2]+mn;				//�r�b�g���]�e�[�u��
		}
		mn=mn/2; mc=2*mc;
	}
}
