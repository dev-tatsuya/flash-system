//�J��Ԃ��v�Z�ɏd�݂�t����ꍇ

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "wingxa.h"

#define DRND(x) ((double)(x)/RAND_MAX*rand())
#define ND 128
#define IG 7

#define INX 400				//�`��window�P��(x����)�̃s�N�Z���T�C�Y
#define INY 400				//�`��window�P��(y����)�̃s�N�Z���T�C�Y

	int nd=ND, ndm=ND-1;		//�v�Z�̈�̈�ӂ̍����������CND-1���`
	int nd2=ND/2, nd2m=ND/2-1;	//(�g�D�̕�����)�^2�F�t�|���G�ϊ����Ŏg�p
	int nd4=ND/4;
	int ig=IG;
	double PI=3.141592, time1;	//�΁C�v�Z�J�E���g��
	double RR=8.3145;			//�K�X�萔
	double V1, V2, Vav;
	double Vx1, Vx2;
	double Q01, Q1, Q02, Q2, Q00;
		
	double ch[ND][ND];			//�g�D���̔Z�x�f�|�^�z��
	double Vh[ND][ND];			//�g�D���̓d�ʃf�|�^�z��
	double Vh2[ND][ND];			//�g�D���̓d�ʃf�|�^�⏕�z��
	double J1[ND][ND], J2[ND][ND];	//����
	double s0h[ND][ND];			//�d�׏�i�����Ə��Łj
	
	double qs;							//�t�|���G�ϊ�(qs:-1)�Ƌt�t�|���G�ϊ�(qs:1)�̋��
	double xi[ND][ND], xr[ND][ND], xif[ND], xrf[ND];		//�t�|���G�ϊ��̎����Ƌ����Ɏg�p����z��
	double s[ND],c[ND];			//sin��cos�̃e�[�u��
	int ik[ND];							//�r�b�g���]�z��

	void shokiha_V();				//�����d�ʏ�ݒ�T�u���|�`��
	void shokiha_S();				//�����d�׏�ݒ�T�u���|�`��
	void graph_V();					//�d�ʏ�\���T�u���|�`��
	void graph_c();
	void graph_J();			//�x�N�g����`��T�u���|�`��
	void datsave();			//�Z�x�f�|�^�ۑ��T�u���|�`��
	void datload0();			//�����g�ǂݍ��ݗp�T�u���|�`��
	void datload();			//�����g�ǂݍ��ݗp�T�u���|�`��
	void table();				//sin��cos�̃e�[�u���쐬�T�u���|�`��
	void fft();					//�P�����e�e�s
	void rcfft();				//�Q�����e�e�s

//********���C���v���O����************************************
int main(void)
{
	int loopief, ief;
	
	double V;
	double c_0;
	double al;									//�v�Z�̈�
	double b1;									//�����u���b�N�T�C�Y

	double s0qrh1[ND][ND],	s0qih1[ND][ND];			//�g�D�̐U���z��
	double dKh[ND][ND];
	double a1_qrh1[ND][ND],	a1_qih1[ND][ND];		//dummy�z��
	double a2_qrh1[ND][ND],	a2_qih1[ND][ND];		//dummy�z��
	
	int   i, j, k, l, ii, jj, kk, iii, jjj;			//����
	int   p, q, m, n;									//����
	int   ip, im, jp, jm;							//����
	
	double K01, K1;
	double K02, K2;
	double K0, Km;

	double nx, ny, nz, alnn;
	double kx, ky, kz;
	double sum1, sum2, fact1;
	double w0;

//********�v�Z��������ѕ����萔�̐ݒ�************************

	printf("loop(6)=  "); scanf(" %d",&loopief); //�����v�Z���[�v��
	//loopief=10;
	
	al=1.0e-02;						//�v�Z�̈��ӂ̒���(m)
	b1=al/ND;							//�����u���b�N��ӂ̒���(m)

	V1=0.0;//��d�ʑ��̓d��(V)
	V2=1.0;//���d�ʑ��̓d��(V)
	Vav=0.5*(V1+V2);

	K01=1.0/1.0e-02; //[1/(��m)]  �J�[�{���u���b�N (c)
	K02=1.0/1.0e+12; //[1/(��m)]  �I���t�B���n�G���X�g�}�[  (1-c)

	//�ݒ�d��[C/(s�Em^3)]
	//Q01=1.0e+7;		//�����ݒ�l
	Q01=1.43e+07;		//���Z�l(No.800�g�D)
	//Q01=2.92e+07;		//���Z�l(No.1000�g�D)

	w0=0.4;//�V�����v�Z���ʂւ̏d�݁A�P�O�̌v�Z���ʂ̏d�݂�1-w0
	//�t�B���[���ł̈ړ��������ꍇ�ɂ́Aw0�̒l��������K�v������B

//**** �������� ********************
	K1=K01/K02;
	K2=K02/K02;

	Q1=-Q01/(V2*K02/b1/b1); 
	Q2=Q01/(V2*K02/b1/b1);

//****** ��̓ǂݍ��݂�sin,cos�e�[�u���̐ݒ� ********************
	datload();
	shokiha_S();
	shokiha_V();	//�����d�ʏ�̐ݒ�
	table();			//sin,cos�e�[�u��

//******** �摜window ***********************************
	gwinsize(INX,INY); ginit(2); gsetorg(0,0);	//�`���window�\��
	graph_c();//�g�D�`�Ԃ̕`��
	//graph_V();//�d�ʌ��z

//******** �X�^�[�g ***************************
//start: ;
	
//********* PF�̕��ϒl�̎Z�o ******************************************************
	sum1=0.0; for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ sum1+=ch[i][j]; } }
  c_0=sum1/nd/nd;
	//printf("c_0= %f \n", c_0);

//******** ���d�� K (���ϒl�ƕϓ���) ******************************************************
	K0=K1*c_0+K2*(1.0-c_0);

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			dKh[i][j]=(K1-K2)*(ch[i][j]-c_0);
		}
	}

//**** �d�ׂ̔����i���Łj�ʂ̃t�|���G�ϊ� s0h ---> s0qh1 ********************************
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

//***** �����v�Z *******************************************************************************************
	for(ief=0;ief<=loopief;ief++){

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			Vh2[i][j]=Vh[i][j];//�⏕�z��ɃR�s�[
		}
	}

//*** �d�ʌ��z��̌v�Z(�����ɒ���) **********
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				ip=i+1; im=i-1; jp=j+1; jm=j-1;
				if(i==ndm){ip=0;} 	if(i==0){im=ndm;}
				if(j==ndm){jp=0;}   if(j==0){jm=ndm;}

				J1[i][j]=0.5*(Vh[ip][j]-Vh[im][j]);
				J2[i][j]=0.5*(Vh[i][jp]-Vh[i][jm]);
			}
		}

//**** �d�ʌ��z*���d���̍��̃t�|���G�ϊ��idKh[][]*J1[][] ---> a1_qrh1�j ********************************
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
	
//**** �d�ʌ��z*���d���̍��̃t�|���G�ϊ� (dKh[][]*J2[][] ---> a2_qrh2�j ********************************
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
	
//***** ���̃X�e�b�v�̓d�ʏ�̌v�Z *************************************
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
				Vh[i][j]=w0*Vh[i][j]+(1.0-w0)*Vh2[i][j];//�d�ݕt������
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

	}//ief��loop


//*** �d�ʌ��z��̌v�Z *****************************************************
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				Vh[i][j]=Vh[i][j]*V2+Vav;//V�̎�����߂��A���ϓd�ʂ�Vav�ɃV�t�g
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
	printf("Vx1, Vx2, ��Vx= %f  %f  %f \n", Vx1, Vx2, Vx2-Vx1);
	printf("V1, V2, ��V= %f  %f  %f \n", V1, V2, V2-V1);

	Q00=Q01*(V2-V1)/(Vx2-Vx1);
	printf("Q01, Q00= %e  %e \n", Q01, Q00);

//*** ���ʂ̕\�� *****************************************************
	graph_J();//�d����
	//graph_V();//�d�ʌ��z��


//*** [���ς̓��d���̌v�Z] ***********************************************
	sum1=0.0; for(j=0;j<=nd2m;j++){ sum1+=J1[nd4][j]; }
	Km=fabs( sum1/nd2*(nd2-3)/(Vx2-Vx1) );
	//Km=fabs( sum1/nd2*(nd2-3)/(V2-V1) );
  //���ς̓��d���̌v�Z(3�������Ă���͍̂��E�̂P�u���b�N�ÂƁA���u���b�N�Q��)
	printf("Km= %e \n", Km);

//************************************************************
	if(keypress()){return 0;}//�L�[�҂����
	//time1=time1+1.;  if(time1<time1max){goto start;}//�ő�J�E���g���ɓ��B�������ǂ����̔��f

	end:;
  return 0;
}
			
//************[�����d�ʏ�]*****************************************
void shokiha_V()
{
	int i, j;
	int ii, jj;
  srand(time(NULL)); // ����������

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

//************[�d�ׂ̔����Ə���]*****************************************
void shokiha_S()
{
	int i, j;
	int ii, jj;
  srand(time(NULL)); // ����������

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

//*******[�d�ʏ�̕`��]**************************************************
void graph_V()
{
	int i, j, ii, jj;
	int col_R, col_G, col_B;
	double col, fact1;
	double c, x, xmax, xmin, y, ymax, ymin, rad0;
	int ixmin=0, iymin=0, igx, igy, irad0;
	int ixmax=INX;
	int iymax=INY;

	gcls(); //��ʃN���A
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

//*******[�g�D�̕`��]**************************************************
void graph_c()
{
	int i, j, ii, jj;
	double col;
	double c, x, xmax, xmin, y, ymax, ymin, rad0;
	int ixmin=0, iymin=0, igx, igy, irad0;
	int ixmax=INX;
	int iymax=INY;

	gcls(); //��ʃN���A
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

//******* �x�N�g����̕\�� ***************************************
void graph_J()
{
	int i, j, k, ii, jj, kk;//����
	double col;//�F
	double x, xmax, xmin, y, ymax, ymin, rad0, dia0;//�K�i�����W�n�̐ݒ�
	int ixmin=0, iymin=0, igx, igy, irad0;//�X�N���[�����W�n�̐ݒ�
	double x1, y1, x2, y2, x3, y3, x4, y4, fact1, fact2, th0;//���֘A�̍��W�ݒ�
	int igx1, igy1, igx2, igy2, igx3, igy3, igx4, igy4;//���֘A�̃X�N���[�����W�̐ݒ�
	int ixmax=INX, iymax=INY;//�`��Window�͈�

  gcls(); //��ʃN���A
	xmin=0.; xmax=1.; ymin=0.; ymax=1.;//�`��̈�i�K�i������Ă���j

	printf("time %f\n",time1);//�v�Z�J�E���g���̕\��
	dia0=1./nd;  
	rad0=dia0/2.;               irad0=(ixmax-ixmin)/(xmax-xmin)*rad0+1;
	//�����u���b�N�̔����̒���	//�X�N���[�����W�n�ɕϊ��i+1�͐��������̐؎̂ĕ␳�j

//**** ���n�̏�𔖂����ĕ`�� ***************************
	for(i=0;i<=nd;i++){
		for(j=0;j<=nd;j++){
			x=rad0+dia0*i;  igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
			y=rad0+dia0*j;  igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
			//���W�v�Z			//�X�N���[�����W�n�ɕϊ�
			ii=i; jj=j;  if(i==nd){ii=0;}  if(j==nd){jj=0;}//�����I���E����
			//col=ch[ii][jj];
			col=1.0-(1.0-ch[ii][jj])*0.5;//��̐F��RGB�ɂĐݒ�i0.5�������Ĕ������Ă���j
			if(col>=1.){col=1.;}  if(col<=0.){col=0.;}
			gcolor((int)(255*col),(int)(255*col),(int)(255*col));//�F�ݒ�
			grect(igx-irad0,igy-irad0,igx+irad0,igy+irad0);//���h��l�p�`�̕`��
		}
	}

//**** ��������ɂĕ`�� ***************************
	//fact1=0.01;//���̒���
	fact1=0.1;//���̒���
	fact2=0.2;//����[�ʒu�̐ݒ�
	th0=0.5*3.14159/2.;//��̕����̊p�x

	for(i=0;i<=nd;i+=4){
		for(j=0;j<=nd;j+=4){
			x=rad0+dia0*i;  igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
			y=rad0+dia0*j;  igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
			//���W�v�Z			//�X�N���[�����W�n�ɕϊ�
			ii=i; jj=j; if(i==nd){ii=0;}  if(j==nd){jj=0;}//�����I���E����

//--- ���̊e���ʂ̍��W�v�Z ---
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
			gcolor(0,0,0);//�F�ݒ�
			gline(igx1, igy1, igx2, igy2);//���������i���̕��i�j
			gline(igx2, igy2, igx3, igy3);//���������i���̕��i�j
			gline(igx2, igy2, igx4, igy4);//���������i���̕��i�j
		}
	}
	swapbuffers();//��ʃX���b�v
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

//********** �P���������t�[���G�ϊ� **************************************
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

//************ �Q���������t�[���G�ϊ� ***********************************
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

//************[��f�|�^�̕ۑ�]************************************
void datsave()
{
	FILE		*stream;
	int 		i, j;

	stream = fopen("test.dat", "a");		//�ۑ��t�@�C������test.dat�Ƃ��Ă���B
	fprintf(stream, "%e\n", time1);		//���Ԃ̏�������
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			fprintf(stream, "%e  %e  %e  ", Vh[i][j], J1[i][j], J2[i][j]);		//��̏�������
		}
	}
	fprintf(stream, "\n");
	fclose(stream);
}
//*********** �g�D�`�ԏ��̓��� **************************
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
			ch[i][j]=(ch[i][j]-cmin)/(cmax-cmin);//�K�i��
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

//*********** �g�D�`�ԏ��̓��� **************************
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
			ch[i][j]=(ch[i][j]-cmin)/(cmax-cmin);//�K�i��
			//ch[i][j]=1.-(ch[i][j]-cmin)/(cmax-cmin);//�K�i��
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
	//		Kh[i][j]=K2+(K1-K2)*ch[i][j];//���d���}�g���b�N�X�̐ݒ�
	//	}
	//}

}


