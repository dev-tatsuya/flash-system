//����
//�g�D�ω�����

#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>
#include <iostream>
// #include "wingxa.h"

using namespace std;

#define DRND(x) ((double)(x)/RAND_MAX*rand())
#define ND 128 //������
#define INXY 400
#define PI 3.1415956535
#define RR 8.3144598
#define FF 9.64853415e+04                //[C/mol] �t�@���f�[�萔
#define VAVA 5300

int i, j;
int nd=ND, ndm=ND-1;        //�g�D�̕���
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
	double V, V_E, V_W, V_N, V_S; 	//�����u���b�N�ɂ�����V�𒆐S�ɁA���̏㉺���E
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

	K01=1.0;      //[S/m]=[1/(��m)]  ZrO2�\��
	K02=1.0e-18;  //[S/m]=[1/(��m)]  �^��̒l��K���ɐݒ�
	K1=K01/K02;  K2=K02/K02;   //��������

	al=3.0;            //�v�Z�̈�̂P��(��m)
	al=1.0e-06*al;     //�v�Z�̈�̂P��(m)
	b1=al/(double)nd;  //�����P�u���b�N�̃T�C�Y

	beta=1.5;
	V1=0.3;       //[V]���d��
	V2=0.01;      //[V]��d��

	time1=0.0;
	time1max=1.0e7;

	sum1=0.0;

	den=6.05e3;
	Cp=0.66e3*den;

	//****** ��̓ǂݍ��� ********************
	datin();
	datin2();
	shokiha_V();	//�����d�ʏ�̐ݒ�

	//*** �J��Ԃ��v�Z�X�^�|�g *******************************************

	cout << "1e^7 times start" << endl;

	for(time1=0.;time1<=time1max;time1+=1.){

		for(i=2;i<=ndm-2;i++){//���[�ƉE�[���Ȃ�
			for(j=0;j<=ndm;j++){
				ip=i+1; im=i-1; jp=j+1; jm=j-1;
				if(i==ndm){ip=ndm-1;} 	if(i==0){im=1;}
				if(j==ndm){jp=0;}       if(j==0){jm=ndm;}

				V=Vh[i][j]; V_E=Vh[ip][j];  V_W=Vh[im][j];  V_N=Vh[i][jp];  V_S=Vh[i][jm];

				K=Kh[i][j];
				K_e=(Kh[ip][j]+K)/2.0;  K_w=(Kh[im][j]+K)/2.0;
				K_n=(Kh[i][jp]+K)/2.0;  K_s=(Kh[i][jm]+K)/2.0;

				Vh[i][j]=beta*(K_e*V_E+K_w*V_W+K_n*V_N+K_s*V_S)/(K_e+K_w+K_n+K_s)+(1.0-beta)*V;
			}
		}

		for(j=0;j<=ndm;j++){ Vh[ndm][j]=Vh[ndm-1][j]=V1;  Vh[0][j]=Vh[1][j]=V2; }

		sum2=0.0;
		for(i=0;i<=ndm;i++){for(j=0;j<=ndm;j++){sum2+=Vh[i][j];}}

		if(fabs(sum1-sum2)<=1.0e-8){break;} 
		else{sum1=sum2;}

	}//time1

	cout << "1e^7 times end" << endl;

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

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			cVa[i][j]=G_T_V(T[i][j]);
		}
		//printf("e");
	}

	// gwinsize(INXY,INXY); ginit(2); gsetorg(0,0);      //�`��Window�\��
	//graph_c();
	//graph_V();
	//graph_delT();
	//graph_yVa();

	// datsave_c();
	// datsave_V();
	// datsave_E();
	// datsave_W();
	// datsave_T();
	// datsave_delT();
	// datsave_yVa();
	// if(keypress()){return 0;}
	cout << "main end" << endl;
	return 0;
}//main


//*********** �g�D�`�ԏ��̓��� **************************
void datin(){
	cout << "datin" << endl;
	
	FILE   *datin0;

 	double ceM, ceP;
	double cmax, cmin; 
	double time2;

    ceM=0.0001;  ceP=0.9999;

	datin0 = fopen("dc/no_change/data/test.dat", "r");

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

	if (time2 != XX) {goto start;}
	fclose(datin0);

	cout << "datin end" << endl;
}


//*********** �g�D�`�ԏ��̓��� **************************
void datin2(){
	cout << "datin2" << endl;

	FILE		*datin1;

	datin1 = fopen("dc/no_change/data/temp3.dat", "r");
	//fscanf(datin0, "%lf", &time1);

	for(i=0;i<vava;i++){
		fscanf(datin1, "%lf\n", &VVa[i]);
	}

	fclose(datin1);

	cout << "datin2 end" << endl;
}

//************[�����d�ʏ�]*************************
void shokiha_V(){
	cout << "shokiha_V" << endl;

	double rnd0; 
  	srand(time(NULL)); // ����������

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			Vh[i][j]=V2;
		}
	}

	for(j=0;j<=ndm;j++){
		Vh[ndm][j]=Vh[ndm-1][j]=V1;
		Vh[0][j]=Vh[1][j]=V2;
	}

	cout << "shokiha_V end" << endl;
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

	stream = fopen("test_c.bin", "ab");	//�������ސ�̃t�@�C����ǋL�����ŃI�[�v��
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

	stream = fopen("test_V.bin", "ab");	//�������ސ�̃t�@�C����ǋL�����ŃI�[�v��
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

//*******[�g�D�̕`��]**************************************************
// void graph_c()
// {
// 	int i, j, ii, jj, A,B;
// 	double col,col_R,col_G,col_B;
// 	double c, x, xmax, xmin, y, ymax, ymin, rad0;
// 	int ixmin=0, iymin=0, igx, igy, irad0;
// 	int ixmax=INXY;
// 	int iymax=INXY;

// 	gcls(); //��ʃN���A
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
// 	swapbuffers();     //��ʃX���b�v
// }

//*******[�g�D�̕`��]**************************************************
// void graph_V()
// {
// 	int i, j, ii, jj;
// 	double col,col_R,col_G,col_B;
// 	double c, x, xmax, xmin, y, ymax, ymin, rad0;
// 	int ixmin=0, iymin=0, igx, igy, irad0;
// 	int ixmax=INXY;
// 	int iymax=INXY;
// 	double A,B,C,D,E,F;

// 	gcls(); //��ʃN���A

// 	A=178.0; B=102.0; C=255.0;   //��
// 	//A=255.0; B=255.0; C=255.0;   //��
// 	D=255.0; E=102.0; F=178.0;   //��

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
// 	swapbuffers();     //��ʃX���b�v
// }


//*******[�g�D�̕`��]**************************************************
// void graph_T()
// {
// 	int i, j, ii, jj;
// 	double col,col_R,col_G,col_B;
// 	double c, x, xmax, xmin, y, ymax, ymin, rad0;
// 	int ixmin=0, iymin=0, igx, igy, irad0;
// 	int ixmax=INXY;
// 	int iymax=INXY;
// 	double A,B,C,D,E,F;

// 	gcls(); //��ʃN���A

// 	D=255.0; E=0.0; F=0.0;   //��
// 	A=0.0; B=0.0; C=255.0;   //��
// 	//A=255.0; B=255.0; C=255.0;   //��

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
// 	swapbuffers();     //��ʃX���b�v
// }

//*******[�g�D�̕`��]**************************************************
// void graph_delT()
// {
// 	int i, j, ii, jj;
// 	double col,col_R,col_G,col_B;
// 	double c, x, xmax, xmin, y, ymax, ymin, rad0;
// 	int ixmin=0, iymin=0, igx, igy, irad0;
// 	int ixmax=INXY;
// 	int iymax=INXY;
// 	double A,B,C,D,E,F;

// 	gcls(); //��ʃN���A

// 	D=255.0; E=0.0; F=0.0;   //��
// 	//A=193.0; B=224.0; C=255.0;   //��
// 	A=255.0; B=255.0; C=255.0;   //��

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
// 	swapbuffers();     //��ʃX���b�v
// }

//*******[�g�D�̕`��]**************************************************
// void graph_yVa()
// {
// 	int i, j, ii, jj;
// 	double col,col_R,col_G,col_B;
// 	double c, x, xmax, xmin, y, ymax, ymin, rad0;
// 	int ixmin=0, iymin=0, igx, igy, irad0;
// 	int ixmax=INXY;
// 	int iymax=INXY;
// 	double A,B,C,D,E,F;

// 	gcls(); //��ʃN���A

// 	A=220.0; B=220.0; C=220.0;   //��
// 	//A=255.0; B=255.0; C=255.0;   //��
// 	// D=242.0; E=135.0; F=18.0;   //��
// 	D=255.0; E=0.0; F=127.0;   //��

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
// 	swapbuffers();     //��ʃX���b�v
// }
