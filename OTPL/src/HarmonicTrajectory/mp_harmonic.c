/**
* @file mp_harmonic.c
* @brif �ú�������г���켣��
*
* ����������ʼʱ��t0,��ֹʱ��t1,��ʼλ��q0����ֹλ��q1�Ͳ�ֵ�����ʱ��ָ�����dT
* ����г���켣��
* @author LiBing
* @date 2018/1/9
*/
#include<stdio.h>
#include<math.h>
#include<mp_harmonic.h>

#define MP_OK 0			/*�ٶȹ滮�ɹ�*/
#define PARA_ERR 1		/*��������Ƿ�*/


#define MIN_ERR  1.0E-6 /*��ֵ�Ƚ���������*/


#define DEBUG 1
int ret = 0;

//�켣�滮�����ӿ�
int mp_harmonic_traj(double *t0,double *t1,double *q0,double *q1,double *dT,mp_para *out)
{
	if ((ret = mp_para_initial(t0, t1,q0, q1, dT, out)) != MP_OK)
	{
		if (DEBUG)
			printf("��������Ƿ�\n");
		return ret;
	}



	if (DEBUG)printf("mp_para_initial ok\n");
	return MP_OK;
}


/*��ʼ���������*/
static int mp_para_initial(double *t0, double *t1, double *q0, double *q1, double *dT, mp_para *out)
{	
	
	if (*t0 < 0.0f || *t1 < 0.0f || *q0<0.0f ||*q1<0.0f || *dT<0.0f)
	{
		ret = PARA_ERR;
		return ret;
	}
	out->t0 = *t0;
	out->t1 = *t1;
	out->T = *t1 - *t0;
	out->q0 = *q0;
	out->q1 = *q1;
	out->h = *q1 - *q0;
	out->dT = *dT;
	return MP_OK;

}