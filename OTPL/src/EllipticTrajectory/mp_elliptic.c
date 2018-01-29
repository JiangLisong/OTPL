/**
* @file mp_elliptic.c
* @brif ���ļ��ĺ���ʵ����Բ�켣�滮��
*
* ����������ʼʱ��t0,��ֹʱ��t1,��ʼλ��q0����ֹλ��q1�Ͳ�ֵ�����ʱ��ָ�����dT
* ������Բ�켣��
* @author:LiBing
* @date  :2018/1/27
*/
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<mp_elliptic.h>

#define MP_OK			0	/*�ٶȹ滮�ɹ�*/
#define PARA_ERR		1	/*��������Ƿ�*/
#define MALLOC_ERR		2   /*�ڴ����ʧ��*/

#define MIN_ERR  1.0E-6 /*��ֵ�Ƚ���������*/
#define TIME_ERR  1.0E-6 /*ʱ�����������С*/
#define M_PI       3.14159265358979323846264338328      /* pi */   

#define INTEGER(x) ((int)(x))    //ȡ��

#define DEBUG 1

/**
@brief �켣�滮�ṹ��

���ڴ洢��������ͺ���������
@para n ��Բ�ĳ��������֮��
@para t0 ��ʼʱ��
@para t1 �յ�ʱ��
@para T  �˶�����ʱ��
@para q0 ��ʼλ��
@para q1 �յ�λ��
@para dT ��ֵ�����ʱ��ָ
@para h  ���λ��
@para vs ��ʼ�ٶ�
@para ve �յ��ٶ�
@para acc ���ٶ�
@para jerk �Ӽ��ٶ�
@para snap �ӼӼ��ٶ�
�˶��滮�����ṹ��
*/
typedef struct {
	double n;		/**<ratio between major and minor ellipse axes */
	double t0;		/**<initial time*/
	double t1;		/**<final time*/
	double T;		/**<duration time*/
	double q0;		/**<initial position*/
	double q1;		/**<final position*/
	double dT;		/**<interpolation cycle*/
	double Tn[3];	/**<duration time  of each continue segments */
	double h;		///<relative displacement
	double vs;		///<initial velocity
	double ve;		///<final velocity
	double acc;     ///<acceleration
	double jerk;	///<jerk
	double snap;	///<snap
	int dataN;		///<numbers of the interpolation data
}mp_para;


int ret = 0;

/*��ʼ���������*/
static int mp_para_initial(double n,double t0, double t1, double q0, double q1, double dT, mp_para *out)
{
	memset(out, 0, sizeof(mp_para));
	if (t0 < 0.0f || t1 < 0.0f || q0<0.0f || q1<0.0f || dT<0.0f || n<1.0E-6)
	{
		ret = PARA_ERR;
		return ret;
	}
	out->n = n;
	out->t0 = t0;
	out->t1 = t1;
	out->T = t1 - t0;
	out->q0 = q0;
	out->q1 = q1;
	out->h = q1 - q0;
	out->dT = dT;
	return MP_OK;

}


/**
@brief ����ʱ��ָ����ָ����ݵ�

����ʱ�䡢λ�ơ��ٶȺͼ��ٶȵ�
@para *in ָ���˶��滮�ṹ�壨mp_para����ָ��
*/
static int mp_data_generate(mp_para *in, double *data[6], int *dataNum)
{
	int T_N = 0;//Total number of interpolation cycles of motion
	int t_N[3] = { 0 };//number of interpolation cycles of each section
	int ret = 0;
	int i = 0;
	double temp1 = 0;
	double temp2 = 0;
	//Temporary data
	double timeTemp = in->t0;
	double DispTemp = in->q0;
	double vecTemp = in->vs;
	double accTemp = in->acc;
	double jerkTemp = in->jerk;

	double tRemain = 0.0;/**<��ǰ�μ�ȥ���������ں�ʣ���ʱ��*/
	double tCompen = 0.0;/**<��Ҫ����һ�����ڲ�����ʱ��*/
	double timeTotal = 0.0;/**<�˶�����ʱ��*/
	double TnSplit[3] = { 0 };/**<���ڴ洢�ָ����ε�ʱ��,ֻ��һ��*/

	int cur_n = 0;
	int next_n = cur_n + 1;
	double alpha = (in->n*in->n - 1) / (in->n*in->n);
	double w = M_PI / in->T;


	for (int i = 0; i < 1; i++)
	{
		TnSplit[i] = in->Tn[i];
		if (in->Tn[i] > TIME_ERR)
		{
			timeTotal = timeTotal + in->Tn[i];
			t_N[i] = INTEGER(in->Tn[i] / in->dT);
		}
	}
	T_N = INTEGER(timeTotal / in->dT);
	if (fabs(timeTotal - T_N*in->dT) >= TIME_ERR)
	{
		T_N = T_N + 1;
	}
	T_N = T_N + 1;
	//�����ڴ�
	for (i = 0; i<6; i++)
	{
		if ((data[i] = (double *)malloc(T_N * sizeof(double))) == NULL)
		{
			ret = MALLOC_ERR;
			printf("�ڴ������� error:%d\n", ret);
			return ret;
		}
		
		else 
		{ 
			memset(data[i], 0, (size_t)T_N); ret = MP_OK;
		}
	}

	data[0][0] = in->t0;
	data[1][0] = in->q0;
	data[2][0] = in->vs;
	data[3][0] = in->acc;
	data[4][0] = in->jerk;
	data[5][0] = in->snap;


	//��������
	if (in->Tn[0]>TIME_ERR)
	{
		//�ֲ�ʱ������
		double tau1 = 0;
		if (t_N[0] >= 1)
		{
			//ʱ����ڻ����һ������ʱ�����м���
			for (int i = cur_n; i< cur_n + t_N[0] + 1; i++)
			{
				tau1 = i*in->dT;
				temp1 = sin(w*tau1);
				temp2 = cos(w*tau1);
				//ʱ�䡢λ�á��ٶȺͼ��ٶ�
				data[0][i] = timeTemp + tau1;
				data[1][i] = DispTemp + in->h/2.0 *(1.0f - temp2/sqrt(1.0f-alpha*temp1*temp1));
				data[2][i] = vecTemp + w*in->h / 2.0f*temp1 / (in->n*in->n) /
					sqrt(pow(1.0f - alpha*temp1*temp1, 3));
				data[3][i] = accTemp + w*w*in->h / 2.0f  * temp2*(1 + 2 * alpha*temp1*temp1)
					/ (in->n*in->n) / sqrt(pow(1.0f - alpha*temp1*temp1, 5));
				data[4][i] = jerkTemp - w*w*w*in->h / 2.0 * temp1*(1 + 2 * alpha*temp1*temp1) / (in->n*in->n) /
					sqrt(pow(1.0f - alpha*temp1*temp1, 5)) + w*w*w*in->h / 2.0f / (in->n*in->n)*temp2*(2.0f*alpha*temp1*temp2
						*sqrt(pow(1.0f - alpha*temp1*temp1, 5)) + 5.0f*alpha*temp1*temp2 * (1 + 2.0f*alpha*temp1*temp1)
						*sqrt(pow(1 - alpha*temp1*temp1, 3))) / pow(1.0f - alpha*temp1*temp1, 5);
			}
			next_n = cur_n + t_N[0] + 1;
		}

		//�����л��㣨�öεĽ����㣩��������
		tau1 = in->Tn[0];
		timeTemp = in->Tn[0];
		temp1 = sin(w*tau1);
		temp2 = cos(w*tau1);
		DispTemp = DispTemp + in->h / 2.0 *(1.0f - cos(w*tau1) / sqrt(1.0f - alpha*temp1*temp1));
		vecTemp = vecTemp + M_PI*in->h / (2.0f*in->T)*temp1 / (in->n*in->n) /
			sqrt(pow(1.0f - alpha*temp1*temp1, 3));
		accTemp = accTemp + M_PI * M_PI*in->h / (2.0 *in->T*in->T) * temp2*(1 + 2 * alpha*temp1*temp1)
			/ (in->n*in->n) / sqrt(pow(1.0f - alpha*temp1*temp1, 5));
		jerkTemp = jerkTemp - w*w*w*in->h / 2.0 * temp1*(1 + 2 * alpha*temp1*temp1) / (in->n*in->n) /
			sqrt(pow(1.0f - alpha*temp1*temp1, 5)) + w*w*w*in->h / 2.0f / (in->n*in->n)*temp2*(2.0f*alpha*temp1*temp2
				*sqrt(pow(1.0f - alpha*temp1*temp1, 5)) + 5.0f*alpha*temp1*temp2 * (1 + 2.0f*alpha*temp1*temp1)
				*sqrt(pow(1 - alpha*temp1*temp1, 3))) / pow(1.0f - alpha*temp1*temp1, 5);

		tRemain = TnSplit[0] - t_N[0] * in->dT;
		if (tRemain > TIME_ERR)
		{
			tCompen = in->dT - tRemain;
			//ʣ��ʱ�䲻��һ�����ڣ���һ��������
			//��ĩ�ٶȽϴ����һ�����ڻ���ڽϴ���ٶ�ͻ��
			data[0][next_n] = timeTemp + tCompen;
			data[1][next_n] = DispTemp;
			data[2][next_n] = (vecTemp - data[2][next_n - 1]) / in->dT;
			data[3][next_n] = (accTemp - data[3][next_n - 1]) / in->dT;
			data[4][next_n] = (jerkTemp - data[3][next_n - 1]) / in->dT;
			next_n = next_n + 1;
		}
		cur_n = next_n;
	}
	in->dataN = cur_n;
	*dataNum = cur_n;
	return MP_OK;
}


/**
@brief ���߹켣�滮�ӿں�����

@para n ��Բ�ĳ��������֮��
@para t0,double����ʼʱ�̡�
@para t1,double����ֹʱ�̡�
@para q0,double����ʼλ�ơ�
@para q1,double���յ�λ�ơ�
@para dT,double��ʱ��ָ�����
@para *data[6],double��������ݵ�ָ�����飬data[0]ָ��ʱ�䣬data[1]ָ��λ�ƣ�data[2]ָ���ٶ�,...������
��Բ�켣ֻ���λ��ǰ���׵����������ٶȣ����ڶ��׵�����ֵ��Ч��
@para *data double��ָ�룬ָ�����ݴ�С����data[i]��Ԫ�ظ�����
@return ����ֵָʾ������ִ�������
@retval 0,����ִ�гɹ���
@retval 1,��������Ƿ���
@retval 2,��̬�ڴ�������
@waring ���øú���������dataָ��ָ����ڴ��ͷţ�ִ�����´��룺
for (int i = 0; i < 6; i++)
{
free(data[i]);
}
*/
int mp_elliptic_traj(double n,double t0, double t1, double q0, double q1, double dT, double *data[6], int *dataNum)
{
	mp_para mp_st;
	if ((ret = mp_para_initial(n,t0, t1, q0, q1, dT, &mp_st)) != MP_OK)
	{
#ifdef DEBUG
		printf("��������Ƿ�\n");
#endif // DEBUG		
		return ret;
	}
	mp_st.Tn[0] = mp_st.T;
	//�������ݵ�
	if ((ret = mp_data_generate(&mp_st, data, dataNum)) != MP_OK)
	{
#ifdef DEBUG
		printf("�������ݵ�ʧ��\n");
#endif // DEBUG

		return ret;
	}

#ifdef DEBUG
	printf("mp_para_initial ok\n");
#endif // DEBUG
	return MP_OK;
}


