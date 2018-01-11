/*harmonic_rajectory*/
#ifndef MP_HAMONIC_H_
#define MP_HAMONIC_H_

#ifdef __cplusplus
extern "C" {
#endif // 

/**
@brief �켣�滮�ṹ��

���ڴ洢��������ͺ���������
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
	double t0;	
	double t1;	
	double T;	
	double q0;	
	double q1;
	double dT;
	double h;
	double vs;
	double ve;
	double acc;
	double jerk;
	double snap;
}mp_para;

int mp_harmonic_traj(double *t0, double *t1, double *q0, double *q1, double *dT, mp_para *out);

#ifdef __cplusplus
}
#endif // __cplusplus

#endif