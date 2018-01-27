/*mp_cycloidal.h*/

#ifndef MP_CYCLOIDAL_H_
#define MP_CYCLOIDA_H_

#ifdef __cplusplus  
extern "C" {
#endif // 


	/**
	@brief ���߹켣�滮�ӿں�����

	@para t0,double����ʼʱ�̡�
	@para t1,double����ֹʱ�̡�
	@para q0,double����ʼλ�ơ�
	@para q1,double���յ�λ�ơ�
	@para dT,double��ʱ��ָ�����
	@para *data[6],double��������ݵ�ָ�����飬data[0]ָ��ʱ�䣬data[1]ָ��λ�ƣ�data[2]ָ���ٶ�,...��
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
	int mp_cycloidal_traj(double t0, double t1, double q0, double q1, double dT, double *data[6], int *dataNum);

#ifdef __cplusplus
}
#endif // __cplusplus

#endif
