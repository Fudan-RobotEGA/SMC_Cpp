#include"slidingmodec.h"
#include<cmath>

SMC YawSMC(20, 120, 0, 0.001, 21, 27, 25000, 0.8, 0.5);

void SMC::SMC_Tick(float angle_now,float angle_vel) //anlge为当前位置(°),ang_vel为角速度(rad/s)
{
	angle = angle_now;
    ang_vel = angle_vel;
	float e_qp;
	float qp = q/p;
	//读取参数
	error = angle - ref;
	ddref = (ref - refl) - dref;
	dref = (ref - refl);

	//误差下限处理
	if (fabs(error) < error_eps)
	{
		u = 0;
		return;
	}

	if (error < 0)
		e_qp = -pow(abs(error), qp);
	else
		e_qp = pow(abs(error), qp);
	// //快速终端滑模
	s = (ang_vel - dref) + C * e_qp;
	ds = -epsilon * Sat(s) - K * s;
	u = -J * (ddref + ds - C * qp * (ang_vel - dref) * e_qp / abs(error));

	if (abs(error) < 1)
	{
		s = C * error + (ang_vel - dref);//smc surface
		u = -J * (ddref - C * (ang_vel - dref) - epsilon * Sat(s) - K * s);
	}

	//控制量限幅
	if (u > u_max)
		u = u_max;
	if (u < -u_max)
		u = -u_max;


	//参数更新
	refl = ref;
}