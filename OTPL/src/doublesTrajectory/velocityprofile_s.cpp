#include <iostream>
#include <stdio.h>
#include <math.h>

#include "velocityprofile_s.h"

#define min(x,y) (x < y ? x : y)

VelocityProfile_S::VelocityProfile_S():
		  Ta(0), Tv(0), Td(0), Tj1(0), Tj2(0), Tj(0), T(0),
		  maxvel(0), minvel(0), 
		  maxacc(0), minacc(0),
		  maxj(0), minj(0),
          startpos(0), endpos(0),
          startvel(0), endvel(0),
          s(0)
{}
		// constructs motion profile class with <maxvel> as parameter of the
		// trajectory.

void VelocityProfile_S::SetProfile(double pos1,double pos2, double vel1, double vel2) {
	startpos = pos1;
	endpos   = pos2;
	startvel = vel1;
	endvel   = vel2;

	if(fabs(endpos-startpos) < 1e-10)
	{
		printf("displacement should not be 0, quit!\n");
		return;
	}

	if((endpos-startpos) > 0.0)
	{
		s = 1;
	} else {
		s = -1;
	}

	startpos = s * pos1;
	endpos   = s * pos2;
	startvel = s * vel1;
	endvel   = s * vel2;

	maxvel = (s+1)/2*maxvel + (s-1)/2*minvel;
	minvel = (s+1)/2*minvel + (s-1)/2*maxvel;
	maxacc = (s+1)/2*maxacc + (s-1)/2*minacc;
	minacc = (s+1)/2*minacc + (s-1)/2*maxacc;
	maxj   = (s+1)/2*maxj + (s-1)/2*minj;
	minj   = (s+1)/2*minj + (s-1)/2*maxj;

	double T1 = sqrt(fabs(endvel - startvel) / maxj);
	double T2 = maxacc / maxj;
	Tj = min(T1, T2);
	if(T1 < T2)
	{
		if((endpos - startpos) < (Tj * (startvel + endvel)))
		{
			printf("displacement is too small, not satisfied the requirement of velocity, quit!\n");
			return;
		}
	} else {
		if((endpos - startpos) < 0.5 * (startvel + endvel) * (Tj + fabs(endvel - startvel) / maxacc))
		{
			printf("displacement is too small, not satisfied the requirement of velocity, quit!\n");
			return;			
		}
	}
	if(((maxvel - startvel)*maxj) < (maxacc*maxacc))
	{
		printf("maximum acceleration is not reached!\n");
		Tj1 = sqrt((maxvel - startvel)/maxj);
		Ta = 2*Tj1;
		a_lima = maxj * Tj1;		
	} else {
		printf("maximum acceleration is reached!\n");
		Tj1 = maxacc/maxj;
		Ta = Tj1 + (maxvel - startvel)/maxacc;
		a_lima = maxacc;
	}

	if(((maxvel - endvel)*maxj) < (maxacc*maxacc))
	{
		printf("minimum acceleration is not reached!\n");
		Tj2 = sqrt((maxvel - endvel)/maxj);
		Td = 2*Tj2;
		a_limd = minj * Tj2;
	} else {
		printf("minimum acceleration is reached!\n");
		Tj2 = maxacc/maxj;
		Td = Tj2 + (maxvel - endvel)/maxacc;
		a_limd = minacc;
	}

	Tv = (endpos - startpos)/maxvel - Ta/2*(1+startvel/maxvel) - Td/2*(1+endvel/maxvel);


	if(Tv > 0.0)
	{
		printf("maxvel is reached!\n");	

		v_lim = maxvel;
		T = Ta + Tv + Td;
		printf("Tj1:%f, Ta:%f, Tv:%f, Tj2:%f, Td:%f, T:%f, v_lim:%f, a_lima:%f, a_limd:%f\n"
		,Tj1, Ta, Tv, Tj2, Td, T
		,v_lim, a_lima, a_limd );		

		return;
	} else {
		printf("maxvel is not reached!\n");
		Tv = 0.0;
		Tj = maxacc / maxj;
		Tj1 = Tj;
		Tj2 = Tj;

		double delta = pow(maxacc, 4) / pow(maxj, 2) + 2 * ( pow(startvel, 2) + pow(endvel, 2) ) + maxacc * (4 * (endpos - startpos) - 2 * maxacc / maxj * (startvel + endvel));
		Ta = (pow(maxacc, 2) / maxj - 2 * startvel + sqrt(delta)) / (2 * maxacc);
		Td = (pow(maxacc, 2) / maxj - 2 * endvel + sqrt(delta)) / (2 * maxacc);

		if(Ta > 2 * Tj && Td > 2 * Tj)
		{
			printf("both of maximum acceleration and minimum acceleration are reached!\n");
			T = Ta + Tv + Td;
			a_lima = maxacc;
			a_limd = minacc;
			v_lim = startvel + (Ta - Tj1) * a_lima;
			
			printf("Tj1:%f, Ta:%f, Tv:%f, Tj2:%f, Td:%f, T:%f, v_lim:%f, a_lima:%f, a_limd:%f\n"
			,Tj1, Ta, Tv, Tj2, Td, T
			,v_lim, a_lima, a_limd );
			return;
		} else {
			printf("at least one segment is not reached!\n");

			double gamma = 0.99;

			while(Ta < 2 * Tj || Td < 2 * Tj)
			{
				if(Ta > 0 && Td > 0)
				{
					printf("please ddeccelerate the acceleration!\n");
					maxacc = gamma * maxacc;
					Tj = maxacc / maxj;
					Tj1 = Tj;
					Tj2 = Tj;
					delta = pow(maxacc, 4) / pow(maxj, 2) + 2 * ( pow(startvel, 2) + pow(endvel, 2) ) + maxacc * (4 * (endpos - startpos) - 2 * maxacc / maxj * (startvel + endvel));
					Ta = (pow(maxacc, 2) / maxj - 2 * startvel + sqrt(delta)) / (2 * maxacc);
					Td = (pow(maxacc, 2) / maxj - 2 * endvel + sqrt(delta)) / (2 * maxacc);					
				} else {
					if(Ta <= 0)
					{
						printf("maximum acceleration is not reached!\n");
						Ta = 0.0;
						Tj1 = 0.0;
						Td = 2 * (endpos - startpos) / (endvel + startvel);
						Tj2 = (maxj * (endpos - startpos) - sqrt(maxj * (maxj * pow((endpos-startpos), 2) + pow((endvel + startvel), 2) * (endvel - startvel)))) / (maxj * (endvel + startvel));								
					} else if(Td <= 0)
					{
						printf("minimum acceleration is not reached!\n");
						Td = 0.0;
						Tj2 = 0.0;
						Ta = 2 * (endpos - startpos) / (endvel + startvel);
						Tj1 = (maxj * (endpos - startpos) - sqrt(maxj * (maxj * pow((endpos-startpos), 2) - pow((endvel + startvel), 2) * (endvel - startvel)))) / (maxj * (endvel + startvel));

					}
					T = Ta + Tv + Td;
					a_lima = maxj * Tj1;
					a_limd = minj * Tj2;
					v_lim = startvel + (Ta - Tj1) * a_lima;
					
					printf("Tj1:%f, Ta:%f, Tv:%f, Tj2:%f, Td:%f, T:%f, v_lim:%f, a_lima:%f, a_limd:%f\n"
					,Tj1, Ta, Tv, Tj2, Td, T
					,v_lim, a_lima, a_limd );
					return;
				}
			}
			T = Ta + Tv + Td;
			a_lima = maxj * Tj1;
			a_limd = minj * Tj2;
			v_lim = startvel + (Ta - Tj1) * a_lima;
			
			printf("Tj1:%f, Ta:%f, Tv:%f, Tj2:%f, Td:%f, T:%f, v_lim:%f, a_lima:%f, a_limd:%f\n"
			,Tj1, Ta, Tv, Tj2, Td, T
			,v_lim, a_lima, a_limd );
			return;
		}
	}
}

void VelocityProfile_S::ExtendDuration(
	double pos1,double pos2, double vel1, double vel2, double newduration) {
	// duration should be longer than originally planned duration
    // Fastest :
	SetProfile(pos1, pos2, vel1, vel2);
    // Must be Slower  :
	double factor = T/newduration;
    if (factor > 1)
        return; // do not exceed max
    maxvel = factor * maxvel;
    minvel = factor * minvel;
    maxacc = factor * factor * maxacc;
    minacc = factor * factor * minacc;
    maxj = factor * factor * factor * maxj;
    minj = factor * factor * factor * minj;
    startvel = factor * startvel;
    endvel = factor * endvel;
}

void VelocityProfile_S::SetMax(double _maxvel, double _maxacc, double _maxjerk)
{
	if(fabs(_maxvel) < 1e-10)
	{
		printf("parameter is wrong, maximum velocity should not be 0!\n");
		return;
	}

	if(fabs(_maxacc) < 1e-10)
	{
		printf("parameter is wrong, maximum acceleration should not be 0!\n");
		return;
	}

	if(fabs(_maxjerk) < 1e-10)
	{
		printf("parameter is wrong, maximum Jerk should not be 0!\n");
		return;
	}	

    maxvel = _maxvel; 
    minvel = - 1.0 * _maxvel;
       
    maxacc = _maxacc; 
    minacc = - 1.0 * _maxacc;  

    maxj = _maxjerk;
    minj = - 1.0 * _maxjerk;
}

double VelocityProfile_S::Pos(double t) const {
	double m_pos = 0.0;;	
	if( t >= 0 && t < Tj1 )
	{
		m_pos = startpos + startvel * t + maxj * pow(t, 3) / 6.0;	
	} else if( t >= Tj1 && t < (Ta - Tj1) )
	{
		m_pos = startpos + startvel * t + a_lima/6.0 * ( 3 * pow(t, 2) - 3 * Tj1 * t + pow(Tj1, 2) );
	} else if( t >= (Ta - Tj1) && t < Ta )
	{
		m_pos = startpos + (v_lim + startvel) * Ta / 2.0 - v_lim * (Ta - t) - minj * pow((Ta - t), 3) / 6.0;
	} else if( t >= Ta && t < (Ta + Tv))
	{
		m_pos = startpos + (v_lim + startvel) * Ta / 2.0 + v_lim * (t - Ta);
	} else if(t >= (T - Td) && t < (T - Td + Tj2))
	{
		m_pos = endpos - ( v_lim + endvel ) * Td / 2.0 + v_lim * (t - T + Td) - maxj * pow((t - T + Td), 3) / 6.0;
	} else if(t >= (T - Td + Tj2 ) && t < (T - Tj2))
	{
		m_pos = endpos - (v_lim + endvel) * Td / 2.0 + v_lim * (t - T + Td) + a_limd / 6.0 * (3 * pow((t - T + Td), 2) - 3 * Tj2 * (t - T + Td) + pow(Tj2, 2));
	} else
	{
		m_pos = endpos - endvel * (T - t) - maxj * pow((T - t), 3) / 6.0;
	}

	m_pos = m_pos * s;
	printf("%f, ", m_pos);
	return m_pos;
}
double VelocityProfile_S::Vel(double t) const {
	double m_vel = 0.0;
	if( t >= 0 && t < Tj1 )
	{
		m_vel = startvel + maxj * pow(t, 2) / 2.0;
	} else if( t >= Tj1 && t < (Ta - Tj1) )
	{
		m_vel = startvel + a_lima *(t - Tj1/2.0);
	} else if( t >= (Ta - Tj1) && t < Ta )
	{
		m_vel = v_lim + minj * pow((Ta - t), 2) / 2.0;
	} else if( t >= Ta && t < (Ta + Tv))
	{
		m_vel = v_lim;
	} else if(t >= (T - Td) && t < (T - Td + Tj2))
	{
		m_vel = v_lim - maxj * pow((t - T + Td), 2) / 2.0;
	} else if(t >= (T - Td + Tj2 ) && t < (T - Tj2))
	{
		m_vel = v_lim + a_limd * ( t - T + Td - Tj2 / 2.0 );
	} else
	{
		m_vel = endvel + maxj * pow((T - t), 2) / 2.0;
	}

	m_vel = m_vel * s;
	printf("%f, ", m_vel);
	return m_vel;
}

double VelocityProfile_S::Acc(double t) const {
	double m_acc = 0.0;
	if( t >= 0 && t < Tj1 )
	{
		m_acc = maxj * t;
	} else if( t >= Tj1 && t < (Ta - Tj1) )
	{
		m_acc = a_lima;
	} else if( t >= (Ta - Tj1) && t < Ta )
	{
		m_acc = - minj * (Ta - t);
	} else if( t >= Ta && t < (Ta + Tv))
	{
		m_acc = 0.0;
	} else if(t >= (T - Td) && t < (T - Td + Tj2))
	{
		m_acc = - maxj * ( t - T + Td );
	} else if(t >= (T - Td + Tj2 ) && t < (T - Tj2))
	{
		m_acc = a_limd;
	} else
	{
		m_acc = - maxj * (T - t);
	}

	m_acc = m_acc * s;
	printf("%f, ", m_acc);
	return m_acc;
}
double VelocityProfile_S::Jerk(double t) const {
	double m_jerk = 0.0;
	if( t >= 0 && t < Tj1 )
	{
		m_jerk = maxj;
	} else if( t >= Tj1 && t < (Ta - Tj1) )
	{
		m_jerk = 0.0;
	} else if( t >= (Ta - Tj1) && t < Ta )
	{
		m_jerk = minj;
	} else if( t >= Ta && t < (Ta + Tv))
	{
		m_jerk = 0.0;
	} else if(t >= (T - Td) && t < (T - Td + Tj2))
	{
		m_jerk = minj;
	} else if(t >= (T - Td + Tj2 ) && t < (T - Tj2))
	{
		m_jerk = 0.0;
	} else
	{
		m_jerk = maxj;
	}

	m_jerk = m_jerk * s;
	printf("%f, \n", m_jerk);
	return m_jerk;	
}


double VelocityProfile_S::Duration() const {
	return T;
}

VelocityProfile_S::~VelocityProfile_S() {}


void VelocityProfile_S::Write(std::ostream& os) const {
	os << "TRAPEZOIDAL[" << maxvel << "," << maxacc <<"]";
}


