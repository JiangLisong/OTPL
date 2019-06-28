#ifndef _VELOCITYPROFILE_S_H
#define _VELOCITYPROFILE_S_H

// #ifdef  __cplusplus
// extern "C"
// {
// #endif//__cplusplus

#include <iostream>
	/*
	 * A double s VelocityProfile implementation.
	 */
class VelocityProfile_S
	{
		// For "running" a motion profile :
		double Ta, Tv, Td, Tj1, Tj2, Tj, T; // duration t of each part

		// specification of the motion profile :
		double maxvel, minvel;
		double maxacc, minacc;
		double maxj, minj;
		double startpos, endpos;
		double startvel, endvel;
		double v_lim, a_lima, a_limd;
		int s;
	public:

		VelocityProfile_S();
		// constructs motion profile class with <maxvel> and <maxacc> as parameters of the
		// trajectory.

		virtual void SetProfile(double pos1,double pos2, double vel1, double vel2);

		virtual void ExtendDuration(
			double pos1, double pos2, double vel1, double vel2, double newduration
		);

        virtual void SetMax(double _maxvel, double _maxacc, double _maxjerk);
		virtual double Duration() const;
		virtual double Pos(double t) const;
		virtual double Vel(double t) const;
		virtual double Acc(double t) const;
		virtual double Jerk(double t) const;
		virtual void Write(std::ostream& os) const;

		virtual ~VelocityProfile_S();
	};

// #ifdef  __cplusplus
// }
// #endif//__cplusplus

#endif
