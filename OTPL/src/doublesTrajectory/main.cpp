#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdint.h>

#include "velocityprofile_s.h"

int main(int argc, const char *argv[])
{
    double q0 = 10.0;
    double q1 = 0.0;
    double v0 = 1.0;
    double v1 = 0.0;
    double v_max = 5.0;
    double a_max = 10.0;
    double j_max = 30.0;

    VelocityProfile_S *pVelocityProfile = new VelocityProfile_S();
    pVelocityProfile->SetMax(v_max, a_max, j_max);
    pVelocityProfile->SetProfile(q0, q1, v0, v1);
    double deltaT = 0.005;
    int32_t i = 0;
    double t = 0.0;
    while(t <= pVelocityProfile->Duration())
    {
        //printf("%f, %f, ", t, pVelocityProfile->Duration());
        pVelocityProfile->Pos(t);
        pVelocityProfile->Vel(t);
        pVelocityProfile->Acc(t);
        pVelocityProfile->Jerk(t);
        i++;
        t = i * deltaT;
    }
    delete pVelocityProfile;
    return 0;
}
