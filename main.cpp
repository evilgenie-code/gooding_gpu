#include <cmath>
#include <cstdio>
#include <algorithm>
#include <iostream>

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))


void getUnit(double* r, double R, double* Unit)
{
    for (int j = 0; j < 3; j++)
    {
        Unit[j] = r[j] / R;
    }
}

void getUnitDouble(double* r, double R, double* Unit)
{
    for (int j = 0; j < 3; j++)
    {
        Unit[j] = r[j] / R;
    }
}

double getVetVal(double* vet1)
{
    double temp = 0.0;

    for (int i = 0; i < 3; i++)
    {
        temp += vet1[i] * vet1[i];
    }

    return sqrt(temp);
}

void cross(double* vet1, double* vet2, double* prod)
{
    prod[0] = (vet1[1] * vet2[2] - vet1[2] * vet2[1]);
    prod[1] = (vet1[2] * vet2[0] - vet1[0] * vet2[2]);
    prod[2] = (vet1[0] * vet2[1] - vet1[1] * vet2[0]);
}

// series approximation to T(x) and its derivatives
// (used for near - parabolic cases)
void sigmax(double y, double sig = 0.0, double dsigdx = 0.0, double d2sigdx2 = 0.0, double d3sigdx3 = 0.0)
{
    double powers[25];
    int i;
    // preload the factors[an]
    // (25 factors is more than enough for 16 - digit accuracy)
    double an[25] = { 4.000000000000000e-001,     2.142857142857143e-001,     4.629629629629630e-002,
                      6.628787878787879e-003,     7.211538461538461e-004,     6.365740740740740e-005,
                      4.741479925303455e-006,     3.059406328320802e-007,     1.742836409255060e-008,
                      8.892477331109578e-010,     4.110111531986532e-011,     1.736709384841458e-012,
                      6.759767240041426e-014,     2.439123386614026e-015,     8.203411614538007e-017,
                      2.583771576869575e-018,     7.652331327976716e-020,     2.138860629743989e-021,
                      5.659959451165552e-023,     1.422104833817366e-024,     3.401398483272306e-026,
                      7.762544304774155e-028,     1.693916882090479e-029,     3.541295006766860e-031,
                      7.105336187804402e-033 };

    for (i = 0; i < 25; i++)
    {
        //powers of y
        powers[i] = pow(y, i + 1);
        sig += 4 / 3 + powers[i] * an[i];
    }

    for (i = 0; i < 24; i++)
    {
        // dsigma / dx
        dsigdx += ((i + 2) * powers[i + 1]) * an[i + 1];
    }

    for (i = 0; i < 23; i++)
    {
        // d2sigma / dx2
        d2sigdx2 = ((i + 3) * (i + 2) * powers[i + 2]) * an[i + 2];
    }

    for (i = 0; i < 23; i++)
    {
        // d3sigma / dx3
        d3sigdx3 = ((i + 4) * (i + 3) * (i + 2) * powers[i + 3]) * an[i + 3];
    }

    dsigdx = dsigdx + an[0];
    d2sigdx2 = d2sigdx2 + 0 + 2 * an[1];
    d3sigdx3 = d3sigdx3 + 0 + 0 + 3 * 2 * an[2];

}

void LancasterBlanchard(double x, double q, double m, double* Tprod)
{
    double E, sig1 = 0.0, sig2 = 0.0, dsigdx1 = 0.0,
            dsigdx2 = 0.0, d2sigdx21 = 0.0, d2sigdx22 = 0.0, d3sigdx31 = 0.0, d3sigdx32 = 0.0,
            y, z, f, g, d;
    // protection against idiotic input
    if (x < -1) //impossible; negative eccentricity
    {
        x = abs(x) - 2;
    }
    if (x == -1) // impossible; offset x slightly
    {
        x = x + 0.00001;
    }

    // compute parameter E
    E = x * x - 1;

    // T(x), T'(x), T''(x)
    if (x == 1) // exactly parabolic; solutions known exactly
    {
        Tprod[0] = 4 / 3 * (1 - pow(q, 3));
        Tprod[1] = 4 / 5 * (pow(q, 5) - 1);
        Tprod[2] = Tprod[1] + 120 / 70 * (1 - pow(q, 7));
        Tprod[3] = 3 * (Tprod[2] - Tprod[1]) + 2400 / 1080 * (pow(q, 9) - 1);
    }
    else if (fabs(x - 1) < 1e-2) // near - parabolic; compute with series
        // evaluate sigma
    {
        sigmax(-E, sig1, dsigdx1, d2sigdx21, d3sigdx31);
        sigmax(-E * q * q, sig2, dsigdx2, d2sigdx22, d3sigdx32);
        Tprod[0] = sig1 - pow(q, 3) * sig2;
        Tprod[1] = 2 * x * (pow(q, 5) * dsigdx2 - dsigdx1);
        Tprod[2] = Tprod[1] / x + 4 * pow(x, 2) * (d2sigdx21 - pow(q, 7) * d2sigdx22);
        Tprod[3] = 3 * (Tprod[2] - Tprod[1] / x) / x + 8 * x * x * (pow(q, 9) * d3sigdx32 - d3sigdx31);
    }
    else // all other cases
    {
        // compute all substitution functions
        y = sqrt(fabs(E));
        z = sqrt(1 + pow(q, 2) * E);
        f = y * (z - q * x);
        g = x * z - q * E;
        d = (E < 0) * (atan2(f, g) + M_PI * m) + (E > 0) * logf(MAX(0, f + g));
        Tprod[0] = 2 * (x - q * z - d / y) / E;
        Tprod[1] = (4 - 4 * pow(q, 3) * x / z - 3 * x * Tprod[0]) / E;
        Tprod[2] = (-4 * pow(q, 3) / z * (1 - pow(q, 2) * pow(x, 2) / pow(z, 2)) - 3 * Tprod[0] - 3 * x * Tprod[1]) / E;
        Tprod[3] = (4 * pow(q, 3) / pow(z, 2) * ((1 - pow(q, 2) * pow(x, 2) / pow(z, 2)) + 2 * pow(q, 2) * x / pow(z, 2) * (z - x)) - 8 * Tprod[1] - 7 * x * Tprod[2]) / E;

    }
}

double getXM(double phr, int revs)
{
    double xMpi = 4 / (3 * M_PI *(2 * revs + 1));

    if (phr < M_PI)
    {
        return xMpi * pow(phr / M_PI, 0.125);
    }
    else if (phr > M_PI)
    {
        return xMpi * (2.0 - pow(2 - (phr / M_PI), 0.125));
    }

    return 0;

}

double halleysMethod(double xM, double q, int revs)
{
    int iterations;
    const double tol = 1e-11;
    double tProd[4] = { 0.0, 0.0, 0.0, 0.0 }, Tp = 0, Tpp, Tppp, xMp;

    while (fabs(Tp) > tol)
    {
        iterations = iterations + 1;

        LancasterBlanchard(xM, q, revs, tProd);

        Tp = tProd[1];
        Tpp = tProd[2];
        Tppp = tProd[3];

        xMp = xM;
        xM = xM - 2 * Tp * Tpp / (2 * pow(Tpp, 2) - Tp * Tppp);

        if (iterations % 7 == 0)
        {
            xM = (xMp + xM) / 2;
        }
//
//        if (iterations > 25)
//        {
//            return ;
//        }
    }

    return xM;
}

void lamdert(const float* r1, const float* r2, float t, int revs, int lw, double mu)
{
    int j, leftbranch = (mu > 0) - (mu < 0);
    double R1 = 0.0, R2 = 0.0, dotProd = 0.0, mcrsProd, theta, c, s, T, q,
            T0 = 0.0, dummy = 0.0;
    double  r1Unit[3], r2Unit[3], crsProd[3], ucrsProd[3], th1Unit[3], th2Unit[3],
            r1Double[3], r2Double[3], Td, phr, v1[3], v2[3], Tprod[4];

    const double tol = 1e-11;

    if (t <= 0)
    {
        printf("ERROR in Lambert Solver: Negative Time in input.\n");
        return;
    }

    for (j = 0; j < 3; j++)
    {
        r1Double[j] = r1[j];
        r2Double[j] = r2[j];
    }

    for (j = 0; j < 3; j++)
    {
        R1 += pow(r1Double[j], 2);
        R2 += pow(r2Double[j], 2);
    }

    R1 = sqrt(R1);
    R2 = sqrt(R2);

    for (j = 0; j < 3; j++)
    {
        r1Unit[j] = r1Double[j] / R1;
        r2Unit[j] = r1Double[j] / R2;
    }

    cross(r1Double, r2Double, crsProd); // перекрестное произведение

    mcrsProd = getVetVal(crsProd); //Величина этого перекрестного произведений

    getUnitDouble(crsProd, mcrsProd, ucrsProd);

    cross(ucrsProd, r1Unit, th1Unit);
    cross(ucrsProd, r2Unit, th2Unit);

    for (j = 0; j < 3; j++)
    {
        dotProd += (r1Double[j] * r2Double[j]);
    }

    theta = acos(MAX(-1, MIN(1, dotProd / (R1 * R2))));

    if (lw == 1)
    {
        theta = theta - 2 * M_PI;
    }

    c = sqrt(pow(R1, 2) + pow(R2, 2) - 2.0 * R1 * R2 * cos(theta));
    s = (R1 + R2 + c) / 2;
    T = sqrt(8.0 * mu / pow(s, 3)) * t;
    q = (sqrt(R1 * R2) / s) * cos(theta / 2);

    double Tprod0[4] = {T0, 0.0, 0.0, 0.0 };

    LancasterBlanchard(0, q, revs, Tprod0);

    T0 = Tprod0[0];

    phr = fmod(2.0 * atan2(1.0 - pow(q, 2), 2.0 * q), 2.0 * M_PI);

    Td = T0 - T;

    // initial output is pessimistic
//    for (int i = 0;i < 3;i++)
//    {
//        v1[i] = 0.0;
//        v2[i] = 0.0;
//    }

    double x01;
    double x0, x02, x03, W, lambda, xMpi, xM, Tp;

    if (revs == 0)
    {
        x01 = T0 * Td / 4 / T;

        if (Td > 0)
        {
            x0 = x01;
        }
        else
        {
            x01 = Td / (4 - Td);
            x02 = -sqrt(-Td / (T + T0 / 2));
            W = x01 + 1.7 * sqrt(2 - phr / M_PI);

            if (W >= 0)
            {
                x03 = x01;
            }
            else
            {
                x03 = x01 + pow(-W, (1.0 / 16)) * (x02 - x01);
            }

            lambda = 1 + x03 * (1 + x01) / 2 - 0.03 * pow(x03, 2) * sqrt(1 + x01);
            x0 = lambda * x03;
        }

        if (x0 < -1)
        {
            return;
        }
    }
    else
    {
        //multi-rev cases
        // determine minimum Tp(x)
        xM = getXM(phr, revs);
        Tp = 1.0;

        // use Halley's method
        xM = halleysMethod(xM, q, revs);

        // xM should be elliptic(-1 < x < 1)
        // (this should be impossible to go wrong)
        if ((xM < -1) || (xM > 1))
        {
            return;
        }

        //corresponding time
        LancasterBlanchard(xM, q, revs, Tprod);

        double TM = Tprod[0];
        double Tpp;

        // T should lie above the minimum T
        if (TM > t)
        {
            return;
        }

        //find two initial values for second solution(again with lambda - type patch)

        // some initial values


        double TmTM = T - TM;
        double  T0mTM = T0 - TM;

        LancasterBlanchard(xM, q, revs, Tprod);
        Tpp = Tprod[2];

        if (leftbranch > 0)
        {
            double x = sqrt(TmTM / (Tpp / 2 + TmTM / pow((1 - xM), 2)));
            W = xM + x;
            W = 4.0 * W / (4.0 + TmTM) + pow((1 - W), 2);
            x0 = x * (1 - (1 + revs + (theta - 0.5)) / (1 + 0.15 * revs) * x * (W*0.5 + 0.03 * x * sqrt(W))) + xM;

            if (x0 > 1)
            {
                return;
            }
        }
        else
        {
            if (Td > 0)
            {
                x0 = xM - sqrt(TM / (Tpp / 2 - TmTM * (Tpp / 2 / T0mTM - 1 / pow(xM, 2))));
            }
            else
            {
                double x00 = Td / (4.0 - Td);
                W = x00 + 1.7 * sqrt(2 * (1 - phr));
                if (W >= 0)
                {
                    x03 = x00;
                } else
                {
                    x03 = x00 - sqrt(pow((-W), 0.125) * (x00 + sqrt(-Td / (1.5 * T0 - Td))));
                }

                W = 4 / (4 - Td);
                lambda = (1 + (1 + revs + 0.24 * (theta - 0.5)) /
                              ((1 + 0.15 * revs) * x03 * (0.5 * W - 0.03 * x03 * sqrt(W))));
                x0 = x03 * lambda;
            }
        }

        if (x0 < -1)
        {
            return;
        }
    }

    //find root of Lancaster & Blancard's function
    //(Halley's method)
    printf("x = %f \n", x0);
    double x = x0;
    double Tx = 99999;
    int iterations = 0;

    while (std::abs(Tx) > tol)
    {
        // compute function value, and first two derivatives

        LancasterBlanchard(x, q, revs, Tprod);
        Tx = Tprod[0];
        Tp = Tprod[1];
        double Tpp = Tprod[2];

        // find the root of the *difference* between the
        // function value[T_x] and the required time[T]
        Tx = Tx - T;

        // new value of x
        double xp = x;
        x = x - (2 * Tx * Tp) / (2 * pow(Tp, 2) - Tx * Tpp);
        if (iterations % 7 == 0)
        {
            x = (xp + x) / 2;
        }

        if (iterations > 25)
        {
            return;
        }

        iterations++;
//        printf("x = %f \n ", x);
    }


    double gamma = sqrt(mu * s / 2);

    double sigma, rho, z;

    if (c == 0)
    {
        sigma = 1;
        rho = 0;
        z = std::abs(x);
    }
    else
    {
        sigma = 2 * sqrt(R1 * R2 / pow(c, 2)) * sin(theta / 2);
        rho = (R1 - R2) / c;
        z = sqrt(1.0 + (pow(q, 2) * (pow(x, 2) - 1)));
    }


    // radial component

    double Vr1 = +gamma * ((q*z - x) - rho * (q*z + x)) / R1;
    double Vr2 = -gamma * ((q*z - x) + rho * (q*z + x)) / R2;

    //tangential component

    double Vtan1 = sigma * gamma * (z + q * x) / R1;
    double Vtan2 = sigma * gamma * (z + q * x) / R2;

    double Vr1vec[3], Vr2vec[3], Vtan1vec[3], Vtan2vec[3];

    for (j = 0; j < 3; j++)
    {
        Vr1vec[j] = Vr1 * r1Unit[j];
        Vr2vec[j] = Vr2 * r2Unit[j];
        Vtan1vec[j] = Vtan1 * th1Unit[j];
        Vtan2vec[j] = Vtan2 * th2Unit[j];

        //Cartesian velocity
        v1[j] = Vtan1vec[j] + Vr1vec[j];
        v2[j] = Vtan2vec[j] + Vr2vec[j];

        printf("v1[%i] = %f, v2[%i] = %f \n", j, v1[j], j, v2[j]);
    }
}

int main()
{
    double AU = 1.49597870691e8;
    double fMSun = 1.32712440018e11;             // km^3/sec^2

    double UnitR = AU;
    double UnitV = sqrt(fMSun / UnitR);          // km/sec
    double UnitT = (UnitR / UnitV) / 86400;         // day

    float mu = 1.0;
    float t = 100.0 / UnitT;
    int revs = 0;
    int lw = 1;
    float r1[3] = { -7.8941608095246896e-01, -6.2501194900473045e-01, 3.5441335698377735e-05 };
    float r2[3] = { 1.3897892184188783e+00, 1.3377137029002054e-01, -3.1287386211010106e-02 };

    lamdert(r1, r2, t, revs, lw, mu);
}
