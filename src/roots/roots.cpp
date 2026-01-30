#include "roots.hpp"
#include <cmath>
#include <limits>

constexpr double min = 1e-6;
constexpr int max = 1e6;

bool bisection(std::function<double(double)> f, double a, double b, double *root)
{
    double fa = f(a);
    double fb = f(b);
    if (fa * fb > 0) return false;

    for (int i = 0; i < max; ++i)
    {
        double mid = (a + b) / 2;
        double fm = f(mid);

        if (std::abs(fm) < min)
        {
            *root = mid;
            return true;
        }

        if ((b - a) / 2 < min)
        {
            *root = mid;
            return true;
        }

        if (fa * fm < 0)
        {
            b = mid;
            fb = fm;
        }
        else
        {
            a = mid;
            fa = fm;
        }
    }
    *root = (a + b) / 2;
    return true;
}

bool regula_falsi(std::function<double(double)> f, double a, double b, double *root)
{
    double fa = f(a);
    double fb = f(b);
    if (fa * fb > 0) return false;
                           
    double c = a;
    for (int i = 0; i < max; ++i)
    {
        c = (a * fb - b * fa) / (fb - fa);
        double fc = f(c);

        if (std::abs(fc) < min)
        {
            *root = c;
            return true;
        }

        if (fa * fc < 0)
        {
            b = c;
            fb = fc;
        }
        else
        {
            a = c;
            fa = fc;
        }
    }
    *root = c;
    return true;
}

bool newton_raphson(std::function<double(double)> f, std::function<double(double)> g, double a, double b, double c, double *root)
{
    double x = c;
    for (int i = 0; i < max; ++i)
    {
        double fx = f(x);
        double gx = g(x);
  
        double x_next = x - fx / gx;

        if (x_next < a)
            return false;
           
        if (x_next > b)
            return false;    

        if (std::abs(x_next - x) < min)
        {
            *root = x_next;
            return true;
        }
        x = x_next;
    }
    return true;
}


bool secant(std::function<double(double)> f, double a, double b, double c, double *root)
{
    double x0 = c;
    double x1 = b;
    double f0 = f(x0);
    double f1 = f(x1);

    for (int i = 0; i < max; ++i)
    {

        double x2 = x1 - f1 * (x1 - x0) / (f1 - f0);

        if (x2 < a)
            return false;

        if (x2 > b)
            return false;    

        if (std::abs(x2 - x1) < min)
        {
            *root = x2;
            return true;
        }

        x0 = x1;
        f0 = f1;
        x1 = x2;
        f1 = f(x1);
    }

    return true;
}