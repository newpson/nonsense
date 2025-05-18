#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <cmath>

#include "evaluator.h"
#include "gnuplot-helpers.h"
#include "math.h"
#include "variable.h"

double eval(const std::string &expression,
            const Evaluator::library_t &library,
            const double x, const double y)
{
    dynamic_cast<Variable &>(*library.at("x")).value = x;
    dynamic_cast<Variable &>(*library.at("y")).value = y;
    return Evaluator::eval(expression, library);
}

double estimate(const std::string &expression,
                const Evaluator::library_t &library,
                const Vector2D &v, const double h)
{
    const double k1 = eval(expression, library, v.x, v.y);
    const double k2 = eval(expression, library, v.x + h/2.0, v.y + h*k1/2.0);
    const double k3 = eval(expression, library, v.x + h/2.0, v.y + h*k2/2.0);
    const double k4 = eval(expression, library, v.x + h, v.y + h*k3);
    return v.y + h/6.0 * (k1 + 2.0*k2 + 2.0*k3 + k4);
}

double inline error(const double y, const double y_half)
{
    return std::abs(y - y_half)/15.0; // (2^s - 1); s = 4
}

void integrate(const std::string &expression,
                                const Evaluator::library_t &library,
                                const double a, const double b,
                                const double ya, double eps = 1e-03, double hmin = 1e-15)
{
    const double hcrit = 2.0*6.0*hmin;
    const double dir = std::copysign(1.0, b - a);
    std::cerr << "info: hcrit = " << hcrit << ", dir = " << dir << std::endl;

    const double length_s = b - a;
    const double length_u = std::abs(length_s);
    double h = 0.0;
    if (length_u < 10.0 * hcrit) {
        if (length_u < hcrit) {
            std::cerr << "error: range is too small OR hmin is too big" << std::endl;
            return;
        }
        h = dir * hcrit;
    } else {
        h = length_s / 10.0;
    }
    std::cerr << "info: init h = " << h << std::endl;

    Vector2D cur(a, ya);
    std::cout << cur << std::endl;

    int i = 0;
    int num_hcrit = 0;
    int num_eps = 0;

    double remain_s = b - cur.x;
    double remain_u = std::abs(remain_s);
    bool big_remain = true;
    double err = 0.0;
    while (dir * remain_s > 0) {
        if (remain_u < 2.0*std::abs(h) && big_remain) {
            std::cerr << "warning: remain is low: " << remain_u << std::endl;
            if (remain_u < hcrit) {
                std::cerr << "error: failed to reach the end of the range: " << remain_u << std::endl;
                break;
            }
            h = remain_s;
            if (remain_u >= 2.0*hcrit && err > eps)
                h /= 2.0;
            big_remain = false;
        }

        double y = estimate(expression, library, cur, h);
        double y_half = estimate(expression, library, cur, h/2.0);
        double y_dhalf = estimate(expression, library, {cur.x + h/2.0, y_half}, h/2.0);
        err = error(y, y_dhalf);

        while (err < eps/16.0 && 2.0*std::abs(h) < remain_u && big_remain) {
            h *= 2.0;
            y_half = y;
            y = estimate(expression, library, cur, h);
            y_dhalf = estimate(expression, library, {cur.x + h/2.0, y_half}, h/2.0);
            err = error(y, y_dhalf);
            std::cerr << "info: *2: h = " << h << ", err = " << err << std::endl;
        }
        while (err > eps && std::abs(h) >= 2.0*hcrit && big_remain) {
            h /= 2.0;
            y = y_half;
            y_half = estimate(expression, library, cur, h/2.0);
            y_dhalf = estimate(expression, library, {cur.x + h/2.0, y_half}, h/2.0);
            err = error(y, y_dhalf);
            std::cerr << "info: /2: h = " << h << ", err = " << err << std::endl;
        }

        cur = {cur.x + h, y};
        if (h == hcrit)
            ++num_hcrit;
        if (err > eps)
            ++num_eps;

        std::cout << cur << std::endl;
        std::cerr << "info: fin: h = " << h << ", err = " << err << std::endl;

        ++i;
        remain_s = b - cur.x;
        remain_u = std::abs(remain_s);
        std::cerr << "info:" << std::endl;
    }

    std::cerr << "info: stats: i = " << i
              << ", num_eps = " << num_eps
              << ", num_hcrit = " << num_hcrit
              << std::endl;
}

int main(int argc, char **argv)
{
    /* Примеры
     * y' = 1 или 2*x или 3*x^2 или 4*x^3, a = 0, b =1, ya = 0; погрешность нулевая
     * y' = 2*y - sin(x), a = 0, b = 1, ya = 0; y = (2*sin(x) + cos(x) - e^(2*x))/5; y^(V) = -(2*sin(x))/5-(cos(x))/5-(64*e^(2*x))/5
     */
    if (argc < 1 + 4) {
        std::cerr << "<expression> <A> <B> <C> [eps [hmin]]" << std::endl;
        return 0;
    }

    Evaluator::library_t library = {
        {"x", std::shared_ptr<Evaluable>(new Variable())},
        {"y", std::shared_ptr<Evaluable>(new Variable())}
    };
    library.insert(Evaluator::default_library.begin(),
                   Evaluator::default_library.end());

    const std::string expression(argv[1]);
    const double a = std::strtod(argv[2], nullptr);
    const double b = std::strtod(argv[3], nullptr);
    const double c = std::strtod(argv[4], nullptr);
    if (argc > 1 + 4) {
        const double eps = std::strtod(argv[5], nullptr);
        if (argc > 1 + 5) {
            const double hmin = std::strtod(argv[6], nullptr);
            integrate(expression, library, a, b, c, eps, hmin);
        } else {
            integrate(expression, library, a, b, c, eps);
        }
    } else {
        integrate(expression, library, a, b, c);
    }
}
