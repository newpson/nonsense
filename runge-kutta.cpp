#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <stack>

#include "evaluator.h"
#include "gnuplot-helpers.h"
#include "math.h"
#include "variable.h"

double eval(const std::string &expression,
            const Evaluator::library_t &library,
            const double x)
{
    dynamic_cast<Variable &>(*library.at("x")).value = x;
    return Evaluator::eval(expression, library);
}

double estimate(const std::string &expression,
                const Evaluator::library_t &library,
                const Vector2D &v, const double h)
{
    const double k1 = eval(expression, library, v.x);
    const double k2 = eval(expression, library, v.x + h/2.0);
    const double k3 = k2;
    const double k4 = eval(expression, library, v.x + h);
    return v.y + h/6.0 * (k1 + 2.0*k2 + 2.0*k3 + k4);
}

double error(const double y, const double y_half)
{
    return std::abs(y - y_half)/15.0; // (2^p - 1); p = 4
}

std::vector<Vector2D> integrate(const std::string &expression,
                                const Evaluator::library_t &library,
                                const double a, const double b,
                                const double ya, double eps, double hmin = 1e-06)
{
    std::vector<Vector2D> points;
    std::stack<double> prohibitors;

    const double hcrit = 2.0*6.0*hmin;
    const double dir = std::copysign(1.0, b - a);
    std::cerr << "info: hcrit = " << hcrit << ", dir = " << dir << std::endl;

    const double length_s = b - a;
    const double length_u = std::abs(length_s);
    double h = 0.0;
    if (length_u < 10.0 * hcrit) {
        if (length_u < hcrit) {
            std::cerr << "error: range is too small OR hmin is too big" << std::endl;
            return points;
        }
        h = dir * hcrit;
    } else {
        h = length_s / 10.0;
    }
    std::cerr << "info: init h = " << h << std::endl;

    Vector2D cur(a, ya);
    // points.push_back(Vector2D(x, y));
    std::cout << cur << std::endl;

    double remain_s = b - cur.x;
    double remain_u = std::abs(remain_s);
    while (dir * remain_s > 0) {
        if (remain_u < 2.0*std::abs(h)) {
            std::cerr << "warning: remain is low" << std::endl;
            if (remain_u < hcrit) { // can't happen (theoretically)
                std::cerr << "error: failed to reach the end of the range" << std::endl;
                return points;
            }
            h = remain_s;
            if (remain_u >= 2.0*hcrit)
                h /= 2.0;
        }

        // TODO approximation stack
        double y = estimate(expression, library, cur, h);
        double y_half = estimate(expression, library, cur, h/2.0);
        while (error(y, y_half) < eps/16.0) {
            h *= 2.0;
            y_half = y;
            y = estimate(expression, library, cur, h);
            std::cerr << "info: *2, h = " << h << ", eps = " << error(y, y_half) << std::endl;
        }
        while (error(y, y_half) > eps && std::abs(h) >= 2.0*hcrit) {
            h /= 2.0;
            y = y_half;
            y_half = estimate(expression, library, cur, h/2.0);
            std::cerr << "info: /2, h = " << h << ", eps = " << error(y, y_half) << std::endl;
        }

        cur.y = y;
        cur.x += h;
        std::cout << cur << std::endl;
        std::cerr << "info: fin, h = " << h << ", eps = " << error(y, y_half) << std::endl;
        // // points.push_back(Vector2D(x, y));

        remain_s = b - cur.x;
        remain_u = std::abs(remain_s);

        std::cerr << "info: ---" << std::endl;
    }

    return points;
}

double antiderivative(const double x)
{
    return -std::sin(2.0*x)/4.0 + x*std::log(x) - x/2.0;
}

int main(int argc, char **argv)
{
    Evaluator::library_t library = {{"x", std::shared_ptr<Evaluable>(new Variable())}};
    library.insert(Evaluator::default_library.begin(),
                   Evaluator::default_library.end());

    const double A = 0.1;
    const double B = 0.9;
    auto points = integrate("sin(x)^2 + log(x)", library, A, B, antiderivative(A), 1e-03);
    // for (const auto &point: points)
    //     std::cout << point << std::endl;
}
