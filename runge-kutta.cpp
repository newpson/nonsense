#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

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

/// Evaluates expression at x+h
double estimate(const std::string &expression,
                const Evaluator::library_t &library,
                const double x, const double y, const double h)
{
    const double k1 = eval(expression, library, x);
    const double k2 = eval(expression, library, x + h/2.0);
    const double k3 = k2;
    const double k4 = eval(expression, library, x + h);
    return y + h/6.0 * (k1 + 2.0*k2 + 2.0*k3 + k4);
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


    const double dir = std::copysign(1.0, b - a);

    double h = 0.0;
    const double signed_length = b - a;
    const double length = std::abs(signed_length);
    if (length < 10.0 * hmin) {
        if (length < hmin) {
            std::cerr << "error: range is less than hmin" << std::endl;
            return points;
        }
        h = dir * hmin;
    } else {
        h = signed_length/10.0;
    }

    double x = a;
    double y = ya;
    // points.push_back(Vector2D(x, y));
    std::cout << Vector2D(x, y) << std::endl;

    double signed_remain = b - x;
    double remain = std::abs(signed_remain);
    while (dir * signed_remain > 0) {
        if (remain < 2.0*std::abs(h)) {
            h = signed_remain;
            if (remain >= 2.0*hmin)
                h /= 2.0;
        }

        double y_next = estimate(expression, library, x, y, h);
        double y_next_half = estimate(expression, library, x, y, h/2.0);
        while (error(y_next, y_next_half) > eps && std::abs(h) >= 2.0*hmin) {
            h /= 2.0;
            y_next = y_next_half;
            y_next_half = estimate(expression, library, x, y, h/2.0);
        }

        std::cerr << "info: h = " << h << ", eps = " << error(y_next, y_next_half) << std::endl;
        y = y_next;
        x += h;
        std::cout << Vector2D(x, y) << std::endl;
        // points.push_back(Vector2D(x, y));

        signed_remain = b - x;
        remain = std::abs(signed_remain);
    }

    return points;
}

double antiderivative(double x)
{
    return -std::sin(2.0*x)/4.0 + x*std::log(x) - x/2.0;
}

int main(int argc, char **argv)
{
    Evaluator::library_t library = {{"x", std::shared_ptr<Evaluable>(new Variable())}};
    library.insert(Evaluator::default_library.begin(),
                   Evaluator::default_library.end());

    const double A = 0.101;
    const double B = 0.100;
    auto points = integrate("sin(x)^2 + log(x)", library, A, B, antiderivative(A), 1e-03, 1e-03);
    // for (const auto &point: points)
    //     std::cout << point << std::endl;
}
