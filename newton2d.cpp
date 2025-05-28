#include <cmath>
#include <iostream>
#include <string>
#include <stdexcept>

#include "evaluator.h"
#include "gnuplot-helpers.h"
#include "variable.h"

#include "Eigen/Dense"

double eval(
    const std::string &expr,
    const Evaluator::library_t &lib,
    const Eigen::Vector2d &x)
{
    dynamic_cast<Variable &>(*lib.at("x")).value = x.x();
    dynamic_cast<Variable &>(*lib.at("y")).value = x.y();
    return Evaluator::eval(expr, lib);
}

double dx(const Eigen::Matrix3d &z_grid, const double h)
{
    const double right = z_grid.coeff(1, 2);
    const double left = z_grid.coeff(1, 0);
    return (right - left) / 2.0 / h;
}

double dy(const Eigen::Matrix3d &z_grid, const double h)
{
    const double top = z_grid.coeff(0, 1);
    const double bottom = z_grid.coeff(2, 1);
    return (top - bottom) / 2.0 / h;
}

double ddy(const Eigen::Matrix3d &z_grid, const double h)
{
    const double bottom = z_grid.coeff(2, 1);
    const double mid = z_grid.coeff(1, 1);
    const double top = z_grid.coeff(0, 1);
    return (bottom - 2.0*mid + top) / h / h;
}

double ddx(const Eigen::Matrix3d &z_grid, const double h)
{
    const double left = z_grid.coeff(1, 0);
    const double mid = z_grid.coeff(1, 1);
    const double right = z_grid.coeff(1, 2);
    return (left - 2.0*mid + right) / h / h;
}

double dxdy(const Eigen::Matrix3d &z_grid, const double h)
{
    const double top_left = z_grid.coeff(0, 0);
    const double top_right = z_grid.coeff(0, 2);
    const double bottom_left = z_grid.coeff(2, 0);
    const double bottom_right = z_grid.coeff(2, 2);
    return (top_right - bottom_right + bottom_left - top_left) / 4.0 / h / h;
}

inline
double dydx(const Eigen::Matrix3d &z_grid, const double h)
{
    return dxdy(z_grid, h);
}

Eigen::Matrix2d inverse(const Eigen::Matrix2d &m)
{
    bool inversible = false;
    Eigen::Matrix2d inverted;
    m.computeInverseWithCheck(inverted, inversible);
    if (!inversible)
        throw std::invalid_argument("matrix is not inversible");
    return inverted;
}

void find_min(
    const std::string &expr,
    const Evaluator::library_t &lib,
    const Eigen::Vector2d x0, const double eps = 1e-06)
{
    const double h = 0.001; // TODO dynamic step depending on epsilon (?)

    Eigen::Vector2d x = x0;
    double err = 0.0;
    int i = 0;
    double rate = 1.0;
    do {
begin:
        Eigen::Matrix3d z_grid;
        /* -h,+h  +0,+h  +h,+h
         * -h,+0  +0,+0  +h,+0
         * -h,-h  +0,-h  +h,-h
         */
        for (int dy = 1; dy >= -1; --dy)
            for (int dx = -1; dx <= 1; ++dx)
                z_grid.coeffRef(1 - dy, dx + 1) = eval(expr, lib, {x.x() + dx*h, x.y() + dy*h});

        Eigen::Matrix2d hessian;
        hessian << ddx(z_grid, h), dxdy(z_grid, h),
                   dydx(z_grid, h), ddy(z_grid, h);
        // std::cerr << "info: hessian:" << std::endl;
        // std::cerr << hessian << std::endl;
        const Eigen::Matrix2d inv_hessian = inverse(hessian);
        const Eigen::Vector2d grad(dx(z_grid, h), dy(z_grid, h));
        const Eigen::Vector2d dir = -inv_hessian * grad;

        // // Armijo–Goldstein condition to make the algorithm more stable
        // const double gap = 0.5;
        // const double z = z_grid.coeff(1, 1);
        // double rate = 2.0;
        // double z_shifted = 0.0;
        // int j = 0;
        // do {
        //     rate /= 2.0;
        //     const Eigen::Vector2d x_shifted = x + rate*dir;
        //     z_shifted = eval(expr, lib, x_shifted);
        //     ++j;
        //     // std::cerr << "shif: " << X_shifted.x() << " " << X_shifted.y() << " " << z_shifted << std::endl;
        //     std::cerr << "info: rate: " << rate << std::endl;
        // } while ((z - z_shifted < rate * (-gap * grad.dot(dir))) && j < 60);

        const Eigen::Vector2d x_next = x + rate*dir;
        if (eval(expr, lib, x_next) > eval(expr, lib, x)) {
            rate /= 2.0;
            goto begin;
        }

        err = (x - x_next).norm();
        x = x_next;
        std::cout << x.x() << " " << x.y() << " " << eval(expr, lib, x) << std::endl;
        ++i;
    } while (err > eps && i < 100);
    std::cerr << "info: i = " << i << std::endl;

    // std::cerr << x << " " << y << " " << eval(expr, lib, x) << std::endl;
}

int main(const int argc, const char **argv)
{
    if (argc < 1 + 3) {
        std::cerr << "<expression> <x0> <y0> [eps]" << std::endl;
        return 0;
    }

    std::shared_ptr<Variable> x(new Variable());
    std::shared_ptr<Variable> y(new Variable());
    Evaluator::library_t library = {{"x", x}, {"y", y}};
    library.insert(Evaluator::default_library.begin(),
                   Evaluator::default_library.end());

    const std::string expression(argv[1]);
    const double x0 = std::strtod(argv[2], nullptr);
    const double y0 = std::strtod(argv[3], nullptr);
    find_min(expression, library, {x0, y0});

    // Vector2D minimum;
    // if (argc > 1 + 2) {
    //     const double eps = std::strtod(argv[3], nullptr);
    //     minimum = find_min(expression, library, x0, eps);
    // } else {
    //     minimum = find_min(expression, library, x0);
    // }

    // std::cerr << minimum << std::endl;

    return 0;
}
