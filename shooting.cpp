#include <iostream>
#include <memory>
#include <string>

#include "evaluator.h"
#include "gnuplot-helpers.h"
#include "math.h"
#include "variable.h"

#include "Eigen/Dense"

double eval(const std::string &expression,
            const Evaluator::library_t &library,
            const double x, const Eigen::Vector3d &u)
{
    dynamic_cast<Variable &>(*library.at("x")).value = x;
    dynamic_cast<Variable &>(*library.at("u")).value = u.coeff(0);
    dynamic_cast<Variable &>(*library.at("u'")).value = u.coeff(1);
    dynamic_cast<Variable &>(*library.at("u''")).value = u.coeff(2);
    return Evaluator::eval(expression, library);
}

Eigen::Vector3d estimate(const std::string &expr,
                         const Evaluator::library_t &lib,
                         const double x, const double h, const Eigen::Vector3d &u)
{
    const double k1 = h * u.coeff(1);
    const double l1 = h * u.coeff(2);
    const double m1 = h * eval(expr, lib, x, u);

    const double k2 = h * (u.coeff(1) + l1/2.0);
    const double l2 = h * (u.coeff(2) + m1/2.0);
    const double m2 = h * eval(expr, lib, x + h/2.0, u + Eigen::Vector3d(k1/2.0, l1/2.0, m1/2.0));

    const double k3 = h * (u.coeff(1) + l2/2.0);
    const double l3 = h * (u.coeff(2) + m2/2.0);
    const double m3 = h * eval(expr, lib, x + h/2.0, u + Eigen::Vector3d(k2/2.0, l2/2.0, m2/2.0));

    const double k4 = h * (u.coeff(1) + l3);
    const double l4 = h * (u.coeff(2) + m3);
    const double m4 = h * eval(expr, lib, x + h, u + Eigen::Vector3d(k3, l3, m3));

    return u + Eigen::Vector3d(
            (k1 + 2.0*k2 + 2.0*k3 + k4) / 6.0,
            (l1 + 2.0*l2 + 2.0*l3 + l4) / 6.0,
            (m1 + 2.0*m2 + 2.0*m3 + m4) / 6.0);
}

Eigen::Vector3d shift_left(const Eigen::Vector3d &v)
{
    return Eigen::Vector3d(v.coeff(1), v.coeff(2), 0.0);
}

Eigen::Vector3d integrate(const std::string &expr,
               const Evaluator::library_t &lib,
               const Eigen::Vector3d &u0,
               const double a, const double b, const int num_intervals,
               const std::string *exact_u,
               const std::string *exact_du,
               const std::string *exact_ddu)
{
    const double h = (a - b) / num_intervals;
    double x = b;
    Eigen::Vector3d u = u0;
    std::cerr << "info: (x, u, u', u'')_b = " << x << " " << u.transpose() << std::endl;
    if (exact_ddu == nullptr)
        for (int i = 0; i < num_intervals; ++i) {
            u = estimate(expr, lib, x, h, u);
            x += h;
            std::cerr << "info: (x, u, u', u'') = " << x << " " << u.transpose() << std::endl;
        }
    else
        for (int i = 0; i < num_intervals; ++i) {
            u = estimate(expr, lib, x, h, u);
            x += h;
            const Eigen::Vector3d exact = {eval(*exact_u, lib, x, {}),
                                           eval(*exact_du, lib, x, {}),
                                           eval(*exact_ddu, lib, x, {})};
            std::cerr << "info: (x, u, u', u'') = " << x << " " << u.transpose() << std::endl;
            std::cerr << "info: diff = " << (u - exact).transpose() << std::endl;
        }
    return u;
}

void shoot(const std::string &expr,
           const Evaluator::library_t &lib,
           const double a, const double b, const int num_intervals,
           const double ua, const double ub, const double dub, const Eigen::Vector2d ddub0,
           const int num_max_iterations, const double eps,
           const std::string *exact_u = nullptr,
           const std::string *exact_du = nullptr,
           const std::string *exact_ddu = nullptr)
{
    int num_iterations;
    Eigen::Vector3d ddub(ddub0.coeff(0), ddub0.coeff(1), 0.0);
    std::cerr << "info: with ddub0 and ddub1" << std::endl;
    Eigen::Vector3d rungea(integrate(expr, lib, Eigen::Vector3d(ub, dub, ddub.coeff(0)), a, b, num_intervals, exact_u, exact_du, exact_ddu).coeff(0) - ua,
                           integrate(expr, lib, Eigen::Vector3d(ub, dub, ddub.coeff(1)), a, b, num_intervals, exact_u, exact_du, exact_ddu).coeff(0) - ua,
                           0.0);
    std::cerr << "info: with secant method" << std::endl;
    for (num_iterations = 0; num_iterations < num_max_iterations; ++num_iterations) {
        std::cerr << "---" << std::endl;
        ddub.coeffRef(2) = ddub.coeff(1) - rungea.coeff(1) * (ddub.coeff(1) - ddub.coeff(0))
                                                 / (rungea.coeff(1) - rungea.coeff(0));
        rungea.coeffRef(2) = integrate(expr, lib, Eigen::Vector3d(ub, dub, ddub.coeff(2)), a, b, num_intervals, exact_u, exact_du, exact_ddu).coeff(0) - ua;

        if (std::abs(rungea.coeff(2)) <= eps) {
            std::cerr << "info: final alpha = " << ddub.coeffRef(2) << std::endl;
            ++num_iterations;
            break;
        }

        ddub = shift_left(ddub);
        rungea = shift_left(rungea);
    }
    std::cerr << "info: num_iterations = " << num_iterations << std::endl;
}

int main(int argc, char **argv)
{
    // 1 1 C
    if (argc < 1 + 11 || (argc > 1 + 11 && argc < 1 + 14)) {
        std::cerr << "<u'''(x, u, u', u'')> <a> <b> <N> <u_a> <u_b> <u'_b> <u''_b (alpha_0)> <u''_b (alpha_1)> <K> <eps> [<u(x)> <u'(x)> <u''(x)>]" << std::endl;
        return 0;
    }

    Evaluator::library_t library = {
        {"x", std::shared_ptr<Evaluable>(new Variable())},
        {"u", std::shared_ptr<Evaluable>(new Variable())},
        {"u'", std::shared_ptr<Evaluable>(new Variable())},
        {"u''", std::shared_ptr<Evaluable>(new Variable())}
    };
    library.insert(Evaluator::default_library.begin(),
                   Evaluator::default_library.end());

    const std::string expression(argv[1]);
    const double a = std::strtod(argv[2], nullptr);
    const double b = std::strtod(argv[3], nullptr);
    const int N = std::strtol(argv[4], nullptr, 10);
    const double ua = std::strtod(argv[5], nullptr);
    const double ub = std::strtod(argv[6], nullptr);
    const double dub = std::strtod(argv[7], nullptr);
    const double ddub_0 = std::strtod(argv[8], nullptr);
    const double ddub_1 = std::strtod(argv[9], nullptr);
    const int K = std::strtol(argv[10], nullptr, 10);
    const double eps = std::strtod(argv[11], nullptr);

    if (argc == 1 + 11) {
        shoot(expression, library, a, b, N, ua, ub, dub, {ddub_0, ddub_1}, K, eps);
    } else {
        const std::string exact_u(argv[12]);
        const std::string exact_du(argv[13]);
        const std::string exact_ddu(argv[14]);
        shoot(expression, library, a, b, N, ua, ub, dub, {ddub_0, ddub_1}, K, eps,
              &exact_u, &exact_du, &exact_ddu);
    }
    // u''' = x^3, [0, 1], u = x^6/120, u' = x^5/20, u'' = x^4/4
    // u''' = -u' * u'', [0, 1], u = ln(1 + x), u' = 1/(x+1), u'' = -1/(x^2 + 2*x + 1), u''' = 2/(x^3 + 3*x^2 + 3*x + 1).
    // shoot("0 - 2 * u' * u''", library,
    //       0.0, 1.0, 10,
    //       0.0, 0.6931471805, 0.5, {-0.4, -0.38},
    //       20, 1e-03);
}
