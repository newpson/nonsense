#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <sstream>
#include <vector>
#include <set>

#include "evaluator.h"
#include "gnuplot-helpers.h"
#include "math.h"
#include "variable.h"

#include "Eigen/Dense"
using namespace Eigen;

class CompareVector2i
{
public:
    bool operator()(const Vector2i &a, const Vector2i &b) const
    {
        return a.x() < b.x() || a.y() < b.y();
    }
};

using vecset = std::set<Vector2i, CompareVector2i>;

double eval(
    const std::string &expr,
    const Evaluator::library_t &lib,
    const double t)
{
    dynamic_cast<Variable &>(*lib.at("t")).value = t;
    return Evaluator::eval(expr, lib);
}

std::vector<Vector2d> subdivide(const std::string &expr_x,
                                const std::string &expr_y,
                                const Evaluator::library_t &lib,
                                const bool close,
                                const double ta = 0.0, const double tb = 1.0, const double h = 0.01)
{
    std::vector<Vector2d> points;
    for (double t = ta; t <= tb; t += h) {
        points.push_back(Vector2d(eval(expr_x, lib, t), eval(expr_y, lib, t)));
    }
    if (close)
        points.push_back(points.front());
    return points;
}

std::vector<Vector2d> read_polygon()
{
    std::vector<Vector2d> points;
    std::string line;
    double x = 0.0;
    double y = 0.0;

    while (std::getline(std::cin, line)) {
        if (line == "e") {
            break;
        } else {
            std::stringstream stream(line);
            stream >> x >> y;
            points.push_back(Vector2d(x, y));
        }
    }
    return points;
}

int winding_number(const Vector2d &target, const std::vector<Vector2d> &frame)
{
    // modification of David G. Alciatore algorithm
    // (number of winds is split into halves and wholes)
    int w = 0;
    int hw = 0;
    for (int i = 0; i < frame.size() - 1; ++i) {
        if (frame[i] == target)
            return 1; // TODO (may be not TODO) modify final sum with this value somehow

        if (frame[i].x()>= target.x()) {
            if (frame[i].y()< target.y()) {
                if (frame[i+1].y()> target.y())
                    ++w;
                else if (frame[i+1].y()== target.y())
                    ++hw;
            } else if (frame[i].y()> target.y()) {
                if (frame[i+1].y()< target.y())
                    --w;
                else if (frame[i+1].y()== target.y())
                    --hw;
            } else if (frame[i].y()== target.y()) {
                if (frame[i+1].y()> target.y())
                    ++hw;
                else if (frame[i+1].y()< target.y())
                    --hw;
            }
        }
    }

    return w + hw/2 + hw%2;
}

void plot_pixel(const Vector2i &p, const Vector2d &h)
{
    std::cout << "set object rect from "
              << p.x() * h.x() << "," << p.y() * h.y()
              << " to "
              << (p.x() + 1.0) * h.x() << "," << (p.y() + 1.0) * h.y()
              << std::endl;
    // std::cerr << p.x*h.x << " " << p.y*h.y << std::endl;
}

Vector2d operator/(const Vector2d &a, const Vector2d &b)
{
    return Vector2d(a.x()/b.x(), a.y()/b.y());
}

Vector2i floorize(const Vector2d &a)
{
    return Vector2i(floor_eps(a.x()), floor_eps(a.y()));
}

void plot_pixels(const Vector2d &p, const Vector2d &h, vecset &pixels)
{
    const Vector2d p_grid = p / h;
    const Vector2i p_cell = floorize(p_grid);
    pixels.insert(p_cell);
    if (are_equal(p_grid.x(), p_cell.x()))
        pixels.insert(Vector2i(p_cell.x() - 1, p_cell.y()));
    if (are_equal(p_grid.y(), p_cell.y()))
        pixels.insert(Vector2i(p_cell.x(), p_cell.y() - 1));
    if (floorize(p_grid) == p_cell)
        pixels.insert(Vector2i(p_cell.x() - 1, p_cell.y() - 1));
}

void plot_line(const Vector2d &A, const Vector2d &B, const Vector2d &h, vecset &pixels)
{
    const Vector2d *begin = &A;
    const Vector2d *end = &B;
    if (B.x() < A.x() || (are_equal(A.x(), B.x()) && B.y() < A.y())) {
        begin = &B;
        end = &A;
    }

    plot_pixels(*begin, h, pixels);
    plot_pixels(*end, h, pixels);

    // implementation of John Amanatides & Andrew Woo ray casting algorithm
    // terms: u + v*t  <=>  origin + guide * t
    Vector2d origin = *begin;
    const Vector2d guide = *end - *begin;
    const int dy = guide.y() < 0 ? -1 : 1;

    const Vector2i end_cell = floorize(*end / h);

    int i = 0;
    while (floorize(origin / h) != end_cell && i < 50) {
        const Vector2i origin_cell = floorize(origin / h);
        const Vector2d origin_aligned = Vector2d(origin_cell.x() * h.x(), origin_cell.y() * h.y());
        const Vector2d grid_node = origin_aligned + Vector2d(h.x(), dy * h.y());
        const Vector2d distance_s = (grid_node - origin) / guide;
        const Vector2d distance_u(std::abs(distance_s.x()), std::abs(distance_s.y()));
        if (distance_u.x() < distance_u.y())
            origin = {grid_node.x(), origin.y() + distance_u.x() * guide.y()};
        else if (distance_u.x() > distance_u.y())
            origin = {origin.x() + distance_u.y() * guide.x(), grid_node.y()};
        else
            origin = grid_node;
        plot_pixels(origin, h, pixels);
        ++i;
    }

}

int main(int argc, char **argv)
{
    std::shared_ptr<Variable> t(new Variable());
    Evaluator::library_t library = {
        {"t", t},
    };
    library.insert(Evaluator::default_library.begin(),
                   Evaluator::default_library.end());

    Vector2d step_space(0.023, 0.067);
    std::vector<Vector2d> path;
    std::vector<Vector2d> frame;

    std::string expr_path_x;
    std::string expr_path_y;
    double path_ta = 0.0;
    double path_tb = 1.0;
    double path_tau = 0.01;
    std::string expr_frame_x;
    std::string expr_frame_y;
    double frame_ta = 0.0;
    double frame_tb = 1.0;
    double frame_tau = 0.01;
    for (int i = 1; i < argc; ++i) {
        const std::string par(argv[i]);
        if (par == "--path") {
            expr_path_x = argv[i+1];
            expr_path_y = argv[i+2];
            path_ta = std::strtod(argv[i+3], nullptr);
            path_tb = std::strtod(argv[i+4], nullptr);
            path_tau = std::strtod(argv[i+5], nullptr);
        }
        else if (par == "--frame") {
            expr_frame_x = argv[i+1];
            expr_frame_y = argv[i+2];
            frame_ta = std::strtod(argv[i+3], nullptr);
            frame_tb = std::strtod(argv[i+4], nullptr);
            frame_tau = std::strtod(argv[i+5], nullptr);
        }
        else if (par == "--cell") {
            step_space.x() = std::strtod(argv[i+1], nullptr);
            step_space.y() = std::strtod(argv[i+2], nullptr);
        }
    }

    path = expr_path_x.empty() ? read_polygon() : subdivide(expr_path_x, expr_path_y, library, false, path_ta, path_tb, path_tau);
    frame = expr_frame_x.empty() ? read_polygon() : subdivide(expr_frame_x, expr_frame_y, library, false, frame_ta, frame_tb, frame_tau);

    vecset pixels;
    std::cout << "unset object" << std::endl;
    for (int i = 0; i < path.size() - 1; ++i) {
        const bool is_inside = (winding_number(path[i], frame) != 0 &&
                                winding_number(path[i+1], frame) != 0);
        if (is_inside)
            plot_line(path[i], path[i+1], step_space, pixels);
        else {
            std::cerr << "Out of bounds at " << path[i] << " -> " << path[i+1] << std::endl;
            break;
        }
    }

    std::cout << "$path << e" << std::endl;
    for (int i = 0; i < path.size(); ++i)
        std::cout << path[i].x() << " " << path[i].y() << std::endl;
    std::cout << "e" << std::endl;
    std::cout << "$frame << e" << std::endl;
    for (int i = 0; i < frame.size(); ++i)
        std::cout << frame[i].x() << " " << frame[i].y() << std::endl;
    std::cout << "e" << std::endl;
    std::cout << "plot $path with lines ls 1, $frame with lines ls 2" << std::endl;

    for (const Vector2i &pixel: pixels) {
        plot_pixel(pixel, step_space);
    }

    return 0;
}
