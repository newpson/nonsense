#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <sstream>
#include <vector>

#include "evaluator.h"
#include "gnuplot-helpers.h"
#include "math.h"
#include "variable.h"

#include "Eigen/Dense"

double eval(
    const std::string &expr,
    const Evaluator::library_t &lib,
    const double t)
{
    dynamic_cast<Variable &>(*lib.at("t")).value = t;
    return Evaluator::eval(expr, lib);
}

std::vector<Vector2D> subdivide(const std::string &expr_x,
                                const std::string &expr_y,
                                const Evaluator::library_t &lib,
                                const bool close,
                                const double ta = 0.0, const double tb = 1.0, const double h = 0.01)
{
    std::vector<Vector2D> points;
    for (double t = ta; t <= tb; t += h) {
        points.push_back(Vector2D(eval(expr_x, lib, t), eval(expr_y, lib, t)));
    }
    if (close)
        points.push_back(points.front());
    return points;
}

std::vector<Vector2D> read_polygon()
{
    std::vector<Vector2D> points;
    std::string line;
    double x = 0.0;
    double y = 0.0;

    while (std::getline(std::cin, line)) {
        if (line == "e") {
            break;
        } else {
            std::stringstream stream(line);
            stream >> x >> y;
            points.push_back(Vector2D(x, y));
        }
    }
    return points;
}

int winding_number(const Vector2D &target, const std::vector<Vector2D> &frame)
{
    // modification of David G. Alciatore algorithm
    // (number of winds is split into halves and wholes)
    int w = 0;
    int hw = 0;
    for (int i = 0; i < frame.size() - 1; ++i) {
        if (frame[i] == target)
            return 1; // TODO (may be not TODO) modify final sum with this value somehow

        if (frame[i].x >= target.x) {
            if (frame[i].y < target.y) {
                if (frame[i+1].y > target.y)
                    ++w;
                else if (frame[i+1].y == target.y)
                    ++hw;
            } else if (frame[i].y > target.y) {
                if (frame[i+1].y < target.y)
                    --w;
                else if (frame[i+1].y == target.y)
                    --hw;
            } else if (frame[i].y == target.y) {
                if (frame[i+1].y > target.y)
                    ++hw;
                else if (frame[i+1].y < target.y)
                    --hw;
            }
        }
    }

    return w + hw/2 + hw%2;
}

void plot_square(const Vector2D &p, const Vector2D &h)
{
    std::cout << "set object rect from "
              << p.x*h.x << "," << p.y*h.y
              << " to "
              << (p.x + 1.0)*h.x << "," << (p.y + 1.0)*h.y
              << std::endl;
    // std::cerr << p.x*h.x << " " << p.y*h.y << std::endl;
}

void plot_point(const Vector2D &p, const Vector2D &h)
{
    const Vector2D p_grid = p / h;
    const Vector2D p_cell = p_grid.floor();
    std::cerr << "origin: " << p << std::endl;
    std::cerr << "guarantee: " << p_cell * h << std::endl;
    plot_square(p_cell, h);
    if (are_equal(p_grid.x, p_cell.x))
        plot_square(Vector2D(p_cell.x - 1.0, p_cell.y), h);
    if (are_equal(p_grid.y, p_cell.y))
        plot_square(Vector2D(p_cell.x, p_cell.y - 1.0), h);
    if (p_grid == p_cell)
        plot_square(p_cell - Vector2D(1.0, 1.0), h);
}

void plot_line(const Vector2D &A, const Vector2D &B, const Vector2D &h)
{
    const Vector2D *begin = &A;
    const Vector2D *end = &B;
    if (B.x < A.x || (are_equal(A.x, B.x) && B.y < A.y)) {
        begin = &B;
        end = &A;
    }

    plot_point(*begin, h);
    plot_point(*end, h);

    // implementation of John Amanatides & Andrew Woo ray casting algorithm
    // terms: u + v*t  <=>  origin + guide * t
    Vector2D origin = *begin;
    const Vector2D guide = *end - *begin;
    const int dy = guide.y < 0 ? -1 : 1;

    const Vector2D end_cell = end->floor(h);

    int i = 0;
    while (origin.floor(h) != end_cell && i < 50) {
        const Vector2D origin_aligned = origin.floor(h) * h;
        const Vector2D grid_node = origin_aligned + Vector2D(h.x, dy * h.y);
        const Vector2D distance_s = (grid_node - origin) / guide;
        const Vector2D distance_u(std::abs(distance_s.x), std::abs(distance_s.y));
        if (distance_u.x < distance_u.y)
            origin = {grid_node.x, origin.y + distance_u.x * guide.y};
        else if (distance_u.x > distance_u.y)
            origin = {origin.x + distance_u.y * guide.x, grid_node.y};
        else
            origin = grid_node;
        plot_point(origin, h);
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

    Vector2D step_space(0.023, 0.067);
    std::vector<Vector2D> path;
    std::vector<Vector2D> frame;

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
            step_space.x = std::strtod(argv[i+1], nullptr);
            step_space.y = std::strtod(argv[i+2], nullptr);
        }
    }

    path = expr_path_x.empty() ? read_polygon() : subdivide(expr_path_x, expr_path_y, library, false, path_ta, path_tb, path_tau);
    frame = expr_frame_x.empty() ? read_polygon() : subdivide(expr_frame_x, expr_frame_y, library, false, frame_ta, frame_tb, frame_tau);

    std::cout << "unset object" << std::endl;
    for (int i = 0; i < path.size() - 1; ++i) {
        const bool is_inside = (winding_number(path[i], frame) != 0 &&
                                winding_number(path[i+1], frame) != 0);
        if (is_inside)
            plot_line(path[i], path[i+1], step_space);
        else {
            std::cerr << "Out of bounds at " << path[i] << " -> " << path[i+1] << std::endl;
            break;
        }
    }

    std::cout << "$path << e" << std::endl;
    for (int i = 0; i < path.size(); ++i)
        std::cout << path[i] << std::endl;
    std::cout << "e" << std::endl;
    std::cout << "$frame << e" << std::endl;
    for (int i = 0; i < frame.size(); ++i)
        std::cout << frame[i] << std::endl;
    std::cout << "e" << std::endl;
    std::cout << "plot $path with lines ls 1, $frame with lines ls 2" << std::endl;

    return 0;
}
