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
    bool going_through = false;
    int w = 0;
    for (int i = 0; i < frame.size() - 1; ++i) {
        if (frame[i].x >= target.x) {
            // if (frame[i].y < target.y && target.y <= frame[i+1].y) {
            //     ++w;
            // }
            // else if (frame[i].y > target.y && target.y >= frame[i+1].y)
            //     --w;
            if (frame[i].y < target.y && target.y < frame[i+1].y) {
                ++w;
            } else if (frame[i].y > target.y && target.y > frame[i+1].y) {
                --w;
            } else if (frame[i].y == target.y && target.y < frame[i+1].y) {
                if (!going_through)
                    ++w;
                going_through = false;
            } else if (frame[i].y == target.y && target.y > frame[i+1].y) {
                if (!going_through)
                    --w;
                going_through = false;
            } else if (frame[i].y < target.y && target.y == frame[i+1].y) {
                ++w;
                going_through = true;
            } else if (frame[i].y > target.y && target.y == frame[i+1].y) {
                --w;
                going_through = true;
            }
        }
        std::cerr << "w = " << w << std::endl;
    }
    return w;
}

void plot_line(const Vector2D &begin, const Vector2D &end, const Vector2D &h)
{
    const Vector2D linevector = end - begin;
    const Vector2D direction(std::copysign(1.0, linevector.x),
                             std::copysign(1.0, linevector.y));
    const double invertor = copysign(1.0, linevector.x * linevector.y); // инверсия осей в чётных квадрантах

    const Vector2D cell_begin(direction.x > 0 ? std::floor(begin.x/h.x)*h.x : std::ceil(begin.x/h.x)*h.x,
                              direction.y > 0 ? std::floor(begin.y/h.y)*h.y : std::ceil(begin.y/h.y)*h.y);
    const Vector2D cell_end(direction.x > 0 ? std::floor(end.x/h.x)*h.x : std::ceil(end.x/h.x)*h.x,
                            direction.y > 0 ? std::floor(end.y/h.y)*h.y : std::ceil(end.y/h.y)*h.y);

    int i = 0;
    Vector2D cell = cell_begin;
    while ((std::abs(cell.x - cell_end.x) > h.x/100.0
            || std::abs(cell.y - cell_end.y) > h.y/100.0) && i < 1000) {
        std::cout << "set object rect from "
                  << cell.x << "," << cell.y
                  << " to "
                  << cell.x + direction.x*h.x << "," << cell.y + direction.y*h.y
                  << std::endl;
        std::cerr << "Going through " << cell << std::endl;
        cell.x += h.x*direction.x;
        cell.y += h.y*direction.y;
        const double D = invertor*((cell.y - begin.y)*linevector.x - (cell.x - begin.x)*linevector.y);
        if (D > 0)
            cell.y -= h.y*direction.y;
        else if (D < 0)
            cell.x -= h.x*direction.x;
        ++i;
    }
    if (i >= 1000)
        std::cerr << "Sampling overflow at " << begin << " -> " << end << std::endl;

    std::cout << "set object rect from "
              << cell.x << "," << cell.y
              << " to "
              << cell.x + direction.x*h.x << "," << cell.y + direction.y*h.y
              << std::endl;
}

int main(int argc, char **argv)
{
    std::shared_ptr<Variable> t(new Variable());
    Evaluator::library_t library = {
        {"t", t},
    };
    library.insert(Evaluator::default_library.begin(),
                   Evaluator::default_library.end());

    Vector2D step_space(0.1, 0.1);
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
