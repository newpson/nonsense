cmake_minimum_required(VERSION 3.14)
project(nonsense CXX)
add_compile_options(-g)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(eval)
add_subdirectory(eigen)

add_executable(app-penalty)
target_sources(app-penalty
    PRIVATE
    penalty.cpp)
target_link_libraries(app-penalty
    PUBLIC
    evaluator
    eigen)

add_executable(app-gradient-descent)
target_sources(app-gradient-descent
    PRIVATE
    gradient-descent.cpp)
target_link_libraries(app-gradient-descent
    PUBLIC
    evaluator
    eigen)

add_executable(app-newton2d)
target_sources(app-newton2d
    PRIVATE
    newton2d.cpp)
target_link_libraries(app-newton2d
    PUBLIC
    evaluator
    eigen)

add_executable(app-newton1d)
target_sources(app-newton1d
    PRIVATE
    newton1d.cpp)
target_link_libraries(app-newton1d
    PUBLIC
    evaluator)

add_executable(app-shooting)
target_sources(app-shooting
    PRIVATE
    shooting.cpp)
target_link_libraries(app-shooting
    PUBLIC
    evaluator)

add_executable(app-runge-kutta)
target_sources(app-runge-kutta
    PRIVATE
    runge-kutta.cpp)
target_link_libraries(app-runge-kutta
    PUBLIC
    evaluator)

add_executable(app-tridiag)
target_sources(app-tridiag
    PRIVATE
    tridiag.cpp)
target_link_libraries(app-tridiag
    PUBLIC
    evaluator)

add_executable(app-polyline)
target_sources(app-polyline
    PRIVATE
    polyline.cpp)
target_link_libraries(app-polyline
    PUBLIC
    evaluator)

add_executable(app-parabola)
target_sources(app-parabola
    PRIVATE
    parabola.cpp)
target_link_libraries(app-parabola
    PUBLIC
    evaluator)

add_executable(app-golden)
target_sources(app-golden
    PRIVATE
    golden.cpp)
target_link_libraries(app-golden
    PUBLIC
    evaluator)

add_executable(app-dichotomy)
target_sources(app-dichotomy
    PRIVATE
    dichotomy.cpp)
target_link_libraries(app-dichotomy
    PUBLIC
    evaluator)

add_executable(app-dichotomy-roots)
target_sources(app-dichotomy-roots
    PRIVATE
    dichotomy-roots.cpp)
target_link_libraries(app-dichotomy-roots
    PUBLIC
    evaluator)
