/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#define CATCH_CONFIG_MAIN // Don't need if linking against libCatch2Main?
// #define CATCH_CONFIG_ENABLE_BENCHMARKING
#include <catch2/catch_all.hpp>


// g++ -std=c++17 -Iinclude -Itests/include -I/usr/local/include -L/usr/local/lib tests/test_main.cpp tests/eos/test_debye.cpp src/eos/debye.cpp -lgsl -lgslcblas -lm -lCatch2Main -lCatch2 -o test

// Macros:
// REQUIRE ( expression )
// CHECK ( expression ) --> check is like require, but will not stop all test
// e.g. 
// CHECK( thisReturnsTrue() );
// REQUIRE( i == 42 );

// for false, do
// REQUIRE_FALSE( expression ) and
// CHECK_FALSE( expression )

// To check exception is called:
// REQUIRE_THROWS( expression ) and
// CHECK_THROWS( expression )

// or specific
// REQUIRE_THROWS_AS( expression, exception type ) and
// CHECK_THROWS_AS( expression, exception type )

// NOTHROW to check no exception thrown

// For floating point:
// #include <catch2/matchers/catch_matchers_floating_point.hpp>

// WithinAbs(double target, double margin)
// WithinRel(FloatingPoint target, FloatingPoint eps)
// WithinULP(FloatingPoint targert, unit64_t maxUlpDiff)

// WithinRel
// REQUIRE_THAT(1.0, WithinRel(1.1, 0.1));
// REQUIRE_THAT(1.1, WithinRel(1.0, 0.1));

// withinULP
// REQUIRE_THAT( -0.f, WithinULP( 0.f, 0 ) );
