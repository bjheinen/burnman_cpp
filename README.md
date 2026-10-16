# burnman_cpp

C++17 translation of core functionality in [BurnMan](https://github.com/geodynamics/burnman).

The aim is a faster core so burnman can be used for large-scale inversions.

Tests can be found in tests/ with benchmarks alongside them (tests/benchmarks/)

## Requirements

| Dependency | Version | Notes |
|---|---|---|
| C++ compiler | C++17 | Built and tested with GCC 13 & 16 |
| [Eigen](https://eigen.tuxfamily.org) | **3.4.x (not 5.x)** | Eigen 5.0 introduces breaking changes |
| [GSL](https://www.gnu.org/software/gsl/) | 2.x | Tested with 2.7.1 |
| [Catch2](https://github.com/catchorg/Catch2) | ≥ 3.4 | Tests and benchmarks only |

The Makefile looks for headers in `/usr/local/include` and `/usr/include/eigen3`,
and libraries in `/usr/local/lib`. If yours are elsewhere, edit `EXTRA_INCLUDE`
and `EXTRA_LIB` at the top of the Makefile.

## Building

```bash
make                               # builds the static library lib/libburnman.a
make test                          # builds the test runner bin/run_tests
./bin/run_tests                    # runs tests
```

To use the library in your own code, add `include/` to the include path and
link against it along with GSL:

```bash
g++ -std=c++17 -Iinclude -I/usr/include/eigen3 code.cpp \
    -Llib -lburnman -lgsl -lgslcblas -lm
```

Build modes are selected with `BUILD_MODE` (default `release`, `-O3`) and apply
to both the library and the test runner:

```bash
make BUILD_MODE=debug        # -g -Og             -> lib/libburnman_debug.a
make test BUILD_MODE=test    # ASan + UBSan       -> bin/run_tests_test
make test PROFILE=1          # adds -pg for gprof -> bin/run_tests_profile
```

`make clean-mode` (with the same flags) removes one mode's build and `make clean` removes everything.

## Benchmarks

Benchmarks use Catch2's benchmarking support and are tagged `[!benchmark]`:

```bash
./bin/run_tests [!benchmark]
```

To save results as XML:

```bash
./bin/run_tests [!benchmark] \
  --reporter XML::out=./benchmark-report.xml \
  --reporter console::out=-::colour-mode=ansi
```

## Licence

GPL-3.0-or-later. See [LICENSE](LICENSE).

burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
