C++ translation of burnman (https://github.com/geodynamics/burnman)

Requires: C++17, Eigen >= 3.4, GSL
Tests/Benchmarks require: Catch2 >= 3.4

Install: make
Build tests/benchmarks: make test

Run benchmarks with: ./bin/run_tests [!benchmark]
Save output to xlm also with:
./bin/run_tests [!benchmark] --reporter XML::out=./benchmark-report.xml --reporter console::out=-::colour-mode=ansi
Use parse_benchmarks.py to process xml output

