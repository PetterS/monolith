# Column Generation

Implements the algorithm described in [_First-order Linear Programming in a Column Generation-Based Heuristic Approach to the Nurse Rostering Problem_](https://www.strandmark.net/papers/first-order-scheduling.pdf) (2020) [doi link](https://doi.org/10.1016/j.cor.2020.104945).

## Running the Experiments in the Paper
The software works on Windows, Linux and MacOS. All dependencies are included in this repo, so it should be fairly easy to build. The latest versions of Visual Studio, Gcc or Clang are strongly recommended.

- Set up your build environment with CMake.
- _Optional_: Run unit tests with `ctest -R colgen`.
- Build the target `shift_scheduling_colgen`.
- Run
  ```
  shift_scheduling_colgen Instance1.txt --num_solutions=5 --use_first_order_solver
  ```
  This should finish in less than a second. If it worked, feel free to continue on the large instances. The `--help` flag will give you a list of options.
  
This software can also be [run online in the browser](https://www.strandmark.net/wasm/shift_scheduling_colgen_page.html).

