# Monolith
Monolith is a monorepo with several optimization projects. Some of the code was originally written for research or as hobby projects in other repositories (e.g. [spii](https://github.com/PetterS/spii) and [easy-IP](https://github.com/PetterS/easy-IP)).

One of the highlights is a state-of-the-art scheduler using column generation, which significantly outperforms all other optimizers at [schedulingbenchmarks.org](http://www.schedulingbenchmarks.org/). **Try it in the browser (wasm) [here](https://www.strandmark.net/wasm/shift_scheduling_colgen_page.html)!**


## Why a monorepo?
 - C++ does not have an ABI. Every compiler, or worse, every flag configuration of every compiler generates potentially incompatible code. I want to use many compilers (MCVC, GCC, Clang) and many settings (debug, release, asan, fuzzers etc.). I also use Emscripten to compile programs to WASM (example [here](https://www.strandmark.net/wasm/glpk.html)).
- Refactoring code becomes much easier if all code with all dependencies is available in one IDE at the same time.
 
## Contact
Petter Strandmark, petter.strandmark@gmail.com

## Modules

### Column generation
The `minimum::linear::colgen` module contains code for solving scheduling problems. It significantly outperforms all other optimizers at [schedulingbenchmarks.org](http://www.schedulingbenchmarks.org/).

Some of the reasons it is fast:
- It uses a first-order LP solver based on papers by Chambolle and Pock.
- The [Ryan-Foster rule](https://strandmark.wordpress.com/2018/01/24/visualizing-the-ryan-foster-rule/) is used to iteratively work towards an integer solution. There is no time to branch and bound for big problems.
- The pricing problem uses highly optimized dynamic programming in a DAG (in `minimum::algorithms`).

### Minimum/AI
The `minimum::ai` module contains code originally from my [Monte-Carlo tree search](https://github.com/PetterS/monte-carlo-tree-search) repository. It is very fast – some years ago it evaluated almost 2 million *complete* games per second when playing 4 in a row.

### Minimum/Linear
The `minimum::linear` module contains code originally from my [easy-IP](https://github.com/PetterS/easy-IP) repository. It is a modelling interface (DSL) for integer programming and supports converting IPs to SAT. It also uses two free solvers: Cbc and Glpk.

`minimum::linear` also contains an implementation of a first-order LP solver based on papers by Chambolle and Pock. Using a first-order LP solver is very important when solving really big scheduling problems.

### Minimum/Nonlinear
The `minimum::nonlinear` module contains code originally from my [spii](https://github.com/PetterS/spii) repository. When a function is defined, it will be automatically differentiated with an autodiff library. This makes the resulting code very fast and numerically stable. The library contains minimization routines like L-BFGS and Newton’s method.

For example, a function of one vector variable is simply defined as:
```
auto lambda = [](auto* x) {
    auto d0 = x[1] - x[0] * x[0];
    auto d1 = 1 - x[0];
    return 100 * d0 * d0 + d1 * d1;
};
```
Having `auto` here instead of e.g. `double` is crucial. Efficient code for first and second-order derivatives is generated with

```
auto term = make_differentiable<2>(lambda);
```

### Minimum/Constrained
Using both `minimum::linear` and `minimum::nonlinear`, this module implements constrained optimization via successive linear programming.

### Minimum/Curvature
The module `minimum::curvature` contains code originally written at [curve_extraction](https://github.com/PetterS/curve_extraction). It computes shortest paths in regular graphs taking curvature into account.

## Building
Everything needed to build the library is included. With CMake and a modern C++ compiler, you can do

    $ cmake /path/to/source
    $ make

## Licence
The license can be found in the `LICENSE` file. Contact me if you need another license.

The software in the `third-party` folder is not written exclusively by me and their respective license files apply.
