# NumToolLibs

This project is a collection of typical functions, types and/or classes in economic & financial research. These methods are developed in several languages (C, Julia, C++, etc.) in which APIs will be unified regularly (but may be not frequently).

We expect to keep this project out of published/compiled libraries until they work stably.

Collaborators: [Zeqi Jin](https://github.com/Jinzeqi) (Fordham U)

## Sources

Nowadays (Spring 2019), the project is expected to cover:

1. **StochInt**: a collection of methods for stochastic integral, control models and stochastic process; including simulation (e.g. BM, GBM, OU, Levy).
2. **NovDist**: a collection of newly proposed distributions & copula, e.g. B, MEGB2; covers type declaration, generating function, PDF/CDF, common moments, estimation, sampling. As an important to current similar libraries.
3. **EasyMetrics**: a collection of decorating functions on existing popular statistical/econometric libraries (e.g. GLM in Julia). This library works to provide more economist-style APIs & result display for scholars.
4. **ConvFuncs**: a collection of conventional "pure" functions in research, e.g. CD function, CRRA, CES, annuity, insurance pricing.
5. **EasyIO**: a collection of input-output methods which significantly simplifies the I/O process. e.g. time tagging, dictionary to a multi-sheet excel file.

## Requirements

We follow some important principles to make our code readable, high-performance and extendable:

1. **Doc Doc Doc**: 

   1. every library is required to have a single documentation to state at least: module name, purpose, member list, brief member introduction, references and last update date. 
   2. Every function/method should have a single documentation to state: what this function to do, input (name, type, definition), output (name, type, definition), algorithm, essential mathematics, references. 
   3. Every class/type should have a single documentation to state: what this class to do, inputs in initialization, essential mathematics, references.
   4. The documentation should be at the beginning of the source file/function/class and should be visible to IDE's documentation manager if applicable.
   5. If applicable, please indicate which book & chapter to read for a specific function's algorithm. For original algorithm, if hard to explain the mathematics in function documentations, please provide a single academic documentation (e.g. a pdf or md file) with your source file, and mark it in the function documentation.

2. **Never too many comments!**:  Because these libraries are not generic-purpose but for academic research, they will become f-word black sh\*t boxes if without detailed comments. We highly recommend contributors to explain ***every*** sentence they write, even it is only a type conversion. In many cases, a simple sentence "pulls one hair and move the whole body" (牵一发而动全身) e.g.

   ```julia
   # julia
   ColIdx = LinRange( 0.0, 1.0, 3 )
   ColIdx = Array{Float64}( ColIdx ) # convert ColIdx to an Array because we then operate on element 
   ColIdx[2] = -1
   ```

3. **Less dependency & More core types** :
   1. We highly suggest contributors to use/import/include as least 3rd-party packages as possible (standard libraries does not matter). This avoids possible incompatibility problems when those packages update. If you do need to import one, please comment your import/using/include sentence to tell readers what this module/package to do in this source file, its version when developing current source file, etc. As many popular tutorials, please avoid to expose 3rd-party package's functions to your current namespace but to use "Pkg.func()".
   2. Meanwhile, we highly suggest contributors to use language-core types/data structures to acquire higher performance, better compatibility and better ability of data exchange.
4. **Static and safe types**: 
   1. We suggest to use static-type languages or use static declaration in dynamic-type languages if applicable. Static types can avoid most mistakes which would not be reported in dynamic-type languages. 
   2. Please try to write type-safe functions: only one type data to return, esp. in pure functions. And, please **explicitly** declare local variables to avoid misunderstanding (esp. when a name like "ColName" is common for every function).
5. **Good naming**: as many popular programming tutorials, we suggest to use meaningful names, e.g. camel for methods, pascal for variables, all-uppercase for macro/constant/global, and all-lowercase for temporary variables such as indices.
6. **Remember to test & demo**: Please provide a "test.xx" file which provides demo callings for every method defined in your source file. Please comment to  explain how to use every function in work flows (e.g. simulate stock prices) and provide performance testing with testing benchmark.



## References

1. Stokey, N. L. . (2008). The economics of inaction: stochastic control models with fixed costs.
2. Press, William H., Flannery, B. P., Teukolsky, S. A., & Vetterling, W. T. (1988). Numerical recipes in C: the art of scientific computing. Numerical Recipes in C The Art of Scientific Computing.