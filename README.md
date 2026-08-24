Main resource: https://repository.lsu.edu/cgi/viewcontent.cgi?article=6692&context=gradschool_theses

2D Incompressible Navier-Stokes Fluid Simulation in Python and C++
- There are two versions, one in pure Python and one in C++. Some parts of the pure Python version were rewritten in C++, thus two different versions exist.
- C++ version is a work in progress

Description: 
- Utilized finite differences, a pressure projection method, and real-time velocity field visualization.

Improvements: 
- Vectorized several nested 'for' loops, improving performance.
- Changed from Gauss-Seidel loop to a Jacobi iteration, improving performance though it will converge more slowly.

