Main resource: https://repository.lsu.edu/cgi/viewcontent.cgi?article=6692&context=gradschool_theses

2D Incompressible Navier-Stokes Fluid Simulation in Python and C++
- There are two versions, one in pure Python and one in C++. Some parts of the pure Python version were rewritten in C++, thus two different versions exist.
- For C++ version, generate data by executing the main.cpp file. It will store this data in a binary file called "data.bin". To visualize this data, run plot.py. 

Description: 
- Utilized finite differences, a pressure projection method, and real-time velocity field visualization.

Improvements: 
- Vectorized several nested 'for' loops in Python version, improving performance.

TO-DO:
- Change Jacobi iteration to red-black Gauss-Seidel method.
- Add a heatmap of velocities for both C++ and pure Python version.
