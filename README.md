## Markov Chain Monte Carlo solver for linear differential equations

Hi everyone! This is my C++ implementation of an MCMC algorithm for solving large systems of differnetioal equations. Theory behind it is described in the paper that is uploaded on the repo. 

This version of the solver is by no means universal, and only deals with the system of equations described in part 4 of the paper. 

I tried to construct functions in a way that it would be easy to follow the steps described in the main body of the work; comments are currently lacking, but hopefully will be added later.

## Building


The build has only been tested for Ubuntu 24.04:
+ [CMake](https://cmake.org/download/) `>= 3.20.2`
+ [make](https://www.gnu.org/software/make/) `4.3` and [GCC](https://gcc.gnu.org/) `11.5.0` 

In order to compile:

From the root directory
```bash
# Create the build directory and change into it
mkdir -p build && cd build

# Configure CMake from the command line passing configure-time options. 
cmake ..

#Build target
make 
```
In order to execute:
```bash
./lin_ode_mcmc
```

If you require a debug build, substitute for 
``` bash
# Configure CMake from the command line passing configure-time options. 
cmake .. -DCMAKE_BUILD_TYPE=Debug

# Build target
make 

# Execute
gdb ./lin_ode_mcmc_db
```

## Collaboration
Feel free to fork the repo and play around with it: e.g. you can build different equation systems, etc.

Don't hesitate to contact me if you have any questions: smmxgr@gmail.com

## Credits

Sergey Ermakov

Tamara Surovikina

## License
 
The MIT License (MIT)

Copyright (c) Max Smilovitskiy

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.