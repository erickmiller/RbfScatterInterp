# RbfScatterInterp

Author: Erick Miller

Purpose: Radial Basis Function based interpolation in unlimited multi-dimensional vector space 

Dependencies:  None. Standard C++ libraries. Algo and math are all self contained.

This is a self contained C++ library for using Radial Basis Function based N-dimensional scattered data interpolation for interpolating between unlimited multi-dimensional vectors in ND space (solver uses Gaussian elimination with partial pivoting, which I found was more stable and smother than thin plate splines which caused more overshooting and tangent space issues)

Using it is easy.  Just include it and initialize the RbfScatterInterp (look at the three versions of the constructor):

All versions of the constructor require you pass in two arguments: 
1) **&pts** Is a reference to a vector of vectors of doubles (these are your ND vectors)
```cpp
vector<vector<double> > &pts
```
The first argument, pts, is your mulit-dimensional vector space.  The number of dimensionss of your vectors space is determined by the length of each vector conained within the vector of vectors.   The number of vectors (or the length of the outside vector) represents how many data points you will be using in the interpolate and thus also defines the spacial stability of your space.  If you imagine wanting to interpolate in stable manner in a 3D space, then you would want to have enough 3D vectors to represent a stable space, thus you'd want at a minimum the equivalent of the Cartesian axis to be defined which would simply be 3 vectors of X, Y and Z.  For an even more stable interpolation, perhaps you'd want to define the space using each row/column (depending on transposition) of a transformation matrix as the vector of your space, so you could capture rotation and scale somewhat more accurately during interpolation.

2) **&values** Is a reference to a vector of doubles
```cpp
vector<double> &values
```
The second argument represents the a vector values, that must have the same number of data points (length) as pts.  Values is simply an array of one dimensional values that will feed the solver along with pts upon initialization, and will then be represented as the variable input for the interpolator to generator interpolated output as the ND space is interpolated through. 

The simple interface easiest to use for initialization:

```cpp                        
        RbfScatterInterp(vector<vector<double> > &pts,  
                                    vector<double> &values );
```

Two alternate versions of the constructor allow slightly more nuanced control by exposing a sigma multiplier as well as a value that's used during initialization to help add stability to the sovler through a matrix regularization hack.  Here are the other two interfaces: 

```cpp
        RbfScatterInterp( vector<vector<double> > &pts,  
            vector<double> &values, double sigma_mult, double reg);
```


```cpp
        RbfScatterInterp(vector<vector<double> > &pts, 
                        vector<double> &values, double sigma_mult);
```

Then, after the class has been created and intialized with the needed numerical data, you can run interpolation of the ND space at any time by passing in a vector of the same length and nature as the values vector.  Typically the data passed to interp will represent some new state or new composition of the values numbers such that passing them in will represent a new point to be interpolated, for example see the interface:

```cpp
double RbfScatterInterp::interp(vector<double> &curr_pts)
```

Ok cool, this code was fun to write and was a bit challenging.   It has been tested very extensively and ran reliably as a core interpolation engine using some pretty large data sets within some relatively critical software in the past, so it should prove to be relatively stable and reliable.  At least any major kinks should have been worked out. This software also has not been tested in years so things may have changed.  Use at your own peril :)  Perhaps could be optimized for speed using some matrix optimizers or CUDA as well. 

Background:
This code was largely inspired by my friend JPLewis' (and Matt) work on n-d pose space, you can find more info at http://scribblethink.org/Work/PSD/
The simple 2D java code was used heavily as reference when writing this, enough so that I think it's fair to call this open source derivative work, despite this library being far more sophisticated in that it is extended into ND multi dimensional space and has a additional couple novel extensions.  I actually wrote this code in 2004, and am open sourcing it now due to me needing to dig it up and use it for some interpolation experiments, and then realizing the original code was released under GNU 2 License, which dictates I also open source this derivative work, so yay there you go :) 
