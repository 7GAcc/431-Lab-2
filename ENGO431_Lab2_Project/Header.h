#ifndef FUNCTIONS_CPP
#define FUNCTIONS_CPP
#define EIGEN_NO_DEBUG
#include "Eigen/Eigen"
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
using namespace std;
using namespace Eigen;
using Eigen::MatrixXd;
using Eigen::VectorXd;


//Reads a text file and converts it into a Matrix, given the file name and number of columns
MatrixXd Matrix_readIn(int col, string filename);
//Converts a matrix in [dms] to a vector in [decimal degrees]
VectorXd DMS_DD(MatrixXd m);
//Outputs a Matrix into a text file, given the name. Overloaded to include units
void MatToText(MatrixXd M, string filename);
void MatToText(MatrixXd M, string filename, vector<string> units);
//Another overload to put units on the first row
void MatToText(MatrixXd M, string filename, vector<string> units, char a);
//Converts a string array into a text file
void VecToText(vector<string> s, string filename);
//Calculates design matrix from fiducial mark comparator coordinates.
//Input matrix must have at least 3 columns, where col 1=fiducial point number, 2= x-coordinate, 3= y-coordinate
MatrixXd A_mat_calc(MatrixXd fid_c_coords);
//Calculates xHat from design matrix, weight matrix, and observations vector
VectorXd xHat_Calc(MatrixXd A, MatrixXd P, VectorXd l);
VectorXd xHat_Calc(MatrixXd A, VectorXd l);
//Calculates non-linear parameters sx,sy,delta,theta from linear parameters a,b,c,d
VectorXd xHat_nonlin_Calc(VectorXd xHat_lin);
//Applies the affine transformation to a set of comparator coordinates using precalculated affine parameters
MatrixXd Affine_Transformation(MatrixXd abcd, VectorXd xf, VectorXd yf, VectorXd translation);
//Calculates principle point, radial lens, decentering lens, and atmospheric refraction corrections using values from the certificate, the heights above the datum, and flying height
MatrixXd Corrections_Calculator(VectorXd x, VectorXd y, double xp, double yp, double k0, double k1, double k2, double k3, double p1, double p2, double c, double H, double h);

#endif
