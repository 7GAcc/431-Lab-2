#include "Header.h"

MatrixXd Matrix_readIn(int col, string filename) {
    ifstream infile;
    MatrixXd result(1, col);
    infile.open(filename, ifstream::in);
    int rows = 1;
    while (!infile.eof()) {
        double a;
        result.conservativeResize(rows, col);
        for (int i = 0; i < col; i++) {
            infile >> a;
            result(rows - 1, i) = a;
        }
        rows++;
    }
    infile.close();
    return result;
}

VectorXd DMS_DD(MatrixXd m) {
    VectorXd result(m.rows());
    for (int i = 0; i < m.rows(); i++) {
        result(i) = m(i, 0) + m(i, 1) / 60 + m(i, 2) / 3600;
    }
    return result;
}

//Outputs Matrix to text file
void MatToText(MatrixXd M, string filename) {
    ofstream outfile;
    outfile.open(filename, ofstream::out); //Output file with specified name created
    double val;
    for (int i = 0; i < M.rows(); i++) {
        for (int j = 0; j < M.cols(); j++) {
            val = M(i, j);
            outfile << val << "\t"; //Stores i-th row and j-th column in text file
        }
        outfile << endl;
    }
    outfile.close();
}
void MatToText(MatrixXd M, string filename, vector<string> units) {
    ofstream outfile;
    outfile.open(filename, ofstream::out); //Output file with specified name created
    double val;
    for (int i = 0; i < M.rows(); i++) {
        for (int j = 0; j < M.cols(); j++) {
            val = M(i, j);
            outfile << val << "\t"; //Stores i-th row and j-th column in text file
        }
        outfile << units[i] << endl;
    }
    outfile.close();
}
void MatToText(MatrixXd M, string filename, vector<string> units, char a) {
    ofstream outfile;
    outfile.open(filename, ofstream::out); //Output file with specified name created
    double val;
    for (int i = 0; i < M.cols(); i++) {
        outfile << units[i] << "\t";
    }
    outfile << endl;
    for (int i = 0; i < M.rows(); i++) {
        if (i != 0)
            outfile << endl;
        for (int j = 0; j < M.cols(); j++) {
            val = M(i, j);
            outfile << fixed << setprecision(9) << val << "\t"; //Stores i-th row and j-th column in text file
        }

    }
    outfile.close();
}
//Converts a string vector to text
void VecToText(vector<string> s, string filename) {
    ofstream outfile;
    outfile.open(filename, ofstream::out); //Output file with specified name created
    for (int i = 0; i < s.size(); i++) {

        outfile << s[i] << endl;
    }
    outfile.close();
}

MatrixXd A_mat_calc(MatrixXd fid_c_coords) {
    MatrixXd result(fid_c_coords.rows() * 2, 6);
    VectorXd x = fid_c_coords.col(1);
    VectorXd y = fid_c_coords.col(2);
    result.fill(0);
    for (int i = 0; i < x.rows(); i++) {
        result(2 * i, 0) = x(i);
        result(2 * i, 1) = y(i);
        result(2 * i, 2) = 1;
        result(2 * i + 1, 3) = x(i);
        result(2 * i + 1, 4) = y(i);
        result(2 * i + 1, 5) = 1;
    }
    return result;
}

VectorXd xHat_Calc(MatrixXd A, MatrixXd P, VectorXd l) {
    VectorXd result(A.cols());
    result = (A.transpose() * P * A).inverse() * A.transpose() * P * l;
    return result;
}
VectorXd xHat_Calc(MatrixXd A, VectorXd l) {
    VectorXd result(A.cols());
    result = (A.transpose() * A).inverse() * A.transpose() * l;
    return result;
}

VectorXd xHat_nonlin_Calc(VectorXd xHat_lin) {
    VectorXd result(4);
    double pi = atan(1) * 4;
    double a, b, c, d;
    a = xHat_lin(0);
    b = xHat_lin(1);
    c = xHat_lin(3);
    d = xHat_lin(4);
    //calculate sx,sy
    result(0) = sqrt(a * a + b * b);
    result(1) = sqrt(b * b + d * d);
    result(2) = atan2(a * b + c * d, a * d - b * c)*180/pi;
    result(3) = atan2(c, a) * 180 / pi;

    return result;
}

MatrixXd Affine_Transformation(MatrixXd abcd, VectorXd xf, VectorXd yf, VectorXd translation) {
    MatrixXd result(xf.rows(), 2);

    VectorXd current_xy(2);
    
    for (int i = 0; i < xf.rows(); i++) {
        current_xy(0) = xf(i);
        current_xy(1) = yf(i);
        result.row(i) = abcd * current_xy + translation;
    }
    return result;
}

MatrixXd Corrections_Calculator(VectorXd x, VectorXd y, double xp, double yp, double k0, double k1, double k2, double k3, double p1, double p2, double c, double H, double h) {
    MatrixXd result(x.rows(), 10);
    double xbar, ybar, r, K, val;
    //converting H and h to km
    H = H / 1000;
    h = h / 1000;
    K = 2410 * H / (H * H - 6 * H + 250) - 2410 / (h * h - 6 * h + 250) * (h / H);
    //Converting K from microradians to radians:
    K = K / pow(10, 6);
    for (int i = 0; i < result.rows(); i++) {
        xbar = x(i) - xp;
        ybar = y(i) - xp;
        r = sqrt(xbar * xbar + ybar * ybar);
        //Principle point offset correction:
        result(i, 0) = -xp;
        result(i, 1) = -yp;
        //Radial lens distortion: 
        result(i, 2) = -1 * abs(xbar * (k0 + k1 * pow(r, 2) + k2 * pow(r, 4) + k3 * pow(r, 6)));
        result(i, 3) = -1 * abs(ybar * (k0 + k1 * pow(r, 2) + k2 * pow(r, 4) + k3 * pow(r, 6)));
        //Decentring lens distortion: 
        result(i, 4) = -1 * abs((p1 * (r * r + 2 * pow(xbar, 2)) + 2 * p2 * xbar * ybar));
        result(i, 5) = -1 * abs((p2 * (r * r + 2 * pow(ybar, 2)) + 2 * p1 * xbar * ybar));
        //Atmospheric Refracrion:
        result(i, 6) = -1 * abs(xbar * K * (1 + r * r / (c * c)));
        result(i, 7) = -1 * abs(ybar * K * (1 + r * r / (c * c)));
        //Sum of the corrections to get the corrected coordinates:
        result(i, 8) = x(i) + result(i, 0) + result(i, 2) + result(i, 4) + result(i, 6);
        result(i, 9) = y(i) + result(i, 1) + result(i, 3) + result(i, 5) + result(i, 7);

    }
    return result;
}