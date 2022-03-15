#include "Header.h"
//Claire Mah, Jad Shehadeh, Vincent Cung, Zian Zahid
//ENGO 431 lab 2
//Winter 2021
//Professor: Dr. Lichti
//TA: Kate Pexman

int main() {
	//Average height from given coordinates (in [m])
	double h = 1089.55;
	//Flying height from lab 1 (courtesy of Claire Mah)
	double H = 751.544;
	//Obtaining the coordinates and standard deviations from text files:
		//Check points and standard deviations in comparator image space
	MatrixXd Check_img27 = Matrix_readIn(5, "in_check_img27.txt");
	MatrixXd Check_img28 = Matrix_readIn(5, "in_check_img28.txt");
		//Control points and standard deviations in comparator image space
	MatrixXd Control_img27 = Matrix_readIn(5, "in_control_img27.txt");
	MatrixXd Control_img28 = Matrix_readIn(5, "in_control_img28.txt");
		//Fiducial mark points and standard deviations in comparator image space
	MatrixXd Fid_img27 = Matrix_readIn(5, "in_fiducial_img27.txt");
	MatrixXd Fid_img28 = Matrix_readIn(5, "in_fiducial_img28.txt");
		//Tie points and standard deviations in comparator image space
	MatrixXd Tie_img27 = Matrix_readIn(5, "in_tie_img27.txt");
	MatrixXd Tie_img28 = Matrix_readIn(5, "in_tie_img28.txt");
		//Test values for part d:
	MatrixXd Fid_Meas_test = Matrix_readIn(3, "in_test_comparator.txt");
	VectorXd lobs_test = Matrix_readIn(1, "in_test_reseau.txt");

		//CalCert: Values from calibration certificate, where:
			/*	
				index 0 = focal length, c, in mm
				index 1-3 = Symmetrical Radial Distortion values K0,K1,K2 (respectively)
				index 4-5 = Decentering distortion values P1, P2 (respectively)
			*/
	VectorXd CalCert = Matrix_readIn(1, "in_calibration_certificate_values.txt");
	double c, k0, k1, k2, k3, p1, p2;
	c = CalCert(0);
	k0 = CalCert(1);
	k1 = CalCert(2);
	k2 = CalCert(3);
	k3 = 0;
	p1 = CalCert(4);
	p2 = CalCert(5);
		//CalCert_coords: Coordinate values from calibration certificate, where:
		/*
			index 0 = principle point coordinates in mm
			index 1-8 = Fiducial mark coordinates in mm
		*/
	MatrixXd CalCert_coords = Matrix_readIn(2, "in_calibration_certificate_values2.txt");
	double xp = CalCert_coords(0, 0);
	double yp = CalCert_coords(0, 1);
		//Converting to meters:
	//CalCert_coords = CalCert_coords / 1000;

	cout << "--------- PART C) ----------" << endl;
	//Part c)
		//Getting the linear parameters a,b,c,d,delta x, delta y:
			//Starting by initializing the observation vector (calibration fiducial coordinates)
			//l_obs stores fiducial mark coordinates in the pattern [xf1, yf1, xf2, yf2...]
	
	VectorXd l_obs(CalCert_coords.rows() * 2 - 2);
	for (int i = 0; i < CalCert_coords.rows()-1; i++) {
		l_obs(2 * i) = CalCert_coords(i + 1, 0);
		l_obs(2 * i + 1) = CalCert_coords(i + 1, 1);
	}
	int n = l_obs.rows();
	MatrixXd A_img27 = A_mat_calc(Fid_img27);
	MatrixXd A_img28 = A_mat_calc(Fid_img28);

	MatrixXd P(n, n);
			//Making P an identity matrix:
	P.fill(0);
	for (int i = 0; i < P.rows(); i++) {
		P(i, i) = 1;
	}
			//Getting xhat and residuals
	VectorXd xHat_lin_img27 = xHat_Calc(A_img27, P, l_obs);
	cout << "linear parameters [a, b, delta x, c, d, delta y] for image 27: " << endl << xHat_lin_img27 << endl << endl;
	VectorXd xHat_lin_img28 = xHat_Calc(A_img28, P, l_obs);
	cout << "linear parameters [a, b, delta x, c, d, delta y] for image 28: " << endl << xHat_lin_img28 << endl << endl;

	VectorXd v_Hat_img27 = A_img27 * xHat_lin_img27 - l_obs;
	VectorXd v_Hat_img28 = A_img28 * xHat_lin_img28 - l_obs;
	cout << "v Hat for image 27: " << endl << v_Hat_img27 << endl << endl;
	cout << "v Hat for image 28: " << endl << v_Hat_img28 << endl << endl;

		//Getting non-linear parameters sx, sy, delta, and theta
	VectorXd xHat_nonlin_img27=xHat_nonlin_Calc(xHat_lin_img27);
	VectorXd xHat_nonlin_img28=xHat_nonlin_Calc(xHat_lin_img28);

	cout << "Non-linear parameters (sx,sy,delta,theta) from img 27: " << endl << xHat_nonlin_img27 << endl << endl;
	cout << "Non-linear parameters (sx,sy,delta,theta) from img 28: " << endl << xHat_nonlin_img28 << endl << endl;
		//Combining the parameter matrices to get the 6 unknowns we need:
	VectorXd xHat_img27(6);
	VectorXd xHat_img28(6);
	for (int i = 0; i < 6; i++) {
		if (i < 4) {
			xHat_img27(i) = xHat_nonlin_img27(i);
			xHat_img28(i) = xHat_nonlin_img28(i);
		}
		else {
			if (i == 4) {
				xHat_img27(i) = xHat_lin_img27(2);
				xHat_img28(i) = xHat_lin_img28(2);
			}
			if (i == 5) {
				xHat_img27(i) = xHat_lin_img27(5);
				xHat_img28(i) = xHat_lin_img28(5);
			}
		}
	}
	cout << "Final img 27 parameters [sx, sy, delta, theta, delta x, delta y]" << endl << xHat_img27 << endl << endl;
	cout << "Final img 28 parameters [sx, sy, delta, theta, delta x, delta y]" << endl << xHat_img28 << endl << endl;

	MatToText(xHat_img27, "out_xHat_img27.txt");
	MatToText(xHat_img28, "out_xHat_img28.txt");
	MatToText(v_Hat_img27, "out_vHat_img27.txt");
	MatToText(v_Hat_img28, "out_vHat_img28.txt");

		//Applying transformation to all coordinates
	MatrixXd abcd_img27(2, 2);
	VectorXd translation_img27(2);
	MatrixXd abcd_img28(2, 2);
	VectorXd translation_img28(2);
	abcd_img27 << xHat_lin_img27(0), xHat_lin_img27(1),
		xHat_lin_img27(3), xHat_lin_img27(4);

	abcd_img28 << xHat_lin_img28(0), xHat_lin_img28(1),
		xHat_lin_img28(3), xHat_lin_img28(4);

	translation_img27 << xHat_img27(4), xHat_img27(5);
	translation_img28 << xHat_img28(4), xHat_img28(5);

	cout << "----Image 27 Points converted to fiducial from comparator (all in mm):----" << endl;
	MatrixXd Check_f_img27 = Affine_Transformation(abcd_img27, Check_img27.col(1), Check_img27.col(2), translation_img27);
	MatrixXd Control_f_img27 = Affine_Transformation(abcd_img27, Control_img27.col(1), Control_img27.col(2), translation_img27);
	MatrixXd Fid_f_img27 = Affine_Transformation(abcd_img27, Fid_img27.col(1), Fid_img27.col(2), translation_img27);
	MatrixXd Tie_f_img27 = Affine_Transformation(abcd_img27, Tie_img27.col(1), Tie_img27.col(2), translation_img27);
	cout << "Check points: " << endl << Check_f_img27 << endl << endl;
	cout << "Control points: " << endl << Control_f_img27 << endl << endl;
	cout << "Fiducial points: " << endl << Fid_f_img27 << endl << endl;
	cout << "Tie points: " << endl << Tie_f_img27 << endl << endl;

	cout << "----Image 28 Points converted to fiducial from comparator (all in mm):----" << endl;
	MatrixXd Check_f_img28 = Affine_Transformation(abcd_img28, Check_img28.col(1), Check_img28.col(2), translation_img28);
	MatrixXd Control_f_img28 = Affine_Transformation(abcd_img28, Control_img28.col(1), Control_img28.col(2), translation_img28);
	MatrixXd Fid_f_img28 = Affine_Transformation(abcd_img28, Fid_img28.col(1), Fid_img28.col(2), translation_img28);
	MatrixXd Tie_f_img28 = Affine_Transformation(abcd_img28, Tie_img28.col(1), Tie_img28.col(2), translation_img28);
	cout << "Check points: " << endl << Check_f_img28 << endl << endl;
	cout << "Control points: " << endl << Control_f_img28 << endl << endl;
	cout << "Fiducial points: " << endl << Fid_f_img28 << endl << endl;
	cout << "Tie points: " << endl << Tie_f_img28 << endl << endl;
	



	//Part d:
	cout << endl << "---------- PART D) ----------" << endl << endl;
	MatrixXd A_test = A_mat_calc(Fid_Meas_test);
	VectorXd xHat_lin_test = xHat_Calc(A_test, lobs_test);
	VectorXd xHat_nonlin_test = xHat_nonlin_Calc(xHat_lin_test);
	VectorXd v_Hat_test = A_test * xHat_lin_test - lobs_test;
	VectorXd xHat_test(6);
	for (int i = 0; i < 6; i++) {
		if (i < 4)
			xHat_test(i) = xHat_nonlin_test(i);
		else
			
			xHat_test(i) = xHat_lin_test(i);
	}
	xHat_test(4) = xHat_lin_test(2);
	xHat_test(5) = xHat_lin_test(5);

	//MatrixXd v_Hat_test_answer=Matrix_readIn()
	cout << "Test linear xHat (a,b,delta x, c, d, delta y): " << endl << xHat_lin_test << endl << endl;
	cout << "Test xHat (sx,sy,delta,theta,delta x, delta y): " << endl << xHat_test << endl << endl;

	cout << "Test residuals: " << endl << v_Hat_test << endl;

	MatrixXd v_Hat_test_ans = Matrix_readIn(1, "in_test_vHat_answers.txt");
	int num_of_correct_residuals = 0;
	for (int i = 0; i < v_Hat_test_ans.rows(); i++) {
		if (abs(v_Hat_test(i) - v_Hat_test_ans(i)) <= 0.0001)
			num_of_correct_residuals++;
	}
	if (num_of_correct_residuals = v_Hat_test_ans.rows())
		cout << "Test residuals match computed residuals. Algorithm is correct" << endl;
	else
		cout << "Test residuals do not match computed residuals. Algorithm is incorrect" << endl;
	MatToText(v_Hat_test, "out_test_vHat.txt");



	cout << endl << "----------PART E)----------" << endl << endl;
	


		//Computing corrections:
	MatrixXd Check_corr_img27 = Corrections_Calculator(Check_f_img27.col(0), Check_f_img27.col(1), xp, yp, k0, k1, k2, k3, p1, p2, c, H, h);
	MatrixXd Check_corr_img28 = Corrections_Calculator(Check_f_img28.col(0), Check_f_img28.col(1), xp, yp, k0, k1, k2, k3, p1, p2, c, H, h);

	MatrixXd Control_corr_img27 = Corrections_Calculator(Control_f_img27.col(0), Control_f_img27.col(1), xp, yp, k0, k1, k2, k3, p1, p2, c, H, h);
	MatrixXd Control_corr_img28 = Corrections_Calculator(Control_f_img28.col(0), Control_f_img28.col(1), xp, yp, k0, k1, k2, k3, p1, p2, c, H, h);

	MatrixXd Fid_corr_img27 = Corrections_Calculator(Fid_f_img27.col(0), Fid_f_img27.col(1), xp, yp, k0, k1, k2, k3, p1, p2, c, H, h);
	MatrixXd Fid_corr_img28 = Corrections_Calculator(Fid_f_img28.col(0), Fid_f_img28.col(1), xp, yp, k0, k1, k2, k3, p1, p2, c, H, h);

	MatrixXd Tie_corr_img27 = Corrections_Calculator(Tie_f_img27.col(0), Tie_f_img27.col(1), xp, yp, k0, k1, k2, k3, p1, p2, c, H, h);
	MatrixXd Tie_corr_img28 = Corrections_Calculator(Tie_f_img28.col(0), Tie_f_img28.col(1), xp, yp, k0, k1, k2, k3, p1, p2, c, H, h);

		//Printing to Console
	cout << "----- THE FOLLOWING VALUES ARE ALL IN [mm] -----" << endl << endl;
	cout << "Corrections for Check points: " << endl;
	cout << " [principle point offset]  [radial lens distortion]       [decentering lens]          [atmospheric]         [Corrected Coordinates]" << endl;
	cout << "Image 27:" << endl;
	cout << Check_corr_img27 << endl;
	cout << "Image 28:" << endl;
	cout << Check_corr_img28 << endl << endl;

	cout << "Corrections for Control points: " << endl;
	cout << " [principle point offset]  [radial lens distortion]       [decentering lens]          [atmospheric]         [Corrected Coordinates]" << endl;
	cout << "Image 27:" << endl;
	cout << Control_corr_img27 << endl;
	cout << "Image 28:" << endl;
	cout << Control_corr_img28 << endl << endl;

	cout << "Corrections for Fiducial points: " << endl;
	cout << " [principle point offset]  [radial lens distortion]       [decentering lens]          [atmospheric]         [Corrected Coordinates]" << endl;
	cout << "Image 27:" << endl;
	cout << Fid_corr_img27 << endl;
	cout << "Image 28:" << endl;
	cout << Fid_corr_img28 << endl << endl;

	cout << "Corrections for Tie points: " << endl;
	cout << " [principle point offset]  [radial lens distortion]       [decentering lens]          [atmospheric]         [Corrected Coordinates]" << endl;
	cout << "Image 27:" << endl;
	cout << Tie_corr_img27 << endl;
	cout << "Image 28:" << endl;
	cout << Tie_corr_img28 << endl << endl;

	MatToText(Check_corr_img27, "out_Check_corrections_img27.txt");
	MatToText(Check_corr_img28, "out_Check_corrections_img28.txt");
	MatToText(Control_corr_img27, "out_Control_corrections_img27.txt");
	MatToText(Control_corr_img28, "out_Control_corrections_img28.txt");
	MatToText(Fid_corr_img27, "out_Ficudial_corrections_img27.txt");
	MatToText(Fid_corr_img28, "out_Fiducial_corrections_img28.txt");
	MatToText(Tie_corr_img27, "out_Tie_corrections_img27.txt");
	MatToText(Tie_corr_img28, "out_Tie_corrections_img28.txt");

	int ex;
	cin >> ex;
	return 0;
}