package de.varylab.discreteconformal.util;

import java.text.NumberFormat;
import java.util.Locale;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import de.jtem.blas.ComplexMatrix;


/**
 * Class to print out colt matrices in a readable form.
 * 
 * @author knoeppel
 * 
 */
public class SimpleMatrixPrintUtility {

	
	/**
	 * Print method for colt matrices
	 * @param M
	 */
	public static void print(DoubleMatrix2D M) {
		
		String tmp;
		for (int i = 0; i < M.rows(); i++) {
			System.out.print("| ");
			for (int j = 0; j < M.columns(); j++) {
				if (M.get(i, j) == 0)
					tmp = " " + 0;
				else {
					tmp = M.get(i, j) >= 0 ? "+" : "";
					tmp += M.get(i, j);
				}
				System.out.print(tmp);
				if (j < M.columns() - 1) {
					System.out.print(" , ");
				} else {
					System.out.println(" |");
				}
			}
		}
		System.out.println();
	}
	
	/**
	 * Print method for colt matrices
	 * @param M
	 * @param noOfFractionDigits
	 */
	public static void print(DoubleMatrix2D M, int noOfFractionDigits) {
		
		NumberFormat numberFormat = NumberFormat
		.getInstance(Locale.ENGLISH);
		
		numberFormat.setMinimumFractionDigits(noOfFractionDigits);
		numberFormat.setMaximumFractionDigits(noOfFractionDigits);
		String tmp;
		for (int i = 0; i < M.rows(); i++) {
			System.out.print("| ");
			for (int j = 0; j < M.columns(); j++) {
				if (M.get(i, j) == 0)
					tmp = " " + numberFormat.format(0.);
				else {
					tmp = M.get(i, j) >= 0 ? "+" : "";
					tmp += numberFormat.format(M.get(i, j));
				}
				System.out.print(tmp);
				if (j < M.columns() - 1) {
					System.out.print(" , ");
				} else {
					System.out.println(" |");
				}
			}
		}
		System.out.println();
	}
	
	/**
	 * Print method for complex matrices
	 * @param M
	 * @param noOfFractionDigits
	 */
	public static void print(ComplexMatrix M, int noOfFractionDigits) {
		System.out.println(toString(M, noOfFractionDigits));
		System.out.println();
	}
	
	/**
	 * Print method for complex matrices
	 * @param M
	 * @param noOfFractionDigits
	 */
	public static String toString(ComplexMatrix M, int noOfFractionDigits) {
		
		NumberFormat numberFormat = NumberFormat.getInstance(Locale.ENGLISH);
		
		numberFormat.setMinimumFractionDigits(noOfFractionDigits);
		numberFormat.setMaximumFractionDigits(noOfFractionDigits);
		
		String newline = System.getProperty("line.separator");
		
		String tmp="";
		for (int i = 0; i < M.getNumRows(); i++) {
			tmp+="| ";
			for (int j = 0; j < M.getNumCols(); j++) {
				tmp += M.get(i, j).re >= 0 ? "+" : "";
				tmp += numberFormat.format(M.get(i, j).re);
				tmp += M.get(i, j).im >= 0 ? "+" : "";
				tmp += numberFormat.format(M.get(i, j).im)+"i";
				if (j < M.getNumCols() - 1) {
					tmp+=" , ";
				} else {
					tmp+=" |"+newline;
				}
			}
		}
		return tmp;
	}
	
	/**
	 * Print method for colt vector
	 * @param V
	 * @param noOfFractionDigits
	 */
	public static void print(DoubleMatrix1D V, int noOfFractionDigits) {

		NumberFormat numberFormat = NumberFormat.getInstance(Locale.ENGLISH);

		numberFormat.setMinimumFractionDigits(noOfFractionDigits);
		numberFormat.setMaximumFractionDigits(noOfFractionDigits);
		String tmp;

		System.out.print("| ");
		for (int i = 0; i < V.size(); i++) {
			if (V.get(i) == 0)
				tmp = " " + numberFormat.format(0.);
			else {
				tmp = V.get(i) >= 0 ? "+" : "";
				tmp += numberFormat.format(V.get(i));
			}
			System.out.print(tmp);
			if (i < V.size() - 1) {
				System.out.print(" , ");
			} else {
				System.out.println(" |");
			}
		}

		System.out.println();
	}

}
