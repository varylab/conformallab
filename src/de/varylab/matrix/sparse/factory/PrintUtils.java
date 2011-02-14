package de.varylab.matrix.sparse.factory;

import java.text.NumberFormat;
import java.util.Locale;

import no.uib.cipr.matrix.Matrix;

public class PrintUtils {

	private static NumberFormat numberFormat = NumberFormat
			.getInstance(Locale.ENGLISH);

	public static void print(double[][] array) {
		print(array, 2);
	}

	public static void print(double[][] array, int noOfFractionDigits) {
		numberFormat.setMinimumFractionDigits(noOfFractionDigits);
		numberFormat.setMaximumFractionDigits(noOfFractionDigits);
		String tmp;
		for (int i = 0; i < array.length; i++) {
			System.out.print("| ");
			for (int j = 0; j < array[i].length; j++) {
				if (array[i][j] == 0)
					tmp = " " + numberFormat.format(0.);
				else {
					tmp = array[i][j] >= 0 ? "+" : "";
					tmp += numberFormat.format(array[i][j]);
				}
				System.out.print(tmp);
				if (j < array[i].length - 1) {
					System.out.print(" , ");
				} else {
					System.out.println(" |");
				}
			}
		}
	}

	public static void print(int[][] array) {
		print(array, 2);
	}

	public static void print(int[][] array, int noOfFractionDigits) {
		numberFormat.setMinimumFractionDigits(noOfFractionDigits);
		numberFormat.setMaximumFractionDigits(noOfFractionDigits);
		String tmp;
		for (int i = 0; i < array.length; i++) {
			System.out.print("| ");
			for (int j = 0; j < array[i].length; j++) {
				if (array[i][j] == 0)
					tmp = " " + numberFormat.format(0.);
				else {
					tmp = array[i][j] >= 0 ? "+" : "";
					tmp += numberFormat.format(array[i][j]);
				}
				System.out.print(tmp);
				if (j < array[i].length - 1) {
					System.out.print(" , ");
				} else {
					System.out.println(" |");
				}
			}
		}
	}

	public static void print(Matrix M) {
		print(M, 2);
	}

	public static void print(Matrix M, int noOfFractionDigits) {
		numberFormat.setMinimumFractionDigits(noOfFractionDigits);
		numberFormat.setMaximumFractionDigits(noOfFractionDigits);
		String tmp;
		for (int i = 0; i < M.numRows(); i++) {
			System.out.print("| ");
			for (int j = 0; j < M.numColumns(); j++) {
				if (M.get(i, j) == 0)
					tmp = " " + numberFormat.format(0.);
				else {
					tmp = M.get(i, j) >= 0 ? "+" : "";
					tmp += numberFormat.format(M.get(i, j));
				}
				System.out.print(tmp);
				if (j < M.numColumns() - 1) {
					System.out.print(" , ");
				} else {
					System.out.println(" |");
				}
			}
		}
		System.out.println();
	}
}
