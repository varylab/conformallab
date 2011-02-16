package de.varylab.matrix.sparse.factory;

import no.uib.cipr.matrix.sparse.CompRowMatrix;
import de.jreality.math.Quaternion;

public class MatrixRepresentationFactory {

	public static double[][] makeDoubleArrayRepresentation(Quaternion Q) {
		double[][] block = new double[4][4];

		block[0][0] = Q.re;
		block[1][1] = Q.re;
		block[2][2] = Q.re;
		block[3][3] = Q.re;

		block[0][1] = Q.x;
		block[1][0] = -Q.x;
		block[2][3] = -Q.x;
		block[3][2] = Q.x;

		block[0][2] = Q.y;
		block[1][3] = Q.y;
		block[2][0] = -Q.y;
		block[3][1] = -Q.y;

		block[0][3] = Q.z;
		block[1][2] = -Q.z;
		block[2][1] = Q.z;
		block[3][0] = -Q.z;

		return block;
	}

	public static CompRowMatrix makeCompRowMatrixRepresentation(Quaternion Q) {
		int size = 0;
		if (Q.re != 0)
			size++;
		if (Q.x != 0)
			size++;
		if (Q.y != 0)
			size++;
		if (Q.z != 0)
			size++;
		int[][] blockid = new int[4][size];
		size = 0;
		if (Q.re != 0) {
			blockid[0][size] = 0;
			blockid[1][size] = 1;
			blockid[2][size] = 2;
			blockid[3][size] = 3;
			size++;
		}
		if (Q.x != 0) {
			blockid[0][size] = 1;
			blockid[1][size] = 0;
			blockid[2][size] = 3;
			blockid[3][size] = 2;
			size++;
		}
		if (Q.y != 0) {
			blockid[0][size] = 2;
			blockid[1][size] = 3;
			blockid[2][size] = 0;
			blockid[3][size] = 1;
			size++;
		}
		if (Q.z != 0) {
			blockid[0][size] = 3;
			blockid[1][size] = 2;
			blockid[2][size] = 1;
			blockid[3][size] = 0;
		}
		CompRowMatrix block = new CompRowMatrix(4, 4, blockid);
		if (Q.re != 0) {
			block.set(0, 0, Q.re);
			block.set(1, 1, Q.re);
			block.set(2, 2, Q.re);
			block.set(3, 3, Q.re);
		}
		if (Q.x != 0) {
			block.set(0, 1, Q.x);
			block.set(1, 0, -Q.x);
			block.set(2, 3, -Q.x);
			block.set(3, 2, Q.x);
		}
		if (Q.y != 0) {
			block.set(0, 2, Q.y);
			block.set(1, 3, Q.y);
			block.set(2, 0, -Q.y);
			block.set(3, 1, -Q.y);
		}
		if (Q.z != 0) {
			block.set(0, 3, Q.z);
			block.set(1, 2, -Q.z);
			block.set(2, 1, Q.z);
			block.set(3, 0, -Q.z);
		}
		return block;
	}

}
