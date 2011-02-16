package de.varylab.matrix.sparse.factory;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.sparse.IterativeSolver;
import no.uib.cipr.matrix.sparse.IterativeSolverNotConvergedException;
import de.jreality.math.Quaternion;

public class QuaternionMatrixUtils {

	/**
	 * Checks whether a given real Matrix is in fact a quaternion valued one.
	 * 
	 * @param M
	 * @return is quaternionic?
	 */
	public static boolean isQuaternionic(Matrix M) {
		if (M.numRows() == 0 || M.numColumns() == 0)
			return true;
		if (M.numRows() % 4 != 0 || M.numColumns() % 4 != 0)
			return false;
		int numQuatRows = M.numRows() / 4;
		int numQuatColumns = M.numColumns() / 4;
		for (int i = 0; i < numQuatRows; i++) {
			for (int j = 0; j < numQuatColumns; j++) {
				// real part
				if (M.get(4 * i + 0, 4 * j + 0) != M.get(4 * i + 1, 4 * j + 1))
					return false;
				if (M.get(4 * i + 1, 4 * j + 1) != M.get(4 * i + 2, 4 * j + 2))
					return false;
				if (M.get(4 * i + 2, 4 * j + 2) != M.get(4 * i + 3, 4 * j + 3))
					return false;
				// i part
				if (M.get(4 * i + 0, 4 * j + 1) != -M.get(4 * i + 1, 4 * j + 0))
					return false;
				if (M.get(4 * i + 1, 4 * j + 0) != M.get(4 * i + 2, 4 * j + 3))
					return false;
				if (M.get(4 * i + 2, 4 * j + 3) != -M.get(4 * i + 3, 4 * j + 2))
					return false;
				// j part
				if (M.get(4 * i + 0, 4 * j + 2) != -M.get(4 * i + 2, 4 * j + 0))
					return false;
				if (M.get(4 * i + 2, 4 * j + 0) != M.get(4 * i + 3, 4 * j + 1))
					return false;
				if (M.get(4 * i + 3, 4 * j + 1) != -M.get(4 * i + 1, 4 * j + 3))
					return false;
				// k part
				if (M.get(4 * i + 0, 4 * j + 3) != -M.get(4 * i + 3, 4 * j + 0))
					return false;
				if (M.get(4 * i + 3, 4 * j + 0) != M.get(4 * i + 1, 4 * j + 2))
					return false;
				if (M.get(4 * i + 1, 4 * j + 2) != -M.get(4 * i + 2, 4 * j + 1))
					return false;
			}
		}
		return true;
	}

	/**
	 * Solves the system M*P=Q, all entries are quaternion valued.
	 * 
	 * @param M
	 * @param Q
	 * @param P
	 * @param solver
	 */
	public static void solve(Matrix M, Quaternion[] Q, Quaternion[] P,
			IterativeSolver solver) {
		DenseVector x = new DenseVector(4 * P.length);
		DenseVector b = new DenseVector(4 * Q.length);
		for (int i = 0; i < Q.length; i++) {
			b.set(4 * i + 0, Q[i].re);
			b.set(4 * i + 1, -Q[i].x);
			b.set(4 * i + 2, -Q[i].y);
			b.set(4 * i + 3, -Q[i].z);
		}
		try {
			solver.solve(M, b, x);
		} catch (IterativeSolverNotConvergedException e) {
			System.err.println("Iterative solver failed to converge");
		}
		for (int i = 0; i < P.length; i++) {
			P[i] = new Quaternion(x.get(4 * i + 0), -x.get(4 * i + 1),
					-x.get(4 * i + 2), -x.get(4 * i + 3));
		}
	}
}
