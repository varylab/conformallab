package de.varylab.discreteconformal.hyperelliptic;

import de.jtem.mfc.field.Complex;
import de.jtem.riemann.constrainedComplex.Constrain;
import de.jtem.riemann.constrainedComplex.ConstrainFactory;
import de.jtem.riemann.constrainedComplex.ConstrainedComplex;
import de.jtem.riemann.surface.BranchPoint;
import de.jtem.riemann.surface.homologie.Transform;

public class CurveUtils {

	public static boolean equal(Curve c1,
			Curve c2) {

		if (c1 == null || c2 == null)
			return false;

		c1.update();
		c2.update();

		for (int i = 0; i < c1.getNumOfBranchPoints(); i++) {
			if (!c1.getBranchPoint(i).equals(c2.getBranchPoint(i)))
				return false;
		}

		if (c1.getGenus() != c2.getGenus())
			return false;

		return true;
	}

	public static BranchPoint[] getUnfixedBranchPoints(Curve curve) {
		int counter = 0;
		for (int i = 0; i < curve.getNumOfBranchPoints(); i++) {
			if (!ConstrainFactory.isFixConstrained(curve.getBranchPoint(i)))
				counter++;
		}
		BranchPoint[] res = new BranchPoint[counter];
		for (int i = 0, j = 0; i < curve.getNumOfBranchPoints(); i++) {
			if (!ConstrainFactory.isFixConstrained(curve.getBranchPoint(i))) {
				res[j] = curve.getBranchPoint(i);
				j++;
			}
		}
		return res;
	}

	public static Complex[] getUnfixedBranchPointCoordsInUpperHalfPlane(
			Curve curve) {
		BranchPoint[] oldBranchPoints = curve.getBranchPoints();
		int counter = 0;
		for (int i = 0; i < oldBranchPoints.length; i++) {
			if (oldBranchPoints[i].getCoords().im >= 0
					&& !ConstrainFactory.isFixConstrained(oldBranchPoints[i]))
				counter++;
		}
		Complex[] c = new Complex[counter];
		int j = 0;
		for (int i = 0; i < oldBranchPoints.length; i++) {
			if (oldBranchPoints[i].getCoords().im >= 0
					&& !ConstrainFactory.isFixConstrained(oldBranchPoints[i])) {
				c[j] = oldBranchPoints[i].getCoords();
				j++;
			}
		}
		return c;
	}

	public static Curve copyCurve(Curve curve) {

		// initialize a new curve

		Curve copy = new Curve(curve.getGenus());
		
		try {
			copyCurve(curve, copy);
		} catch (Exception e) {
			System.err.println(e.getMessage());
			e.printStackTrace();
		}

		return copy;

	}

	public static void copyCurve(Curve curve,
			Curve copy) {

		copy.update();

		// handle the constrained branch points

		boolean[] wasAlreadyRecognized = new boolean[curve
				.getNumOfBranchPoints()];
		for (int i = 0; i < curve.getNumOfBranchPoints(); i++) {
			Constrain c = curve.getBranchPoint(i).getConstrain();
			ConstrainedComplex[] z = ConstrainFactory
					.getConstrainedComplexOf(c);
			int[] id = new int[z.length];
			boolean ok = true;
			for (int j = 0; j < z.length; j++) {
				id[j] = getIdOfBranchPoint(curve, z[j]);
				if (id[j] < 0)
					ok = false;
				if (wasAlreadyRecognized[id[j]])
					ok = false;
			}
			if (ok) {
				for (int j = 0; j < id.length; j++) {
					int ID = id[j];
					copy.getBranchPoint(ID).deconstrain();
					copy.getBranchPoint(ID).setCoords(
							curve.getBranchPoint(ID).getCoords());
					wasAlreadyRecognized[ID] = true;
				}
				if (ConstrainFactory.getConstrainTypeOf(c) == ConstrainFactory.IMAGINARY) {
					ConstrainFactory.setImaginary(copy
							.getBranchPoint(id[0]));
				}
				if (ConstrainFactory.getConstrainTypeOf(c) == ConstrainFactory.FIX) {
					ConstrainFactory.setFix(copy.getBranchPoint(id[0]));
				}
				if (ConstrainFactory.getConstrainTypeOf(c) == ConstrainFactory.REAL) {
					ConstrainFactory.setReal(copy.getBranchPoint(id[0]));
				}
				if (ConstrainFactory.getConstrainTypeOf(c) == ConstrainFactory.IMAGINARYSYMMETRY) {
					ConstrainFactory.setImaginarySymmetry(copy
							.getBranchPoint(id[0]), copy
							.getBranchPoint(id[1]));
				}
				if (ConstrainFactory.getConstrainTypeOf(c) == ConstrainFactory.POINTSYMMETRY) {
					ConstrainFactory.setPointSymmetry(copy
							.getBranchPoint(id[0]), copy
							.getBranchPoint(id[1]));
				}
				if (ConstrainFactory.getConstrainTypeOf(c) == ConstrainFactory.REALSYMMETRY) {
					ConstrainFactory.setRealSymmetry(copy
							.getBranchPoint(id[0]), copy
							.getBranchPoint(id[1]));
				}
				if (ConstrainFactory.getConstrainTypeOf(c) == ConstrainFactory.UNITCIRCLESYMMETRY) {
					ConstrainFactory.setUnitCircleSymmetry(copy
							.getBranchPoint(id[0]), copy
							.getBranchPoint(id[1]));
				}
			}
		}

		// copy the transform

		Transform t = curve.getTransform() == null ? null : (Transform) curve
				.getTransform().clone();

		copy.getTransform().assign(t);

		copy.outdate();
		copy.update();

	}

	public static int getIdOfBranchPoint(Curve curve,
			ConstrainedComplex z) {
		for (int i = 0; i < curve.getNumOfBranchPoints(); i++) {
			if (z == curve.getBranchPoint(i))
				return i;
		}
		return -1;
	}
}
