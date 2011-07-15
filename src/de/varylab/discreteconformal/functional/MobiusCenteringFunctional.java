package de.varylab.discreteconformal.functional;

import static de.jreality.math.Pn.HYPERBOLIC;
import static de.jreality.math.Pn.innerProduct;
import static de.jtem.halfedgetools.functional.FunctionalUtils.addRowToHessian;
import static de.jtem.halfedgetools.functional.FunctionalUtils.addVectorToGradient;
import static java.lang.Math.log;
import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.generic.Position4d;
import de.jtem.halfedgetools.functional.DomainValue;
import de.jtem.halfedgetools.functional.Energy;
import de.jtem.halfedgetools.functional.Functional;
import de.jtem.halfedgetools.functional.Gradient;
import de.jtem.halfedgetools.functional.Hessian;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;

public class MobiusCenteringFunctional implements Functional<CoVertex, CoEdge, CoFace> {

	private AdapterSet
		aSet = null;
	private double[]
	    g1 = new double[4],
		g2 = new double[4];
	
	public MobiusCenteringFunctional(AdapterSet a) {
		this.aSet = a;
	}
	
	@Override
	public <HDS extends HalfEdgeDataStructure<CoVertex, CoEdge, CoFace>> void evaluate(
		HDS hds, 
		DomainValue d, 
		Energy E, 
		Gradient G, 
		Hessian H
	) {
		double[] x = {d.get(0), d.get(1), d.get(2), d.get(3)};
		System.out.println(Pn.norm(x, Pn.EUCLIDEAN));
		double xx = innerProduct(x, x, HYPERBOLIC);
		if (E != null) {
			E.setZero();
		}
		if (G != null) {
			G.setZero();
		}
		if (H != null) {
			H.setZero();
		}
		for (CoVertex v : hds.getVertices()) {
			double[] p = aSet.getD(Position4d.class, v);
			double xp = innerProduct(p, x, HYPERBOLIC);
			if (E != null) {
				E.add(-log(xp / xx));
			}
			if (G != null) {
				Rn.times(g1, 2/xx, x);
				Rn.times(g2, -1/xp, p);
				g1[3] *= -1;
				g2[3] *= -1;
				addVectorToGradient(G, 0, g1);
				addVectorToGradient(G, 0, g2);
			}
			if (H != null) {
				for (int i = 0; i < 4; i++) {
					Rn.times(g1, p[i] / (xp*xp), p);
					addRowToHessian(H, i, g1);
					double sign = i == 3 ? -1 : 1;
					H.add(i, i, sign * 2 / xx - 4*x[i]*x[i]/(xx*xx));
				}
			}
		}
	}

	@Override
	public boolean hasHessian() {
		return true;
	}

	@Override
	public <HDS extends HalfEdgeDataStructure<CoVertex, CoEdge, CoFace>> int getDimension(HDS hds) {
		return 4;
	}

	@Override
	public <HDS extends HalfEdgeDataStructure<CoVertex, CoEdge, CoFace>> int[][] getNonZeroPattern(HDS hds) {
		int[][] nnz = new int[4][4];
		for (int i = 0; i < nnz.length; i++) {
			for (int j = 0; j < nnz[i].length; j++) {
				nnz[i][j] = j;
			}
		}
		return nnz;
	}

}
