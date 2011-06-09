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
		DomainValue x, 
		Energy E, 
		Gradient G, 
		Hessian H
	) {
		double[] p = {x.get(0), x.get(1), x.get(2), x.get(3)};
		System.out.println(Pn.norm(p, Pn.EUCLIDEAN));
		double pp = innerProduct(p, p, HYPERBOLIC);
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
			double[] vp = aSet.getD(Position4d.class, v);
			double vpp = innerProduct(vp, p, HYPERBOLIC);
			if (E != null) {
				E.add(-log(vpp / pp));
			}
			if (G != null) {
				Rn.times(g1, 2/pp, p);
				Rn.times(g2, -1/vpp, vp);
				g1[3] *= -1;
				g2[3] *= -1;
				addVectorToGradient(G, 0, g1);
				addVectorToGradient(G, 0, g2);
			}
			if (H != null) {
				for (int i = 0; i < 4; i++) {
					Rn.times(g1, p[i], p);
					Rn.times(g2, vp[i], vp);
					Rn.times(g1, pp*pp, g1);
					Rn.times(g2, -vpp*vpp, g2);
					g1[3] *= -1;
					g2[3] *= -1;
					addRowToHessian(H, i, g1);
					addRowToHessian(H, i, g2);
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
