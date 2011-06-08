package de.varylab.discreteconformal.functional;

import static de.jtem.halfedgetools.functional.FunctionalUtils.addRowToHessian;
import static de.jtem.halfedgetools.functional.FunctionalUtils.addVectorToGradient;
import static de.jtem.halfedgetools.functional.FunctionalUtils.subtractRowFromHessian;
import static de.jtem.halfedgetools.functional.FunctionalUtils.subtractVectorFromGradient;
import static java.lang.Math.log;
import static java.lang.Math.sqrt;
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
		double pp = dot4(p, p);
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
			double vpp = dot4(vp, p);
			if (E != null) {
				E.add(log(vpp / sqrt(pp)));
			}
			if (G != null) {
				Rn.copy(g1, p);
				Rn.copy(g2, vp);
				g1[3] *= pp;
				g2[3] *= vpp;
				addVectorToGradient(G, 0, g1);
				subtractVectorFromGradient(G, 0, g2);
			}
			if (H != null) {
				for (int i = 0; i < 4; i++) {
					Rn.times(g1, p[i], p);
					Rn.times(g2, vp[i], vp);
					g1[3] *= pp*pp/2;
					g2[3] *= vpp*vpp;
					addRowToHessian(H, i, g1);
					subtractRowFromHessian(H, i, g2);
				}
				H.add(3, 3, pp);
			}
		}
	}

	public double dot4(double[] u, double[] v) {
		return u[3]*v[3] - u[2]*v[2] - u[1]*v[1] - u[2]*v[2];
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
