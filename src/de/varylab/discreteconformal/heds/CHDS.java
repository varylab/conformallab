package de.varylab.discreteconformal.heds;

import static java.lang.Math.PI;
import static java.lang.Math.atan2;
import static java.lang.Math.exp;
import static java.lang.Math.log;
import static java.lang.Math.sqrt;
import geom3d.Point;
import geom3d.Triangle;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.varylab.discreteconformal.math.Lob;

public class CHDS extends HalfEdgeDataStructure<CVertex, CEdge, CFace> {

	public CHDS() {
		super(CVertex.class, CEdge.class, CFace.class);
	}

	
	public void prepareData(final Vector theta) {
		// set initial lambdas
		for (final CEdge e : getPositiveEdges()) {
			final double l = e.getLength();
			e.setLambda(2.0 * log(l));
			e.getOppositeEdge().setLambda(e.getLambda());
		}
		// set thetas
		for (final CVertex v : getVertices()) {
			v.setTheta(theta.get(v.getIndex()));
		}
		// set initial face energy
		for (final CFace f : getFaces()) {
			final CEdge 
				e1 = f.getBoundaryEdge(),
				e2 = e1.getNextEdge(),
				e3 = e1.getPreviousEdge();
			final Point 
				p1 = e1.getTargetVertex().getPosition(),
				p2 = e2.getTargetVertex().getPosition(),
				p3 = e3.getTargetVertex().getPosition();
			final double 
				a1 = Triangle.angleAt(p2, p1, p3),
				a2 = Triangle.angleAt(p3, p2, p1),
				a3 = Triangle.angleAt(p1, p3, p2);
			final double 
				E1 = a1*e1.getLambda() + a2*e2.getLambda() + a3*e3.getLambda(),
				E2 = Lob.valueAt(a1) + Lob.valueAt(a2) + Lob.valueAt(a3);
			f.setEnergy(E1 + E2);
		}
	}
	
	
	private void cotEntries(final Vector u, final CFace f, final double[] cotE, final double[] cotV) {
		final CEdge
			e1 = f.getBoundaryEdge(),
			e2 = e1.getNextEdge(),
			e3 = e1.getPreviousEdge();
		final CVertex 
			v1 = e1.getTargetVertex(),
			v2 = e2.getTargetVertex(),
			v3 = e3.getTargetVertex();
		final double 
			u1 = u.get(v1.getIndex()),
			u2 = u.get(v2.getIndex()),
			u3 = u.get(v3.getIndex());
		final double 
			x12 = e2.getLambda() + u1 + u2,
			x23 = e3.getLambda() + u2 + u3,
			x31 = e1.getLambda() + u3 + u1;
		final double 
			xmean = (x12 + x23 + x31) / 3;
		final double 
			l12 = exp(x12 - xmean),
			l23 = exp(x23 - xmean),
			l31 = exp(x31 - xmean);
		final double
			t31 = +l12+l23-l31,
			t23 = +l12-l23+l31,
			t12 = -l12+l23+l31;
		double 
			cot1 = 0.0,
			cot2 = 0.0,
			cot3 = 0.0;
		if (t31 > 0 && t23 > 0 && t12 > 0) {
			final double
				l123 = l12 + l23 + l31,
				denom = sqrt(t12 * t23 * t31 * l123) * 2;
			cot1 = (t23*l123 - t31*t12) / denom;
			cot2 = (t31*l123 - t12*t23) / denom;
			cot3 = (t12*l123 - t23*t31) / denom;
		}
		cotE[0] = cot2;
		cotE[1] = cot3;
		cotE[2] = cot1;
		cotV[0] = cot2 + cot3;
		cotV[1] = cot3 + cot1;
		cotV[2] = cot1 + cot2;
	}
	
	
	
	private double triangleEnergyAndAlphas(final Vector u, final CFace f, final double[] a123) {
		final CEdge 
			e1 = f.getBoundaryEdge(),
			e2 = e1.getNextEdge(),
			e3 = e1.getPreviousEdge();
		final CVertex 
			v1 = e1.getTargetVertex(),
			v2 = e2.getTargetVertex(),
			v3 = e3.getTargetVertex();
		a123[0] = 0.0;
		a123[1] = 0.0;
		a123[2] = 0.0;
		final double 
			u1 = u.get(v1.getIndex()),
			u2 = u.get(v2.getIndex()),
			u3 = u.get(v3.getIndex());
		final double 
			umean = (u1+u2+u3)/3;
		final double 
			x12 = e2.getLambda() + (u1+u2 - 2*umean),
			x23 = e3.getLambda() + (u2+u3 - 2*umean),
			x31 = e1.getLambda() + (u3+u1 - 2*umean);
		final double 
			l12 = exp(x12),
			l23 = exp(x23),
			l31 = exp(x31);
		final double 
			t31 = +l12+l23-l31,
			t23 = +l12-l23+l31,
			t12 = -l12+l23+l31;
		if (t31 > 0 && t23 > 0 && t12 > 0) {
			final double 
				l123 = l12 + l23 + l31,
				denom = sqrt(t12 * t23 * t31 * l123);
			a123[0] = 2 * atan2(t12 * t31, denom);
			a123[1] = 2 * atan2(t23 * t12, denom);
			a123[2] = 2 * atan2(t31 * t23, denom);
		} else if (t31 <= 0) {
			a123[1] = 0.0;
		} else if (t23 <= 0) {
			a123[0] = 0.0;
		} else if (t12 <= 0) {
			a123[2] = 0.0;
		}
		final double 
			E1 = a123[0]*x23 + a123[1]*x31 + a123[2]*x12,
			E2 = Lob.valueAt(a123[0]) + Lob.valueAt(a123[1]) + Lob.valueAt(a123[2]),
			E3 = - PI * umean - f.getEnergy();
		return E1 + E2 + E3; 
	}
	
	
	
	
	public void conformalEnergy(final Vector u, final double[] E, final Vector G, final Matrix H) {
		// Vertex Energy
		if (E != null) 
			E[0] = 0.0;
		if (G != null)
			G.zero();
		if (H != null)
			H.zero();
		for (final CVertex v : getVertices()) {
			if (E != null)
				E[0] += v.getTheta() * u.get(v.getIndex());
			if (G != null)
				G.add(v.getIndex(), v.getTheta());
		}
		if (E != null)
			E[0] *= 0.5;
		// Face Energy
		final double[] a123 = {0, 0, 0};
		for (final CFace t : getFaces()) {
			final CEdge 
				e1 = t.getBoundaryEdge(),
				e2 = e1.getNextEdge(),
				e3 = e1.getPreviousEdge();
			final CVertex 
				v1 = e1.getTargetVertex(),
				v2 = e2.getTargetVertex(),
				v3 = e3.getTargetVertex();
			final double e = triangleEnergyAndAlphas(u, t, a123);
			if (E != null)
				E[0] += e;
			if (G != null) {
				G.add(v1.getIndex(), -a123[0]);
				G.add(v2.getIndex(), -a123[1]);
				G.add(v3.getIndex(), -a123[2]);
			}
			if (H != null) {
				final double[] 
				     cotE = {0, 0, 0},
				     cotV = {0, 0, 0};
				cotEntries(u, t, cotE, cotV);
				// edge hessian
				H.add(v1.getIndex(), v3.getIndex(), cotE[0]);
				H.add(v3.getIndex(), v1.getIndex(), cotE[0]);
				H.add(v2.getIndex(), v1.getIndex(), cotE[1]);
				H.add(v1.getIndex(), v2.getIndex(), cotE[1]);
				H.add(v2.getIndex(), v3.getIndex(), cotE[2]);
				H.add(v3.getIndex(), v2.getIndex(), cotE[2]);
				// vertex hessian
				H.add(v1.getIndex(), v1.getIndex(), cotV[0]);
				H.add(v2.getIndex(), v2.getIndex(), cotV[1]);
				H.add(v3.getIndex(), v3.getIndex(), cotV[2]);
			}
		}
	}
	
	
	
}
