package de.varylab.discreteconformal.heds;

import static de.varylab.discreteconformal.math.Lob.lob;
import static java.lang.Math.PI;
import static java.lang.Math.atan2;
import static java.lang.Math.exp;
import static java.lang.Math.log;
import static java.lang.Math.sqrt;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.sparse.SparseVector;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.varylab.discreteconformal.math.optimization.Optimizable;

public class CHDS extends HalfEdgeDataStructure<CVertex, CEdge, CFace> implements Optimizable{

	
	private Integer 
		dim = 0;
	private boolean
		texCoordinatesValid = false;
	
	public CHDS() {
		super(CVertex.class, CEdge.class, CFace.class);
	}

	
	public boolean isVariable(CVertex v) {
		return v.getSolverIndex() >= 0;
	}
	
	
	/**
	 * Compute algorithm invariant data. Boundary is the natural mesh boundary.
	 */
	public void prepareInvariantData() {
		HashSet<CVertex> b = new HashSet<CVertex>();
		b.addAll(HalfEdgeUtils.boundaryVertices(this));
		prepareInvariantData(b);
	}
	
	/**
	 * Compute algorithm invariant data
	 * @param boundary the boundary vertices which do not belong to the solver system
	 */
	public void prepareInvariantData(Collection<CVertex> boundary) {
		// set initial lambdas
		for (final CEdge e : getPositiveEdges()) {
			final double l = e.getLength();
			e.setLambda(log(l));
			e.getOppositeEdge().setLambda(e.getLambda());
		}
		// set thetas and solver indices
		dim = 0;
		for (final CVertex v : getVertices()) {
			if (boundary.contains(v)) {
				v.setTheta(0.0);
				v.setSolverIndex(-1);
			} else {
				v.setTheta(2 * PI);
				v.setSolverIndex(dim++);
			}
		}
		Vector u = new SparseVector(numVertices());
		Map<CEdge, Double> aMap = new HashMap<CEdge, Double>();
		for (final CFace f : getFaces()) {
			f.setEnergy(triangleEnergyAndAlphas(u, f, aMap));
		}
	}
	
	
	protected void triangleHessian(final Vector u, final CFace f, final double[] cotE, final double[] cotV) {
		final CEdge
			e1 = f.getBoundaryEdge(),
			e2 = e1.getNextEdge(),
			e3 = e1.getPreviousEdge();
		final CVertex 
			v1 = e1.getTargetVertex(),
			v2 = e2.getTargetVertex(),
			v3 = e3.getTargetVertex();
		final double 
			u1 = isVariable(v1) ? u.get(v1.getSolverIndex()) : 0.0,
			u2 = isVariable(v2) ? u.get(v2.getSolverIndex()) : 0.0,
			u3 = isVariable(v3) ? u.get(v3.getSolverIndex()) : 0.0;
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
	
	
	
	protected double triangleEnergyAndAlphas(final Vector u, final CFace f, Map<CEdge, Double> alphas) {
		final CEdge 
			e1 = f.getBoundaryEdge(),
			e2 = e1.getNextEdge(),
			e3 = e1.getPreviousEdge();
		final CVertex 
			v1 = e1.getTargetVertex(),
			v2 = e2.getTargetVertex(),
			v3 = e3.getTargetVertex();
		double 
			a1 = 0.0,
			a2 = 0.0,
			a3 = 0.0;
		final double 
			u1 = isVariable(v1) ? u.get(v1.getSolverIndex()) : 0.0,
			u2 = isVariable(v2) ? u.get(v2.getSolverIndex()) : 0.0,
			u3 = isVariable(v3) ? u.get(v3.getSolverIndex()) : 0.0;
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
			a1 = 2 * atan2(t12 * t31, denom);
			a2 = 2 * atan2(t23 * t12, denom);
			a3 = 2 * atan2(t31 * t23, denom);
		} else if (t31 <= 0) {
			a2 = PI;
		} else if (t23 <= 0) {
			a1 = PI;
		} else if (t12 <= 0) {
			a3 = PI;
		}
		final double 
			E1 = a1*x23 + a2*x31 + a3*x12,
			E2 = lob(a1) + lob(a2) + lob(a3),
			E3 = - PI * umean - f.getEnergy();
		alphas.put(e1, a2);
		alphas.put(e2, a3);
		alphas.put(e3, a1);
		return E1 + E2 + E3; 
	}
	
	
	
	
	public Map<CEdge, Double> calculateAlphas(final Vector u) {
		HashMap<CEdge, Double> a = new HashMap<CEdge, Double>();
		for (CFace f : getFaces())
			triangleEnergyAndAlphas(u, f, a);
		return a;
	}
	
	
	
	
	public void conformalEnergy(final Vector u, final double[] E, final Vector G, final Matrix H, Map<CEdge, Double>... a) {
		Map<CEdge, Double> aMap = a.length != 0 ? a[0] : new HashMap<CEdge, Double>();
		aMap.clear();
		// Vertex Energy
		if (E != null) 
			E[0] = 0.0;
		if (G != null)
			G.zero();
		if (H != null)
			H.zero();
		for (final CVertex v : getVertices()) {
			if (!isVariable(v))
				continue;
			if (E != null)
				E[0] += v.getTheta() * u.get(v.getSolverIndex());
			if (G != null)
				G.add(v.getSolverIndex(), v.getTheta());
		}
		// Face Energy
		for (final CFace t : getFaces()) {
			final CEdge 
				e1 = t.getBoundaryEdge(),
				e2 = e1.getNextEdge(),
				e3 = e1.getPreviousEdge();
			final CVertex 
				v1 = e1.getTargetVertex(),
				v2 = e2.getTargetVertex(),
				v3 = e3.getTargetVertex();
			final int
				v1i = v1.getSolverIndex(),
				v2i = v2.getSolverIndex(),
				v3i = v3.getSolverIndex();
			final double e = triangleEnergyAndAlphas(u, t, aMap);
			if (E != null)
				E[0] += e;
			if (G != null) {
				if (isVariable(v1))
					G.add(v1i, -aMap.get(e3));
				if (isVariable(v2))
					G.add(v2i, -aMap.get(e1));
				if (isVariable(v3))
					G.add(v3i, -aMap.get(e2));
			}
			if (H != null) {
				final double[] 
				     cotE = {0, 0, 0},
				     cotV = {0, 0, 0};
				triangleHessian(u, t, cotE, cotV);
				// edge hessian
				if (isVariable(v1) && isVariable(v3)) {
					H.add(v1i, v3i, -cotE[0]);
					H.add(v3i, v1i, -cotE[0]);
				}
				if (isVariable(v2) && isVariable(v1)) {
					H.add(v2i, v1i, -cotE[1]);
					H.add(v1i, v2i, -cotE[1]);
				}
				if (isVariable(v3) && isVariable(v2)) {
					H.add(v2i, v3i, -cotE[2]);
					H.add(v3i, v2i, -cotE[2]);
				}
				// vertex hessian
				if (isVariable(v1))
					H.add(v1i, v1i, cotV[0]);
				if (isVariable(v2))
					H.add(v2i, v2i, cotV[1]);
				if (isVariable(v3))
					H.add(v3i, v3i, cotV[2]);
			}
		}
	}


	public Double evaluate(Vector x, Vector gradient, Matrix hessian) {
		double[] E = new double[]{0.0};
		conformalEnergy(x, E, gradient, hessian);
		return E[0];
	}


	public Double evaluate(Vector x, Vector gradient) {
		double[] E = new double[]{0.0};
		conformalEnergy(x, E, gradient, null);
		return E[0];
	}


	public Double evaluate(Vector x, Matrix hessian) {
		double[] E = new double[]{0.0};
		conformalEnergy(x, E, null, hessian);
		return E[0];
	}


	public Double evaluate(Vector x) {
		double[] E = new double[]{0.0};
		conformalEnergy(x, E, null, null);
		return E[0];
	}


	public Integer getDomainDimension() {
		return dim;
	}
	
	public void setDomainDimension(Integer dim) {
		this.dim = dim;
	}


	public boolean isTexCoordinatesValid() {
		return texCoordinatesValid;
	}


	public void setTexCoordinatesValid(boolean texCoordinatesValid) {
		this.texCoordinatesValid = texCoordinatesValid;
	}
	
}
