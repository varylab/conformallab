package de.varylab.discreteconformal.functional;

import de.jreality.math.Rn;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.functional.DomainValue;
import de.jtem.halfedgetools.functional.Energy;
import de.jtem.halfedgetools.functional.Functional;
import de.jtem.halfedgetools.functional.FunctionalUtils;
import de.jtem.halfedgetools.functional.Gradient;
import de.jtem.halfedgetools.functional.Hessian;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Position;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Variable;

public class ElectrostaticSphereFunctional <
	V extends Vertex<V, E, F>,
	E extends Edge<V, E, F>,
	F extends Face<V, E, F>
> implements Functional<V, E, F> {

	private Variable<V, E>
		var = null;
	private Position<V>
		pos = null;
	
	public ElectrostaticSphereFunctional(Variable<V, E> variable, Position<V> position) {
		this.var = variable;
		this.pos = position;
	}
	
	@Override
	public <
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void evaluate(
		HDS hds,
		DomainValue x,
		Energy E,
		Gradient G,
		Hessian H
	) {
		if (E != null) { E.setZero(); }
		if (G != null) { G.setZero(); }
		int n = this.getDimension(hds);
		int nsq = n * n;
		final double EXP = 1.0;
		double[] vPos = new double[3];
		double[] wPos = new double[3];
		double[] dir = new double[3];
		double[] g = new double[3];
		for (V v : hds.getVertices()) {
			int vi = var.getVarIndex(v);
			this.getPosition(v, x, vPos);
			// electrostatic term
			for (V w : hds.getVertices()) {
				if (v == w) continue;
				int wi = var.getVarIndex(w);
				this.getPosition(w, x, wPos);
				double exp = !var.isVariable(v) || !var.isVariable(w) ? EXP : -EXP;
				Rn.subtract(dir, wPos, vPos);
				double dsq = Rn.innerProduct(dir, dir);
				if (E != null) {
					E.add(Math.pow(dsq, exp));
				}
				if (G != null) {
					Rn.times(g, exp * 2 * Math.pow(dsq, exp - 1), dir);
					if (var.isVariable(v)) {
						FunctionalUtils.subtractVectorFromGradient(G, vi * 3, g);
					}
					if (var.isVariable(w)) {
						FunctionalUtils.addVectorToGradient(G, wi * 3, g);
					}
				}
			}
			// spherical term
			if (!var.isVariable(v)) { continue; }
			double vdotv = Rn.innerProduct(vPos, vPos) - 1;
			if (E != null) {
				double ds = vdotv * vdotv;
				E.add(ds * nsq);
			}
			if (G != null) {
				double factor = 4 * nsq * vdotv;
				Rn.times(g, factor, vPos);
				FunctionalUtils.addVectorToGradient(G, vi * 3, g);
			}
		}
	}
	
	
	private void getPosition(V v, DomainValue x, double[] pos) {
		if (var.isVariable(v)) {
			int vi = this.var.getVarIndex(v);
			pos[0] = x.get(vi * 3 + 0);
			pos[1] = x.get(vi * 3 + 1);
			pos[2] = x.get(vi * 3 + 2);
		} else {
			double[] vPos = this.pos.getPosition(v);
			System.arraycopy(vPos, 0, pos, 0, 3);
		}
	}
	
	
	@Override
	public boolean hasGradient() {
		return true;
	}
	
	@Override
	public boolean hasHessian() {
		return false;
	}

	@Override
	public <
		HDS extends HalfEdgeDataStructure<V,E,F>
	> int getDimension(HDS hds) {
		int dim = 0;
		for (V v : hds.getVertices()) {
			if (var.isVariable(v)) { dim++;	}
		}
		return dim * 3;
	}

	@Override
	public <
		HDS extends HalfEdgeDataStructure<V, E, F>
	> int[][] getNonZeroPattern(HDS hds) {
		return new int[getDimension(hds)][0];
	}
	
}
