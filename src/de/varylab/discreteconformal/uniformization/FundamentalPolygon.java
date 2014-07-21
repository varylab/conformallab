package de.varylab.discreteconformal.uniformization;

import static de.varylab.discreteconformal.uniformization.FundamentalPolygonUtility.context;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import de.jreality.math.Pn;
import de.jtem.discretegroup.core.DiscreteGroup;
import de.jtem.discretegroup.core.DiscreteGroupConstraint;
import de.jtem.discretegroup.core.DiscreteGroupElement;
import de.varylab.discreteconformal.math.RnBig;

public class FundamentalPolygon {
	
	private List<FundamentalEdge>
		edgeList = new LinkedList<FundamentalEdge>();
	
	
	public FundamentalPolygon() {
	}
	
	public FundamentalPolygon(List<FundamentalEdge> edgeList) {
		this.edgeList.addAll(edgeList);
	}

	
	/**
	 * The number of edges in this polygon
	 * @return
	 */
	public int getLength() {
		return edgeList.size();
	}
	
	/**
	 * Calculates the number of vertices this fundamental polygon
	 * has. Identified vertices are counted only once.
	 * @return
	 */
	public int getNumVertices() {
		Set<FundamentalVertex> vSet = new HashSet<FundamentalVertex>();
		for (FundamentalEdge e : edgeList) {
			vSet.add(e.end);
			vSet.add(e.start);
		}
		return vSet.size();
	}
	
	/**
	 * The positions of vertices in cyclic order. Here every vertex has a unique 
	 * position in hyperbolic 2-space embedded in RP3 (x,y,0,w)
	 * @return
	 */
	public List<BigDecimal[]> getVertexPositions() {
		List<BigDecimal[]> result = new ArrayList<BigDecimal[]>(edgeList.size());
		for (FundamentalEdge e : edgeList) {
			result.add(e.startPosition);
		}
		return result;
	}
	
	
	/**
	 * Calculates the relation around a vertex of the polygon
	 * @param root
	 * @return
	 */
	public boolean checkRelation() {
		FundamentalEdge start = edgeList.get(0);
		FundamentalEdge active = start;
		BigDecimal[] Tbig = new BigDecimal[16];
		RnBig.setIdentityMatrix(Tbig);
		do {
			RnBig.times(Tbig, Tbig, active.motionBig, FundamentalPolygonUtility.context);
			active = active.partner.nextEdge;
		} while (active != start);
		System.out.println("\nRelation: \n" + RnBig.matrixToString(Tbig));
		BigDecimal[] I = new BigDecimal[16];
		RnBig.setIdentityMatrix(I);
		for (int j = 0; j < I.length; j++) {
			if (Math.abs(Tbig[j].doubleValue() - I[j].doubleValue()) > 1E-8) {
				return false;
			}
		}
		return true;
	}
	

	/**
	 * Returns the edges of this polygon in the order of
	 * appearance
	 * @return
	 */
	public List<FundamentalEdge> getEdges() {
		return Collections.unmodifiableList(edgeList);
	}
	
	/**
	 * Returns all vertices of this polygon. Identified vertices are not
	 * returned.
	 * @return
	 */
	public Set<FundamentalVertex> getVertices() {
		Set<FundamentalVertex> vSet = new HashSet<FundamentalVertex>();
		for (FundamentalEdge e : edgeList) {
			vSet.add(e.start);
			vSet.add(e.end);
		}
		return vSet;
	}
	
	/**
	 * Calculates the valence of the given vertex in the universal cover
	 * @param v
	 * @return
	 */
	public int getValence(FundamentalVertex v) {
		int valence = 0;
		for (FundamentalEdge e : edgeList) {
			if (e.start == v) valence++;
		}
		return valence / 2;
	}
	
	
	
	public FundamentalVertex getMaxValenceVertex() {
		int max = 0;
		FundamentalVertex maxV = null;
		for (FundamentalVertex v : getVertices()) {
			int valence = getValence(v);
			if (valence > max) {
				max = valence;
				maxV = v;
			}
		}
		assert maxV != null;
		return maxV;
	}
	
	
	public void normalizeEdgeList() {
		List<FundamentalEdge> normalizedList = new LinkedList<FundamentalEdge>();
		FundamentalEdge a = edgeList.get(0);
		do {
			normalizedList.add(a);
			a = a.nextEdge;
		} while (a != edgeList.get(0));
		edgeList = normalizedList;
	}
	
	
	@Override
	public String toString() {
		StringBuffer sb = new StringBuffer();
		sb.append("Fundamental Polygon edges: " + edgeList.size());
		sb.append("\n");
		for (FundamentalEdge fe : edgeList) {
			sb.append(fe + ": ");
			sb.append(fe.start.index + " <-> " + fe.end.index);
			sb.append(" -> " + fe.partner);
			sb.append(";  " + fe.prevEdge.index + ":-:" + fe.nextEdge.index);
			sb.append("\n");
		}
		sb.append("genus: " + getGenus() + "\n");
		sb.append("---------------------------------");
		return sb.toString();
	}
	
	
	public int getGenus() {
		int X = getNumVertices() - edgeList.size() / 2 + 1;
		return 1 - X/2;
	}
	
	
	public DiscreteGroup getDiscreteGroup() {
		int dim = edgeList.size();
		DiscreteGroupElement[] gArr = new DiscreteGroupElement[dim];
		DiscreteGroup G = new DiscreteGroup();
		for (FundamentalEdge e : edgeList) {
			DiscreteGroupElement g = new DiscreteGroupElement();
			g.setMatrix(e.motion);
			gArr[e.index] = g;
		}
		G.setGenerators(gArr);
		if (getGenus() == 0) {
			G.setMetric(Pn.ELLIPTIC);
		} else if (getGenus() == 1) {
			G.setMetric(Pn.EUCLIDEAN);
		} else {
			G.setMetric(Pn.HYPERBOLIC);
		}
		G.calculateGenerators();
		assert !G.isFree() : "a fuchsian uniformization group has at least one relation";
		assert !G.isFinite() : "these groups should tesselate hyperbolic space and hence are not finite";
		return G;
	}

	
	public List<DiscreteGroupElement> createGoupElements(final int maxDrawElements, final double maxDrawDistance) {
		DiscreteGroupConstraint constraint = new DiscreteGroupConstraint() {
			@Override
			public void update() {
			}
			@Override
			public void setMaxNumberElements(int arg0) {
			}
			@Override
			public int getMaxNumberElements() {
				return maxDrawElements;
			}
			@Override
			public boolean acceptElement(DiscreteGroupElement s) {
				double mean = 0.0;
				for (FundamentalEdge e : getEdges()) {
					BigDecimal[] p = e.startPosition;
					BigDecimal[] T = RnBig.toBig(null, s.getMatrix());
					p = RnBig.matrixTimesVector(null, T, p, context);
					double[] pd = RnBig.toDouble(null, p);
					mean += Math.sqrt(Pn.norm(pd, Pn.HYPERBOLIC));
				}
				mean /= getLength();
				return mean < maxDrawDistance;
			}
		};
		DiscreteGroup G = getDiscreteGroup();
		G.setConstraint(constraint);
		G.generateElements();
		DiscreteGroupElement[] elements = G.getElementList();
		List<DiscreteGroupElement> elementList = Arrays.asList(elements);
		return elementList;
	}
		
}