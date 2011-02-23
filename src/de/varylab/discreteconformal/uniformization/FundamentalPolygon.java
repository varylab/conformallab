package de.varylab.discreteconformal.uniformization;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

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
		return Tbig[1].compareTo(BigDecimal.valueOf(1E-5)) < 0;
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
	
}