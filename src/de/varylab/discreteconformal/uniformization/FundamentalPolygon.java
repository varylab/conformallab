package de.varylab.discreteconformal.uniformization;

import static de.jreality.math.Pn.HYPERBOLIC;

import java.math.BigDecimal;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import de.jreality.math.Matrix;
import de.jreality.math.Pn;
import de.varylab.discreteconformal.math.PnBig;
import de.varylab.discreteconformal.math.RnBig;

public class FundamentalPolygon {
	
	public List<FundamentalEdge>
		edgeList = new LinkedList<FundamentalEdge>();
	
	public int getLength() {
		return edgeList.size();
	}
	
	
	public int getNumVertices() {
		Set<FundamentalVertex> vSet = new HashSet<FundamentalVertex>();
		for (FundamentalEdge e : edgeList) {
			vSet.add(e.end);
			vSet.add(e.start);
		}
		return vSet.size();
	}
	
	/**
	 * Transforms a given position around the polygon by applying
	 * consecutive motions. The order of the result is clockwise
	 * @param root
	 * @return
	 */
	public List<BigDecimal[]> getOrbit(double[] root) {
		BigDecimal[] pos1 = RnBig.toBig(null, root);
		BigDecimal[] pos2 = RnBig.toBig(null, root);
		
		Map<FundamentalEdge, BigDecimal[]> posMap = new HashMap<FundamentalEdge, BigDecimal[]>();
		FundamentalEdge start1 = edgeList.get(0);
		FundamentalEdge active1 = start1;
		FundamentalEdge start2 = edgeList.get(edgeList.size() - 1);
		FundamentalEdge active2 = start2; 
		while (posMap.keySet().size() < edgeList.size()) {
			posMap.put(active2.nextEdge, pos2.clone());
			posMap.put(active1, pos1.clone());
			RnBig.matrixTimesVector(pos1, active1.partner.motionBig, pos1, UniformizationUtility.context);
			RnBig.matrixTimesVector(pos2, active2.partner.motionBig, pos2, UniformizationUtility.context);
			PnBig.normalize(pos1, pos1, HYPERBOLIC, UniformizationUtility.context);
			PnBig.normalize(pos2, pos2, HYPERBOLIC, UniformizationUtility.context);
			active1 = active1.partner.nextEdge;
			active2 = active2.partner.prevEdge;
		};
		List<BigDecimal[]> result = new LinkedList<BigDecimal[]>();
		for (FundamentalEdge e : edgeList) {
			BigDecimal[] posBig = posMap.get(e);
			result.add(posBig);
		}
		return result;
	}
	
	
	/**
	 * Calculates the relation around a vertex of the polygon
	 * @param root
	 * @return
	 */
	public List<double[]> getDualOrbit(double[] root) {
		List<double[]> result = new LinkedList<double[]>();
		FundamentalEdge start = edgeList.get(0);
		FundamentalEdge active = start;
		Matrix T = new Matrix();
		BigDecimal[] Tbig = new BigDecimal[16];
		RnBig.setIdentityMatrix(Tbig);
		System.out.println("------------------------ orbit calculation");
		do {
			System.out.print(active.index + ", ");
			double[] pos = T.multiplyVector(root);
			Pn.normalize(pos, pos, HYPERBOLIC);
			result.add(pos);
			// apply in the opposite order to get the relation
			T.multiplyOnRight(active.motion);
			RnBig.times(Tbig, Tbig, active.motionBig, UniformizationUtility.context);
			active = active.partner.nextEdge;
		} while (active != start);
		System.out.println("\nDual orbit transform: \n" + T);
		System.out.println("\nBig dual orbit transform: \n" + RnBig.matrixToString(Tbig));
		return result;
	}
	
	
	/**
	 * Constructs a minimal fundamental polygon from a general polygon
	 * The result will have only one vertex and 4g edges
	 * @return
	 */
	public FundamentalPolygon getMinimal() {
		// canonical polygon construction --------------
		System.out.println("Constructing minimal polygon...");
//		FundamentalVertex fRoot = edgeList.get(0).start;
//		List<FundamentalEdge> newPoly = new LinkedList<FundamentalEdge>();
//		Set<FundamentalEdge> contracted = new TreeSet<FundamentalEdge>();
//		for (FundamentalEdge e : edgeList) {
//			if (contracted.contains(e)) {
//				continue;
//			}
//			FundamentalVertex badVertex = e.end;
//			if (badVertex != fRoot) { // edge contraction
//				System.out.println("contracting edge " + e + ", #source " + e.sourceEdgeCount);
//				contracted.add(e);
//				contracted.add(e.partner);
//				for (FundamentalEdge fe : edgeList) {
//					if (fe.start == badVertex) {
//						fe.start = fRoot;
//					}
//					if (fe.end == badVertex) {
//						fe.end = fRoot;
//					}
//				}
//			} else {
//				newPoly.add(e);
//			}
//		}
//		edgeList = newPoly;
//		// linkage
//		for (int i = 0; i < newPoly.size(); i++) {
//			FundamentalEdge prev = newPoly.get(i - 1 < 0 ? newPoly.size() - 1 : i - 1);
//			FundamentalEdge next = newPoly.get((i + 1) % newPoly.size());
//			FundamentalEdge act = newPoly.get(i);
//			prev.nextEdge = act;
//			act.index = i;
//			act.prevEdge = prev;
//			act.nextEdge = next;
//			next.prevEdge = act;
//		}
//		FundamentalPolygon result = new FundamentalPolygon();
//		result.edgeList = newPoly;
//		return result;
		
		
		FundamentalVertex rootVertex = edgeList.get(0).start;
		List<FundamentalEdge> newPoly = new LinkedList<FundamentalEdge>(edgeList);
		
		FundamentalEdge cEdge = findShortestContractableEdge(newPoly, rootVertex);
		while (cEdge != null) {
			System.out.println("contracting edge " + cEdge + " #src " + cEdge.sourceEdgeCount);
			FundamentalEdge prev1 = cEdge.prevEdge;
			FundamentalEdge next1 = cEdge.nextEdge;
			FundamentalEdge prev2 = cEdge.partner.prevEdge;
			FundamentalEdge next2 = cEdge.partner.nextEdge;
			prev1.nextEdge = next1;
			next1.prevEdge = prev1;
			prev2.nextEdge = next2;
			next2.prevEdge = prev2;
			next1.start = rootVertex;
			prev2.end = rootVertex;
			newPoly.remove(cEdge);
			newPoly.remove(cEdge.partner);
			FundamentalVertex v = cEdge.end;
			for (FundamentalEdge fe : edgeList) { // move vertex
				if (fe.start == v) {
					fe.start = rootVertex;
				}
				if (fe.end == v) {
					fe.end = rootVertex;
				}
			}
			cEdge = findShortestContractableEdge(newPoly, rootVertex);
		}
		int index = 0; // set edge indices
		for (FundamentalEdge e : newPoly) {
			e.index = index++;
		}
		FundamentalPolygon result = new FundamentalPolygon();
		result.edgeList = newPoly;
		return result;
	}
	
	
	/**
	 * Finds the fundamental edge with the smallest number of source edges 
	 * in the given polygon. This edge will have root vertex as start vertex
	 * @param root 
	 * @return the contractable fundamental edge with the smallest number of source
	 * edges of null if the polygon does not contain any contractable edges.
	 */
	private FundamentalEdge findShortestContractableEdge(List<FundamentalEdge> edges, FundamentalVertex root) {
		FundamentalEdge result = null;
		int length = Integer.MAX_VALUE;
		for (FundamentalEdge e : edges) {
			if (e.start == root && e.end != root) {
				if (e.sourceEdgeCount < length) {
					length = e.sourceEdgeCount;
					result = e;
				}
			}
		}
		return result;
	}
	
	
	
	
	/**
	 * Constructs a canonical fundamental polygon. The order of the 
	 * edges is aba-1b-1.... This polygon can be non-convex
	 * @return
	 */
	public FundamentalPolygon getNaiveCanonical() {
		FundamentalPolygon min = getMinimal();
		FundamentalEdge a = min.edgeList.get(0);
		Set<FundamentalEdge> normalized = new HashSet<FundamentalEdge>();
		int g = getGenus();
		while (normalized.size() < 4*g) {
			if (normalized.contains(a)) {
				a = a.nextEdge;
				System.out.println("skipping normalized edge " + a);
				continue;
			}
			System.out.println("Collecting handle ---------------");
			FundamentalEdge b = findLinkedEdge(a);
			bringTogether(a, b);
			getDualOrbit(new double[] {0,0,0,1});
			enumerateFrom(a);
			bringTogether(b, a.partner);
			getDualOrbit(new double[] {0,0,0,1});
			enumerateFrom(a);
			bringTogether(a.partner, b.partner);
			getDualOrbit(new double[] {0,0,0,1});
			enumerateFrom(a);
			normalized.add(a);
			normalized.add(b);
			normalized.add(a.partner);
			normalized.add(b.partner);
			a = b.partner.nextEdge;
		};
		// assemble the polygon list
		List<FundamentalEdge> canList = new LinkedList<FundamentalEdge>();
		int index = 0;
		FundamentalEdge first = a;
		do {
			a.index = index++;
			canList.add(a);
			a = a.nextEdge;
		} while (a != first);
		FundamentalPolygon r = new FundamentalPolygon();
		r.edgeList = canList;
		return r;
	}
	
	
	/**
	 * Tries to use a heuristics to minimize the group operations
	 * whine transforming the polygon. This uses the notion of distance
	 * or cost of a linked pair
	 * @return
	 */
	public FundamentalPolygon getFastCanonical() {
		FundamentalPolygon P = getMinimal();
		FundamentalEdge[] link = findMinimalLinkableEdges(P);
		System.out.println("Constructing fast canonical polygon...");
		if (link == null)  {
			System.out.println("polygon is canonical");
			return P;
		}
		FundamentalEdge lastHandle = null;
		while (link != null) {
			FundamentalEdge a = link[0];
			FundamentalEdge b = link[1];
			System.out.println("Collecting handle --- " + a + " :: " + b);
			bringTogether(a, b);
			getDualOrbit(new double[] {0,0,0,1});
			enumerateFrom(a);
			bringTogether(b, a.partner);
			getDualOrbit(new double[] {0,0,0,1});
			enumerateFrom(a);
			bringTogether(a.partner, b.partner);
			getDualOrbit(new double[] {0,0,0,1});
			enumerateFrom(a);
			link = findMinimalLinkableEdges(P);
			lastHandle = a;
		};
		if (lastHandle == null) {
			lastHandle = P.edgeList.get(0);
		}
		// assemble the polygon list
		List<FundamentalEdge> canList = new LinkedList<FundamentalEdge>();
		int index = 0;
		FundamentalEdge a = lastHandle;
		do {
			a.index = index++;
			canList.add(a);
			a = a.nextEdge;
		} while (a != lastHandle);
		FundamentalPolygon r = new FundamentalPolygon();
		r.edgeList = canList;
		System.out.println("Result of normalization ");
		enumerateFrom(lastHandle);
		return r;
	}
	
	
	
	
	/**
	 * Enumerates the polygon starting from the given edge clockwise
	 * @param e
	 */
	private void enumerateFrom(FundamentalEdge e) {
		System.out.println("Polygon Enumeration -----------------------------");
		FundamentalEdge act = e;
		do {
			System.out.print(act + ": ");
			System.out.print(act.start.index + " <-> " + act.end.index);
			System.out.print(" -> " + act.partner);
			System.out.print("\n");
			act = act.nextEdge;
		} while (act != e);
	}
	
	
	/**
	 * Finds some linked edge
	 * @param a
	 * @return
	 */
	private FundamentalEdge findLinkedEdge(FundamentalEdge a) {
		Set<FundamentalEdge> checkSet = new TreeSet<FundamentalEdge>();
		for (FundamentalEdge e = a.nextEdge; e != a.partner; e = e.nextEdge) {
			checkSet.add(e);
		}
		FundamentalEdge b = null;
		for (FundamentalEdge e : checkSet) {
			if (!checkSet.contains(e.partner)) {
				b = e;
				break;
			}
		}
		return b;
	}
	
	
	/**
	 * Finds the linked pair with shortest distance that is not already minimally linked
	 * @param P
	 * @return a linked pair or null if no pair is linkable
	 */
	private FundamentalEdge[] findMinimalLinkableEdges(FundamentalPolygon P) {
		// exclude handles
		Set<FundamentalEdge> handles = new HashSet<FundamentalEdge>();
		for (FundamentalEdge e : P.edgeList) {
			FundamentalEdge l = findMinimalLinkedEdge(e);
			int cost = calculateLinkedPairCost(e, l);
			if (cost == 3) { // is a handle
				handles.add(e);
				handles.add(l);
				handles.add(e.partner);
				handles.add(l.partner);
			}
		}
		
		// collect linkable edges
		Set<FundamentalEdge> linkable = new HashSet<FundamentalEdge>(P.edgeList);
		linkable.removeAll(handles);
		if (linkable.isEmpty()) return null;
		
		FundamentalEdge[] result = new FundamentalEdge[2];
		int min = Integer.MAX_VALUE;
		for (FundamentalEdge e : linkable) {
			FundamentalEdge l = findMinimalLinkedEdge(e);
			int cost = calculateLinkedPairCost(e, l);
			if (cost < min) { 
				min = cost;
				result[0] = e;
				result[1] = l;
			}
		}
		return result;
	}
	
	
	private int calculateLinkedPairCost(FundamentalEdge a, FundamentalEdge b) {
		int cost = 0;
		cost += calculateDistanceInPolygon(a, b);
		cost += calculateDistanceInPolygon(b, a.partner);
		cost += calculateDistanceInPolygon(a.partner, b.partner);
		return cost;
	}
	
	
	
	/**
	 * Finds an edge whose partner is close to the partner of a
	 * @param a
	 * @return
	 */
	private FundamentalEdge findMinimalLinkedEdge(FundamentalEdge a) {
		Set<FundamentalEdge> checkSet = new TreeSet<FundamentalEdge>();
		for (FundamentalEdge e = a.nextEdge; e != a.partner; e = e.nextEdge) {
			checkSet.add(e);
		}
		FundamentalEdge b = null;
		int min = Integer.MAX_VALUE;
		for (FundamentalEdge e : checkSet) {
			if (!checkSet.contains(e.partner)) {
				int dist = calculateDistanceInPolygon(a.partner, e.partner);
				if (dist < min) {
					min = dist;
					b = e;	
				}
			}
		}
		return b;
	}
	
	
	/**
	 * Calculates the edge distance of two edges when going in next direction
	 * @param e1
	 * @param e2
	 * @return
	 */
	private int calculateDistanceInPolygon(FundamentalEdge e1, FundamentalEdge e2) {
		int i = 0;
		while (e1 != e2) {
			i++;
			e1 = e1.nextEdge;
		}
		return i;
	}
	
	
	/**
	 * Performs a cut and paste operation on the given edges
	 * such that they are next to each other afterwards
	 * @param a
	 * @param b
	 */
	private void bringTogether(FundamentalEdge a, FundamentalEdge b) {
		System.out.println("bring together " + a + " - " + b);
		if (a.nextEdge == b) return; // already together
		Set<FundamentalEdge> cSet = new TreeSet<FundamentalEdge>();
		for (FundamentalEdge e = a.nextEdge; e != b; e = e.nextEdge) {
			cSet.add(e);
		}
		FundamentalEdge aiPrev = a.partner.prevEdge;
		FundamentalEdge c1 = a.nextEdge;
		FundamentalEdge cn = b.prevEdge;
		
		Matrix A = a.motion;
		Matrix Ainv = a.partner.motion;
		BigDecimal[] ABig = a.motionBig;
		BigDecimal[] ABiginv = a.partner.motionBig;
		System.out.println("A = \n" + A);
		
		for (FundamentalEdge c : cSet) {
			System.out.println(c.index + " = " + c.index + " " + a.partner.index);
			c.motion.multiplyOnLeft(Ainv);
			RnBig.times(c.motionBig, ABiginv, c.motionBig, UniformizationUtility.context);
			System.out.println(c.partner.index + " = " + a.index + " " + c.partner.index);
			c.partner.motion.multiplyOnRight(A);
			RnBig.times(c.partner.motionBig, c.partner.motionBig, ABig, UniformizationUtility.context);
		}
		// move first connection
		aiPrev.nextEdge = c1;
		c1.prevEdge = aiPrev;
		cn.nextEdge = a.partner;
		a.partner.prevEdge = cn;
		// bring together
		a.nextEdge = b;
		b.prevEdge = a;
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