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
	public FundamentalPolygon getCanonical() {
		FundamentalPolygon min = getMinimal();
		FundamentalEdge a = min.edgeList.get(0);
		int g = getGenus();
		for (int i = 0; i < g; i++) {
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