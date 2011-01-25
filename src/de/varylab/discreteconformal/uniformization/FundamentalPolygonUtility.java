package de.varylab.discreteconformal.uniformization;

import static de.jreality.math.Pn.HYPERBOLIC;
import static de.jtem.halfedge.util.HalfEdgeUtils.incomingEdges;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import de.jreality.math.Matrix;
import de.jreality.math.P2;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.math.P2Big;
import de.varylab.discreteconformal.math.RnBig;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;

public class FundamentalPolygonUtility {

	public static MathContext 
		context = new MathContext(50);
	

	public static FundamentalPolygon copyPolygon(
		FundamentalPolygon P,
		Map<FundamentalVertex, FundamentalVertex> oldNewVertexMap,
		Map<FundamentalEdge, FundamentalEdge> oldNewEdgeMap
	) {
		if (oldNewVertexMap == null) {
			oldNewVertexMap = new HashMap<FundamentalVertex, FundamentalVertex>();
		}
		if (oldNewEdgeMap == null) {
			oldNewEdgeMap = new HashMap<FundamentalEdge, FundamentalEdge>();
		}
		List<FundamentalEdge> newEdges = new ArrayList<FundamentalEdge>();
		// new nodes
		for (FundamentalEdge e : P.getEdges()) {
			FundamentalEdge newE = new FundamentalEdge(e.index);
			oldNewEdgeMap.put(e, newE);
			newE.motion = new Matrix(e.motion);
			newE.motionBig = e.motionBig.clone();
			newE.sourceEdgeCount = e.sourceEdgeCount;
			newE.startPosition = e.startPosition.clone();
			FundamentalVertex start = oldNewVertexMap.get(e.start);
			if (start == null) {
				start = new FundamentalVertex(e.start.index);
				oldNewVertexMap.put(e.start, start);
			}
			FundamentalVertex end = oldNewVertexMap.get(e.end);
			if (end == null) {
				end = new FundamentalVertex(e.end.index);
				oldNewVertexMap.put(e.end, end);
			}
			newE.start = start;
			newE.end = end;
			newEdges.add(newE);
		}
		// linkage
		for (FundamentalEdge e : P.getEdges()) {
			FundamentalEdge newE = oldNewEdgeMap.get(e);
			FundamentalEdge next = oldNewEdgeMap.get(e.nextEdge);
			FundamentalEdge prev = oldNewEdgeMap.get(e.prevEdge);
			FundamentalEdge partner = oldNewEdgeMap.get(e.partner);
			newE.nextEdge = next;
			newE.prevEdge = prev;
			newE.partner = partner;
		}
		return new FundamentalPolygon(newEdges);
	}
	
	public static FundamentalPolygon copyPolygon(FundamentalPolygon P) {
		return copyPolygon(P, null, null);
	}
	
	/**
	 * Enumerates the polygon starting from the given edge clockwise
	 * @param e
	 */
	public static void printPolygonFrom(FundamentalEdge e) {
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
	 * Performs a cut and paste operation on the given edges
	 * such that they are next to each other afterwards. It uses
	 * the transformation a to move the cut part
	 * @param a
	 * @param b
	 * @return the number of matrix products performed
	 */
	private static int bringTogetherViaA(FundamentalEdge a, FundamentalEdge b) {
		System.out.println("bring together via a " + a + " - " + b);
		int cost = 0;
		if (a.nextEdge == b) { // already together
			return cost; 
		}
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
		// move start and end
		RnBig.matrixTimesVector(a.partner.startPosition, ABiginv, b.startPosition, context);
		// move c edges
		for (FundamentalEdge c : cSet) {
			System.out.println(c.index + " = " + c.index + " " + a.partner.index);
			c.motion.multiplyOnLeft(Ainv);
			RnBig.times(c.motionBig, ABiginv, c.motionBig, FundamentalPolygonUtility.context);
			System.out.println(c.partner.index + " = " + a.index + " " + c.partner.index);
			c.partner.motion.multiplyOnRight(A);
			RnBig.times(c.partner.motionBig, c.partner.motionBig, ABig, FundamentalPolygonUtility.context);
			RnBig.matrixTimesVector(c.startPosition, ABiginv, c.startPosition, context);
			cost += 2;
		}
		// move first connection
		aiPrev.nextEdge = c1;
		c1.prevEdge = aiPrev;
		cn.nextEdge = a.partner;
		a.partner.prevEdge = cn;
		// bring together
		a.nextEdge = b;
		b.prevEdge = a;
		return cost; 
	}
	
	
	/**
	 * Performs a cut and paste operation on the given edges
	 * such that they are next to each other afterwards. It uses
	 * the transformation b to move the cut part
	 * @param a
	 * @param b
	 * @return the number of matrix products performed
	 */
	private static int bringTogetherViaB(FundamentalEdge a, FundamentalEdge b) {
		System.out.println("bring together via b " + a + " - " + b);
		int cost = 0;
		if (a.nextEdge == b) { // already together
			return cost; 
		}
		Set<FundamentalEdge> cSet = new TreeSet<FundamentalEdge>();
		for (FundamentalEdge e = a.nextEdge; e != b; e = e.nextEdge) {
			cSet.add(e);
		}
		FundamentalEdge biNext = b.partner.nextEdge;
		FundamentalEdge c1 = a.nextEdge;
		FundamentalEdge cn = b.prevEdge;
		
		Matrix B = b.motion;
		Matrix Binv = b.partner.motion;
		BigDecimal[] BBig = b.motionBig;
		BigDecimal[] BBiginv = b.partner.motionBig;
		
		b.startPosition = c1.startPosition.clone();
		for (FundamentalEdge c : cSet) {
			System.out.println(c.index + " = " + c.index + " " + a.partner.index);
			c.motion.multiplyOnLeft(Binv);
			RnBig.times(c.motionBig, BBiginv, c.motionBig, FundamentalPolygonUtility.context);
			System.out.println(c.partner.index + " = " + a.index + " " + c.partner.index);
			c.partner.motion.multiplyOnRight(B);
			RnBig.times(c.partner.motionBig, c.partner.motionBig, BBig, FundamentalPolygonUtility.context);
			RnBig.matrixTimesVector(c.startPosition, BBiginv, c.startPosition, context);
			cost += 2;
		}
		// move first connection
		biNext.prevEdge = cn;
		cn.nextEdge = biNext;
		c1.prevEdge = b.partner;
		b.partner.nextEdge = c1;
		// bring together
		a.nextEdge = b;
		b.prevEdge = a;
		return cost; 
	}
	
	
	
	/**
	 * Brings together linked edges a and b to form a canonical handle.
	 * The cost of this method is 6C + 4D + 2E
	 * @param a
	 * @param b
	 * @return the cost of the operation
	 */
	private static int gatherHandle(FundamentalEdge a, FundamentalEdge b) {
		int cost = 0;
		cost += bringTogetherViaA(a, b);
		cost += bringTogetherViaA(b, a.partner);
		cost += bringTogetherViaA(a.partner, b.partner);
		return cost;
	}
	
	/**
	 * Brings together linked edges a and b to form a canonical handle.
	 * The cost of this method is 2C + 4D + 2E
	 * @param a
	 * @param b
	 * @return the cost of the operation
	 */
	private static int gatherHandleFast(FundamentalEdge a, FundamentalEdge b) {
		int cost = 0;
		cost += bringTogetherViaB(a, b);
		cost += bringTogetherViaA(b, a.partner);
		cost += bringTogetherViaA(a.partner, b.partner);
		return cost;
	}

	
	
	/**
	 * Finds some linked edge
	 * @param a
	 * @return
	 */
	private static FundamentalEdge findLinkedEdge(FundamentalEdge a) {
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
	private static FundamentalEdge[] findMinimalLinkableEdges(FundamentalPolygon P) {
		// exclude handles
		List<FundamentalEdge> edgeList = P.getEdges();
		Set<FundamentalEdge> handles = new HashSet<FundamentalEdge>();
		for (FundamentalEdge e : edgeList) {
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
		Set<FundamentalEdge> linkable = new HashSet<FundamentalEdge>(edgeList);
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
	
	
	private static int calculateLinkedPairCost(FundamentalEdge a, FundamentalEdge b) {
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
	private static FundamentalEdge findMinimalLinkedEdge(FundamentalEdge a) {
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
	private static int calculateDistanceInPolygon(FundamentalEdge e1, FundamentalEdge e2) {
		int i = 0;
		while (e1 != e2) {
			i++;
			e1 = e1.nextEdge;
		}
		return i;
	}
	
	
	
	
	/**
	 * Finds the fundamental edge with the smallest number of source edges 
	 * in the given polygon. This edge will have root vertex as start vertex
	 * @param root 
	 * @return the contractable fundamental edge with the smallest number of source
	 * edges of null if the polygon does not contain any contractable edges.
	 */
	private static FundamentalEdge findShortestContractableEdge(List<FundamentalEdge> P, FundamentalVertex root) {
		FundamentalEdge result = null;
		int length = Integer.MAX_VALUE;
		for (FundamentalEdge e : P) {
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
	 * Tries to use a heuristics to minimize the group operations
	 * while transforming the polygon. This uses the notion of distance
	 * or cost of a linked pair.
	 * Here maybe no edge is preserved
	 * @param P a minimal polygon
	 * @return
	 */
	public static FundamentalPolygon canonicalize(FundamentalPolygon P) {
		System.out.println("Vertices P: " + P.getVertices());
		FundamentalPolygon R = copyPolygon(P);// minimize(P, root);
		System.out.println("Vertices R: " + R.getVertices());
		FundamentalEdge[] link = findMinimalLinkableEdges(R);
		if (link == null)  {
			System.out.println("polygon is canonical");
			return R;
		}
		FundamentalEdge lastHandle = null;
		int cost = 0;
		while (link != null) {
			FundamentalEdge a = link[0];
			FundamentalEdge b = link[1];
			System.out.println("Collecting handle --- " + a + " :: " + b);
			cost += gatherHandleFast(a, b);
			link = findMinimalLinkableEdges(R);
			lastHandle = a;
		};
		if (lastHandle == null) {
			lastHandle = R.getEdges().get(0);
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
		R = new FundamentalPolygon(canList);
		System.out.println("Vertices R: " + R.getVertices());
		return R;
	}
	
	
	/**
	 * Constructs a canonical fundamental polygon. The order of the 
	 * edges is aba-1b-1.... This polygon can be non-convex 
	 * The first edge of the polygon is preserved, thus one can proceed drawing 
	 * the fundamental polygon by using the start vertex of the first edge
	 * @return
	 */
	public static FundamentalPolygon canonicalizeNaive(FundamentalPolygon P, FundamentalVertex root) {
		FundamentalPolygon R = minimize(P, root); // we need the polygon to have only one vertex
		FundamentalEdge a = R.getEdges().get(0);
		Set<FundamentalEdge> normalized = new HashSet<FundamentalEdge>();
		int g = R.getGenus();
		int cost = 0;
		while (normalized.size() < 4*g) {
			if (normalized.contains(a)) {
				a = a.nextEdge;
				System.out.println("skipping normalized edge " + a);
				continue;
			}
			System.out.println("Collecting handle ---------------");
			FundamentalEdge b = findLinkedEdge(a);
			cost += gatherHandle(a, b);
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
		R = new FundamentalPolygon(canList);
		System.out.println("Result of normalization cost = " + cost);
		R.checkRelation();
		return R;
	}
	
	
	/**
	 * Constructs a minimal fundamental polygon from a general polygon
	 * The result will have only one vertex and 4g edges
	 * @return
	 */
	public static FundamentalPolygon minimize(FundamentalPolygon P, FundamentalVertex root) {
		Map<FundamentalVertex, FundamentalVertex> vertMap = new HashMap<FundamentalVertex, FundamentalVertex>();
		System.out.println("Vertices P: " + P.getVertices());
		FundamentalPolygon R = copyPolygon(P, vertMap, null);
		System.out.println("Vertices R: " + R.getVertices());
		FundamentalVertex newRoot = vertMap.get(root);
		
		List<FundamentalEdge> edgeList = new ArrayList<FundamentalEdge>(R.getEdges());
		FundamentalEdge cEdge = findShortestContractableEdge(edgeList, newRoot);
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
			next1.start = newRoot;
			prev2.end = newRoot;
			edgeList.remove(cEdge);
			edgeList.remove(cEdge.partner);
			FundamentalVertex v = cEdge.end;
			for (FundamentalEdge fe : edgeList) { // move vertex
				if (fe.start == v) {
					fe.start = newRoot;
				}
				if (fe.end == v) {
					fe.end = newRoot;
				}
			}
			// move vertex positions by iterating around the root vertex
			FundamentalEdge e = next1;
			e.startPosition = cEdge.startPosition;
			while (e != prev2.partner) {
				BigDecimal[] rootPos = e.startPosition;
				BigDecimal[] motion = e.partner.motionBig;
			 	BigDecimal[] newStart = RnBig.matrixTimesVector(null, motion, rootPos, context); 
				e.partner.nextEdge.startPosition = newStart;
				e = e.partner.nextEdge;
			}
			
			cEdge = findShortestContractableEdge(edgeList, newRoot);
		}
		int index = 0; // set edge indices
		for (FundamentalEdge e : edgeList) {
			e.index = index++;
		}
		return new FundamentalPolygon(edgeList);
	}
	
	
	
	
	/**
	 * Uses the information from cutInfo to create the generators and 
	 * vertex positions of a fundamental polygon
	 * @param cutInfo
	 * @return
	 */
	public static FundamentalPolygon constructFundamentalPolygon(
		CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo
	) {
		List<FundamentalEdge> edgeList = new ArrayList<FundamentalEdge>();
		// find max valence branch
		Set<CoVertex> branchSet = cutInfo.getBranchSet();
		int maxBranch = 0; 
		CoVertex root = branchSet.iterator().next();
		for (CoVertex branch : branchSet) {
			int branchNumber = cutInfo.getCopies(branch).size();
			if (branchNumber > maxBranch) {
				maxBranch = branchNumber;
				root = branch;
			}
		}
		
		// find the start boundary edge
		CoEdge rootEdge = null;
		for (CoEdge e : incomingEdges(root)) {
			if (e.getOppositeEdge().getLeftFace() == null) {
				rootEdge = e.getOppositeEdge();
				break;
			}
		}
		assert rootEdge != null;
		
		// create fundamental vertices and partition the branch set
		Map<CoVertex, FundamentalVertex> branchMap = new HashMap<CoVertex, FundamentalVertex>();
		Set<CoVertex> tmpSet = new HashSet<CoVertex>(branchSet);
		int index = 0;
		for (CoVertex v : branchSet) {
			if (!tmpSet.contains(v)) {
				continue;
			}
			Set<CoVertex> copies = cutInfo.getCopies(v);
			FundamentalVertex fv = new FundamentalVertex(index++);
			for (CoVertex copy : copies) {
				branchMap.put(copy, fv);
			}
			tmpSet.removeAll(copies);
		}
		
		Set<CoEdge> visited = new HashSet<CoEdge>();
		CoEdge eActive = rootEdge;
		CoEdge firstOfSegment = rootEdge.getOppositeEdge();
		
		FundamentalVertex lastFunV = branchMap.get(root);
		FundamentalEdge lastFunE = null;
		FundamentalEdge firstFunE = null;
		index = 0;
		Map<CoEdge, FundamentalEdge> funEdgeMap = new HashMap<CoEdge, FundamentalEdge>();
		// circle around the polygon
		int sourceEdgeCount = 1;
		double[] s1 = new double[3], s2 = new double[3], t1 = new double[3], t2 = new double[3];
		BigDecimal[] s1b = new BigDecimal[3], s2b = new BigDecimal[3], t1b = new BigDecimal[3], t2b = new BigDecimal[3];
		while (!visited.contains(eActive)) {
			CoVertex vTarget = eActive.getTargetVertex();
			if (branchSet.contains(vTarget)) {
				FundamentalVertex start = lastFunV;
				FundamentalVertex end = branchMap.get(vTarget);
				// hyperbolic motion
				CoEdge coEdge = cutInfo.edgeCutMap.get(eActive.getOppositeEdge());
				double[] lastTargetPoint = firstOfSegment.getTargetVertex().T;
				double[] lastStartPoint = cutInfo.edgeCutMap.get(firstOfSegment).getStartVertex().T;
				double[] actTargetPoint = vTarget.T;
				double[] actStartPoint = coEdge.getTargetVertex().T;
				
				// double precision isometry
				double[] T = P2.makeDirectIsometryFromFrames(null, 
					P2.projectP3ToP2(s1, lastStartPoint), 
					P2.projectP3ToP2(s2, actStartPoint), 
					P2.projectP3ToP2(t1, lastTargetPoint), 
					P2.projectP3ToP2(t2, actTargetPoint), 
					HYPERBOLIC
				);
				Matrix A = new Matrix(P2.imbedMatrixP2InP3(null, T));
	
				// arbitrary precision isometry
				BigDecimal[] TBig = P2Big.makeDirectIsometryFromFrames(null, 
					P2Big.projectP3ToP2(s1b, RnBig.toBig(null, lastStartPoint)), 
					P2Big.projectP3ToP2(s2b, RnBig.toBig(null, actStartPoint)), 
					P2Big.projectP3ToP2(t1b, RnBig.toBig(null, lastTargetPoint)), 
					P2Big.projectP3ToP2(t2b, RnBig.toBig(null, actTargetPoint)), 
					HYPERBOLIC,
					FundamentalPolygonUtility.context
				);
				BigDecimal[] ABig = P2Big.imbedMatrixP2InP3(null, TBig);
				BigDecimal[] startPos = RnBig.toBig(null, lastTargetPoint);
				FundamentalEdge fEdge = new FundamentalEdge(
					index++, 
					sourceEdgeCount, 
					start, 
					end, 
					startPos,
					A, 
					ABig
				);
				funEdgeMap.put(firstOfSegment, fEdge);
				edgeList.add(fEdge);
				// identification
				FundamentalEdge partner = funEdgeMap.get(coEdge);
				if (partner != null) {
					fEdge.partner = partner;
					partner.partner = fEdge;
				}
				
				// linkage
				fEdge.prevEdge = lastFunE;
				if (lastFunE != null) {
					lastFunE.nextEdge = fEdge;
				}
				lastFunE = fEdge;
				lastFunV = end;
				if (firstFunE == null) {
					firstFunE = fEdge;
				}
				firstOfSegment = eActive.getNextEdge().getOppositeEdge();
				sourceEdgeCount = 0;
			}
			
			sourceEdgeCount++;
			visited.add(eActive);
			eActive = eActive.getNextEdge();
		}
		assert lastFunE != null;
		assert firstFunE != null;
		lastFunE.nextEdge = firstFunE;
		firstFunE.prevEdge = lastFunE;
		return new FundamentalPolygon(edgeList);
	}
	
}
