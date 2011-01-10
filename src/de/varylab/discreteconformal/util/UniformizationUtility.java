package de.varylab.discreteconformal.util;

import static de.jreality.math.Pn.HYPERBOLIC;
import static de.jtem.halfedge.util.HalfEdgeUtils.incomingEdges;
import static java.lang.Double.isNaN;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.acos;
import static java.lang.Math.cos;
import static java.lang.Math.cosh;
import static java.lang.Math.exp;
import static java.lang.Math.sin;
import static java.lang.Math.sinh;
import static java.lang.Math.sqrt;
import static java.util.Collections.sort;

import java.math.BigDecimal;
import java.math.MathContext;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;
import java.util.TreeSet;

import no.uib.cipr.matrix.Vector;
import de.jreality.math.Matrix;
import de.jreality.math.P2;
import de.jreality.math.Pn;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Node;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.impl.StringAdapter;
import de.jtem.halfedgetools.adapter.type.Label;
import de.jtem.halfedgetools.util.TriangulationException;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;

public class UniformizationUtility {

	private static int 
		log = 3;
	public static MathContext 
		context = new MathContext(40);
	
	public static void normalize(Matrix M) {
		double max = 0;
		for (double d : M.getArray()) {
			if (Math.abs(d) > max) {
				max = d;
			}
		}
		M.times(1 / max);
	}
	
	
	private static void log(int log, String msg) {
		if (UniformizationUtility.log >= log) {
			System.out.println(msg); 
		}
	}
	
	
	public static class FundamentalVertex {
		
		public int 
			index = 0;
		
		public FundamentalVertex(int index) {
			this.index = index;
		}

		@Override
		public String toString() {
			return "FunVertex " + index;
		}
		
	}
	
	
	public static class FundamentalEdge implements Comparable<FundamentalEdge> {
		
		public int 
			index = 0;
		public FundamentalVertex 
			start = new FundamentalVertex(0),
			end = new FundamentalVertex(0);
		public Matrix
			motion = new Matrix(); // identification motion
		public BigDecimal[]
		    motionBig = null;
		public FundamentalEdge
			prevEdge = null,
			nextEdge = null;
		public FundamentalEdge
			partner = null;

		public FundamentalEdge(int index) {
			this.index = index;
		}

		public FundamentalEdge(
			int index,
			FundamentalVertex start, 
			FundamentalVertex end,
			Matrix motion,
			BigDecimal[] motionBig
		) {
			this(index);
			this.start = start;
			this.end = end;
			this.motion = motion;
			this.motionBig = motionBig;
		}
		
		@Override
		public String toString() {
			return "FunEdge " + index;
		}

		@Override
		public int compareTo(FundamentalEdge o) {
			return index - o.index;
		}

	}
	
	
	public static class FundamentalPolygon {
		
		public List<FundamentalEdge>
			edgeList = new LinkedList<FundamentalEdge>();
		
		public int getLength() {
			return edgeList.size();
		}
		
		
		public int getNumVertices() {
			Set<FundamentalVertex> vSet = new HashSet<UniformizationUtility.FundamentalVertex>();
			for (FundamentalEdge e : edgeList) {
				vSet.add(e.end);
				vSet.add(e.start);
			}
			return vSet.size();
		}
		
		
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
				RnBig.matrixTimesVector(pos1, active1.partner.motionBig, pos2, context);
				RnBig.matrixTimesVector(pos2, active2.partner.motionBig, pos2, context);
				PnBig.normalize(pos1, pos1, HYPERBOLIC, context);
				PnBig.normalize(pos2, pos2, HYPERBOLIC, context);
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
				RnBig.times(Tbig, Tbig, active.motionBig, context);
//				P3.orthonormalizeMatrix(T.getArray(), T.getArray(), 1E-5, HYPERBOLIC);
//				if (P3.isometryIsUnstable(T.getArray(), HYPERBOLIC)) {
//					System.out.println("Isometry is unstable \n" + T);
//				}
//				normalize(T);
				active = active.partner.nextEdge;
			} while (active != start);
			System.out.println("\nDual orbit transform: \n" + T);
			System.out.println("\nBig dual orbit transform: \n" + RnBig.matrixToString(Tbig));
			return result;
		}
		
		
		public FundamentalPolygon getMinimal() {
			// canonical polygon construction --------------
			FundamentalVertex fRoot = edgeList.get(0).start;
			List<FundamentalEdge> newPoly = new LinkedList<FundamentalEdge>();
			Set<FundamentalEdge> deleted = new TreeSet<FundamentalEdge>();
			for (FundamentalEdge e : edgeList) {
				if (deleted.contains(e)) {
					continue;
				}
				FundamentalVertex badVertex = e.end;
				if (badVertex != fRoot) {
					deleted.add(e);
					deleted.add(e.partner);
					for (FundamentalEdge fe : edgeList) {
						if (fe.start == badVertex) {
							fe.start = fRoot;
						}
						if (fe.end == badVertex) {
							fe.end = fRoot;
						}
					}
				} else {
					newPoly.add(e);
				}
			}
			edgeList = newPoly;
			// linkage
			for (int i = 0; i < newPoly.size(); i++) {
				FundamentalEdge prev = newPoly.get(i - 1 < 0 ? newPoly.size() - 1 : i - 1);
				FundamentalEdge next = newPoly.get((i + 1) % newPoly.size());
				FundamentalEdge act = newPoly.get(i);
				prev.nextEdge = act;
				act.index = i;
				act.prevEdge = prev;
				act.nextEdge = next;
				next.prevEdge = act;
			}
			FundamentalPolygon result = new FundamentalPolygon();
			result.edgeList = newPoly;
			return result;
		}
		
		
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
			List<FundamentalEdge> canList = new LinkedList<UniformizationUtility.FundamentalEdge>();
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
//			do {
//				System.out.println(act.motion);
//				act = act.nextEdge;
//			} while (act != e);
		}
		
		private FundamentalEdge findLinkedEdge(FundamentalEdge a) {
			Set<FundamentalEdge> checkSet = new TreeSet<UniformizationUtility.FundamentalEdge>();
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
			Set<FundamentalEdge> cSet = new TreeSet<UniformizationUtility.FundamentalEdge>();
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
				RnBig.times(c.motionBig, ABiginv, c.motionBig, context);
				System.out.println(c.partner.index + " = " + a.index + " " + c.partner.index);
				c.partner.motion.multiplyOnRight(A);
				RnBig.times(c.partner.motionBig, c.partner.motionBig, ABig, context);
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
	
	
	public static FundamentalPolygon constructFundamentalPolygon(
		CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo
	) {
 		FundamentalPolygon poly = new FundamentalPolygon();
		
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
//		System.out.println("brachSet: " + branchSet);
		for (CoVertex v : branchSet) {
			if (!tmpSet.contains(v)) {
				continue;
			}
			Set<CoVertex> copies = cutInfo.getCopies(v);
//			System.out.println("copies for " + index + ": " + copies);
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
				System.out.println("Motion:\n" + A);

				// arbitrary precision isometry
				BigDecimal[] TBig = P2Big.makeDirectIsometryFromFrames(null, 
					P2Big.projectP3ToP2(s1b, RnBig.toBig(null, lastStartPoint)), 
					P2Big.projectP3ToP2(s2b, RnBig.toBig(null, actStartPoint)), 
					P2Big.projectP3ToP2(t1b, RnBig.toBig(null, lastTargetPoint)), 
					P2Big.projectP3ToP2(t2b, RnBig.toBig(null, actTargetPoint)), 
					HYPERBOLIC,
					context
				);
				BigDecimal[] ABig = P2Big.imbedMatrixP2InP3(null, TBig);
				System.out.println("Big Motion:\n" + RnBig.matrixToString(ABig));
				
				
				FundamentalEdge fEdge = new FundamentalEdge(index++, start, end, A, ABig);
				funEdgeMap.put(firstOfSegment, fEdge);
				poly.edgeList.add(fEdge);
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
			}
			
			visited.add(eActive);
			eActive = eActive.getNextEdge();
		}
		assert lastFunE != null;
		assert firstFunE != null;
		lastFunE.nextEdge = firstFunE;
		firstFunE.prevEdge = firstFunE;
		
		// try dual orbit of the raw polygon
		System.out.println("Cutted Polygon:\n" + poly);
		poly.getDualOrbit(new double[] {0,0,0,1});
		
		// reduce to minimal polygon
		poly = poly.getMinimal();
		System.out.println("Minimal Polygon:\n" + poly);
		poly.getDualOrbit(new double[] {0,0,0,1});
		
//		// construct the canonical polygon
//		poly = poly.getCanonical();
//		System.out.println("Canonical Polygon:\n" + poly);
//		poly.getDualOrbit(new double[] {0,0,0,1});
		return poly;
	}
	
	
	
	
	
	
	
	public static class UAdapter {
		
		private Vector
			u = null;
		
		public UAdapter(Vector u) {
			this.u = u;
		}
		
		public double getU(CoVertex v) {
			if (v.getSolverIndex() < 0 || v.getSolverIndex() > u.size() - 1) {
				return 0.0;
			} else {
				return u.get(v.getSolverIndex());
			}
		}
		
		public void setU(Vector u) {
			this.u = u;
		}
		
	}
	
	
	public static Map<CoEdge, Double> getLengthMap(		
		CoHDS hds, 
		UAdapter u
	) {
		Map<CoEdge, Double> lMap = new HashMap<CoEdge, Double>();
		for (CoEdge e : hds.getPositiveEdges()) {
			double l = getNewLength(e, u);
			lMap.put(e, l);
			lMap.put(e.getOppositeEdge(), l);
		}
		return lMap;
	}
	

	public static void toFundamentalPolygon(
		CoHDS hds, 
		CoVertex root,
		UAdapter u
	) {
		Map<CoEdge, Double> lMap = getLengthMap(hds, u);
		checkLengthsAndAngles(hds, lMap);
		int tries = 0;
		while (hds.numVertices() > 1 && tries++ < 10) {
			int flipSuccess = 0;
			int flipFailed = 0;
			
			List<CoVertex> vList = new ArrayList<CoVertex>(hds.getVertices());
			Collections.sort(vList, new ValenceComparator());
//			Collections.shuffle(vList);
			
			for (CoVertex v : vList) {
				if (hds.numVertices() == 1) {
					break;
				}
				Set<CoEdge> starUnique = new HashSet<CoEdge>(incomingEdges(v));
				List<CoEdge> star = new LinkedList<CoEdge>(starUnique);
				sort(star, new EdgeAngleDivideComparator());
				
				int degree = star.size();
				int flipTries = degree * 3; 
				for (int flips = 0; flips < flipTries; flips++) {
					CoEdge flipEdge = star.get(flips % degree);
					if (!isFlippable(flipEdge)) {
						continue;
					}
					boolean flipped = flipLengthAndAlphas(flipEdge, lMap);
					if (flipped) {
						flipSuccess++;
					} else {
						flipFailed++;
					}
					starUnique = new HashSet<CoEdge>(incomingEdges(v));
					star = new LinkedList<CoEdge>(starUnique);
					sort(star, new EdgeAngleDivideComparator());
					degree = star.size();
					if (degree == 3) {
						break;
					}
				}
				if (degree == 3) {
					removeTrivalentVertex(v);
				} else {
					System.out.println("could not remove vertex " + v + " degree " + degree);
				}
			}
			System.out.println("try " + tries + ": vertices " + hds.numVertices() + ", successfull flips " + flipSuccess + ", flips failed " + flipFailed);
		}
		log(3, "tries: " + tries);
		for (CoVertex v : hds.getVertices()) {
			log(3, "vertex angle sum " + v + ": " + (2*PI - getAngleSum(v)));
		}
		for (CoFace f : hds.getFaces()) {
			log(3, "face angle sum " + f + ": " + (PI - getAngleSum(f)));
		}
		checkLengthsAndAngles(hds, lMap);
		log(3, hds.toString());
	}
	
	
	private static boolean isFlippable(CoEdge e) {
		double alpha = 0.0;
		alpha += e.getNextEdge().getAlpha();
		alpha += e.getOppositeEdge().getPreviousEdge().getAlpha();
		double beta = 0.0;
		beta += e.getPreviousEdge().getAlpha();
		beta += e.getOppositeEdge().getNextEdge().getAlpha();
		log(4, "check flippable: alpha:" + (alpha - PI) + " beta:" + (beta - PI));
		return alpha < PI && beta < PI;
	}
	
	
	public static boolean flipLengthAndAlphas(
		CoEdge e, 
		Map<CoEdge, Double> lMap
	) {  
		CoEdge cEdge = e;
		CoEdge cPEdge = e.getOppositeEdge();
		CoEdge aEdge = cEdge.getNextEdge();
		CoEdge aPEdge = cPEdge.getPreviousEdge();
		CoEdge bEdge = cEdge.getPreviousEdge();
		CoEdge bPEdge = cPEdge.getNextEdge();

		double a = lMap.get(aEdge);
		double aP = lMap.get(aPEdge);
		double b = lMap.get(bEdge);
		double bP = lMap.get(bPEdge);
		double c = lMap.get(e);
//		double alpha = aEdge.getAlpha();
//		double alphaP = aPEdge.getAlpha();
//		double beta = bEdge.getAlpha();
//		double betaP = bPEdge.getAlpha();
//		double gamma = cEdge.getAlpha();
//		double gammaP = cPEdge.getAlpha();
		double alpha = calculateAlpha(a, b, c);
		double beta = calculateAlpha(b, c, a);
		double alphaP = calculateAlpha(aP, c, bP);
		double betaP = calculateAlpha(bP, c, aP);
//		double gamma = calculateAlpha(c, a, b);
//		double gammaP = calculateAlpha(a, aP, bP);
		double betaSum = beta + betaP;
		double alphaSum = alpha + alphaP;
		double coshC = cosh(a)*cosh(aP) - sinh(a)*sinh(aP)*cos(betaSum);
		double cNew = arcosh(coshC);
//		double cPCheck = arcosh(cosh(b)*cosh(bP) - sinh(b)*sinh(bP)*cos(alphaSum));
//		
//		// check compatibility of both ways of calculation
//		if (abs(cP - cPCheck) > 1E-5) {
//			log(4, "lengths differ too much");
//			return false;
//		}  
//		
//		// take the mean
//		double l = (cP + cPCheck) / 2;
//		if (Double.isNaN(l) || Double.isInfinite(l)) {
//			log(4, "We have a NaN/Inf length here -> will not flip");
//			return false;
//		}
		
		// new angles
		double alphaNew = calculateAlpha(a, aP, cNew);
		double betaNew = calculateAlpha(b, bP, cNew);
		double alphaPNew = calculateAlpha(aP, a, cNew);
		double betaPNew = calculateAlpha(bP, b, cNew);
		double gammaA = calculateAlpha(cNew, a, aP);
		double gammaB = calculateAlpha(cNew, b, bP);
//		if (Double.isNaN(alphaNew) || Double.isNaN(alphaPNew)
//				|| Double.isNaN(betaNew) || Double.isNaN(betaPNew)) {
//			log(4, "We have a new NaN/0.0 angle -> no flip possible");
//			return false;
//		}
//
//		// check angle calculations by comparison with the known sums
//		double checkGammaP = alphaNew + betaNew;
//		double checkGamma = alphaPNew + betaPNew;
//		if (abs(checkGamma - gamma) > 1E-5 || abs(checkGammaP - gammaP) > 1E-5) {
//			log(4, "Error calculating new angles -> no flip");
//			return false;
//		}
//		
//		// check if the triangles are hyperbolic 
//		double checkTop = betaSum + alphaNew + alphaPNew;
//		double checkBottom = alphaSum + betaNew + betaPNew;
//		if (checkTop >= PI || checkBottom >= PI) {
//			log(4, "no hyperbolic triangle after flip -> no flip");
//			return false;
//		}

		boolean check1 = checkTriangle(a, aP, cNew, alphaNew, alphaPNew, betaSum);
		boolean check2 = checkTriangle(b, bP, cNew, betaNew, betaPNew, alphaSum);
		if (!(check1 && check2)) {
			return false; 
		}

		// flip the edge
		flip(e);
		
		// insert new angles and lengths
		aEdge.setAlpha(alphaNew);
		bEdge.setAlpha(betaNew);
		aPEdge.setAlpha(alphaPNew);
		bPEdge.setAlpha(betaPNew);
		cEdge.setAlpha(gammaA);
		cPEdge.setAlpha(gammaB); 
		lMap.put(e, cNew);
		lMap.put(e.getOppositeEdge(), cNew);

		return true;
	}
	
	
	public static boolean checkLengthsAndAngles(CoHDS hds, Map<CoEdge, Double> lMap) {
		for (CoFace f : hds.getFaces()) {
			if (!checkTriangle(f, lMap)) {
				System.out.println("checkLengthsAndAngles failed for face " + f);
				return false;
			}
		}
		System.out.println("checkLengthsAndAngles succeeded");
		return true;
	}
	
	
	public static boolean checkTriangle(CoFace f, Map<CoEdge, Double> lMap) {
		CoEdge aEdge = f.getBoundaryEdge();
		CoEdge bEdge = aEdge.getNextEdge();
		CoEdge cEdge = bEdge.getNextEdge();
		double alpha = aEdge.getAlpha();
		double beta = bEdge.getAlpha();
		double gamma = cEdge.getAlpha();
		double a = lMap.get(aEdge);
		double b = lMap.get(bEdge);
		double c = lMap.get(cEdge);
		return checkTriangle(a, b, c, alpha, beta, gamma);
	}
	
	
	public static boolean checkTriangle(
		double a, 
		double b, 
		double c, 
		double alpha,
		double beta,
		double gamma
	) {
		double eps = 1E-10;
		if (isNaN(a) || isNaN(b) || isNaN(c) || isNaN(alpha) || isNaN(beta) || isNaN(gamma)) {
			return false;
		}
		// check law of cosines
		double checkA = cosh(a) - cosh(b)*cosh(c) + sinh(b)*sinh(c)*cos(alpha);
		double checkB = cosh(b) - cosh(a)*cosh(c) + sinh(a)*sinh(c)*cos(beta);
		double checkC = cosh(c) - cosh(a)*cosh(b) + sinh(a)*sinh(b)*cos(gamma);
		if (abs(checkA) > eps || abs(checkB) > eps || abs(checkC) > eps) {
			return false;
		}
		// check law of cosines II
		checkA = cos(alpha) + cos(beta)*cos(gamma) - sin(beta)*sin(gamma)*cosh(a);
		checkB = cos(beta) + cos(alpha)*cos(gamma) - sin(alpha)*sin(gamma)*cosh(b);
		checkC = cos(gamma) + cos(alpha)*cos(beta) - sin(alpha)*sin(beta)*cosh(c);
		if (abs(checkA) > eps || abs(checkB) > eps || abs(checkC) > eps) {
			return false;
		}
		return true;
	}
	
	
	
	
	private static double getAngleSum(CoVertex v) {
		double sum = 0.0;
		for (CoEdge e : incomingEdges(v)) {
			sum += e.getPreviousEdge().getAlpha();
		}
		return sum;
	}
	
	
	private static double getAngleSum(CoFace f) {
		double sum = 0.0;
		CoEdge e = f.getBoundaryEdge();
		sum += e.getAlpha();
		e = e.getNextEdge();
		sum += e.getAlpha();
		e = e.getNextEdge();
		sum += e.getAlpha();
		return sum;
	}
	
	
	
	private static double calculateAlpha(double a, double b, double c) {
		return acos((cosh(b)*cosh(c) - cosh(a)) / (sinh(b)*sinh(c)));
	}
	
	
	private static double getNewLength(
		CoEdge e, 
		UAdapter u
	) {
		CoVertex v1 = e.getStartVertex();
		CoVertex v2 = e.getTargetVertex();
		Double u1 = v1.getSolverIndex() >= 0 ? u.getU(v1) : 0.0; 
		Double u2 = v2.getSolverIndex() >= 0 ? u.getU(v2) : 0.0;
		Double la = e.getLambda();
		Double lambdaNew = la + u1 + u2;
		return 2 * arsinh( exp(lambdaNew / 2) );
	}
	
	
	
	private static void removeTrivalentVertex(CoVertex v) {
		log(4, "remove vertex " + v);
		HalfEdgeDataStructure<CoVertex, CoEdge, CoFace> hds = v.getHalfEdgeDataStructure();
		CoEdge in1 = v.getIncomingEdge();
		CoEdge out1 = in1.getOppositeEdge();
		CoEdge outer1 = in1.getPreviousEdge();
		CoEdge in2 = out1.getPreviousEdge();
		CoEdge out2 = in2.getOppositeEdge();
		CoEdge outer2 = in2.getPreviousEdge();
		CoEdge in3 = out2.getPreviousEdge();
		CoEdge out3 = in3.getOppositeEdge();
		CoEdge outer3 = in3.getPreviousEdge();
		
		outer1.setAlpha(out1.getAlpha() + in3.getAlpha());
		outer2.setAlpha(out2.getAlpha() + in1.getAlpha());
		outer3.setAlpha(out3.getAlpha() + in2.getAlpha());
		
		CoFace f1 = in1.getLeftFace();
		CoFace f2 = in2.getLeftFace();
		CoFace f3 = in3.getLeftFace();

		hds.removeVertex(v);
		hds.removeEdge(in1);
		hds.removeEdge(in2);
		hds.removeEdge(in3);
		hds.removeEdge(out1);
		hds.removeEdge(out2);
		hds.removeEdge(out3);
		hds.removeFace(f2);
		hds.removeFace(f3);
		
		outer1.setLeftFace(f1);
		outer2.setLeftFace(f1);
		outer3.setLeftFace(f1);
		
		outer1.linkNextEdge(outer2);
		outer2.linkNextEdge(outer3);
		outer3.linkNextEdge(outer1);
	}
	
	
	
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> void flip(E e) {
		F leftFace = e.getLeftFace();
		F rightFace = e.getRightFace();
		if (leftFace == rightFace)
			return;
		E a1 = e.getOppositeEdge().getNextEdge();
		E a2 = a1.getNextEdge();
		E b1 = e.getNextEdge();
		E b2 = b1.getNextEdge();
		
		V v1 = e.getStartVertex();
		V v2 = e.getTargetVertex();
		V v3 = a1.getTargetVertex();
		V v4 = b1.getTargetVertex();

		//new connections
		e.linkNextEdge(a2);
		e.linkPreviousEdge(b1);
		e.getOppositeEdge().linkNextEdge(b2);
		e.getOppositeEdge().linkPreviousEdge(a1);
		e.setTargetVertex(v3);
		e.getOppositeEdge().setTargetVertex(v4);
		
		a2.linkNextEdge(b1);
		b2.linkNextEdge(a1);
		
		//set faces
		b2.setLeftFace(rightFace);
		a2.setLeftFace(leftFace);
		b2.setTargetVertex(v1);
		a2.setTargetVertex(v2);
		a1.setTargetVertex(v3);
		b1.setTargetVertex(v4);
	}
	
	
	
	private static double arsinh(double x) {
		double r = x + sqrt(x*x + 1);
		return Math.log(r);
	}
	
	private static double arcosh(double x) {
		double r = x + sqrt(x*x - 1);
		return Math.log(r);
	}
	
	
	/**
	 * Checks whether this edge is locally delaunay
	 * @param edge the edge to check
	 * @return the check result
	 */
	public static  boolean isDelaunay(CoEdge edge) {
		Double gamma = edge.getNextEdge().getAlpha();
		Double delta = edge.getOppositeEdge().getNextEdge().getAlpha();
		return gamma + delta <= Math.PI - 1E-1;
	}
	
	
	/**
	 * Returns the positive edges belonging to the kite of the edge 
	 */
	public static List<CoEdge> getPositiveKiteBorder(CoEdge edge){
		LinkedList<CoEdge> result = new LinkedList<CoEdge>();
		CoEdge e1 = edge.getNextEdge();
		CoEdge e2 = edge.getPreviousEdge();
		CoEdge e3 = edge.getOppositeEdge().getNextEdge();
		CoEdge e4 = edge.getOppositeEdge().getPreviousEdge();
		if (!e1.isPositive())
			e1 = e1.getOppositeEdge();
		if (!e2.isPositive())
			e2 = e2.getOppositeEdge();
		if (!e3.isPositive())
			e3 = e3.getOppositeEdge();
		if (!e4.isPositive())
			e4 = e4.getOppositeEdge();
		result.add(e1);
		result.add(e2);
		result.add(e3);
		result.add(e4);
		return result;
	}
	
	
	/**
	 * Constructs the delaunay triangulation of the given structure.
	 * @param graph must be a triangulation
	 * @throws TriangulationException if the given graph is no triangulation or 
	 * if the triangle inequality doesn't hold for some triangle
	 * @return the number of flips performed
	 */
	public static int constructDelaunay(
		CoHDS graph,
		Map<CoEdge, Double> lMap
	) {
		int result = 0;
		Set<CoEdge> markSet = new HashSet<CoEdge>();
		Stack<CoEdge> stack = new Stack<CoEdge>();
		for (CoEdge positiveEdge : graph.getPositiveEdges()){
			markSet.add(positiveEdge);
			stack.push(positiveEdge);
		}
		while (!stack.isEmpty()){
			CoEdge ab = stack.pop();
			markSet.remove(ab);
			if (!isDelaunay(ab) && isFlippable(ab)){
				flipLengthAndAlphas(ab, lMap);
				for (CoEdge xy : getPositiveKiteBorder(ab)){
					if (!markSet.contains(xy)){
						markSet.add(xy);
						stack.push(xy);
					}
				}
				result++;
			}
		}
		return result;
	}	
	
	
	
	@Label
	public static class AngleSumAdapter extends StringAdapter {

		public AngleSumAdapter() {
			super(true, false);
		}
		
		@Override
		public <
			V extends Vertex<V, E, F>,
			E extends Edge<V, E, F>,
			F extends Face<V, E, F>,
			N extends Node<V, E, F>
		> String get(N n, AdapterSet a) {
			return "" + n.getIndex();
		}
		
		@Override
		public <N extends Node<?, ?, ?>> boolean canAccept(Class<N> nodeClass) {
			return true;
		}

		@Override
		public double getPriority() {
			return 1;
		}
		
	}
	
	@Label
	public static class AlphaAdapter extends StringAdapter {

		public AlphaAdapter() {
			super(true, false);
		}

		private DecimalFormat
			df = new DecimalFormat("#.####");
		
		@Override
		public <
			V extends Vertex<V, E, F>,
			E extends Edge<V, E, F>,
			F extends Face<V, E, F>,
			N extends Node<V, E, F>
		> String get(N n, AdapterSet a) {
			CoEdge ce = (CoEdge)n;
			return n.getIndex() + " α:" + df.format(ce.getAlpha());
		}


		@Override
		public <N extends Node<?, ?, ?>> boolean canAccept(Class<N> nodeClass) {
			return CoEdge.class.isAssignableFrom(nodeClass);
		}


		@Override
		public double getPriority() {
			return 1;
		}
		
	}
	

	private static class EdgeAngleDivideComparator implements Comparator<CoEdge> {

		@Override
		public int compare(CoEdge o1, CoEdge o2) {
			double left1 = o1.getPreviousEdge().getAlpha();
			double right1 = o1.getOppositeEdge().getNextEdge().getAlpha();
			double left2 = o2.getPreviousEdge().getAlpha();
			double right2 = o2.getOppositeEdge().getNextEdge().getAlpha();
			return (left1 + right1) - (left2 + right2) < 0 ? -1 : 1;
		}
		
	}
	
	
	private static class ValenceComparator implements Comparator<CoVertex> {
		
		@Override
		public int compare(CoVertex o1, CoVertex o2) {
			int v1 = HalfEdgeUtils.incomingEdges(o1).size();
			int v2 = HalfEdgeUtils.incomingEdges(o2).size();
			return v1 - v2;
		}
		
	}
	
}
