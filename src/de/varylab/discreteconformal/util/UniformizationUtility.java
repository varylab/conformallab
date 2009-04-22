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
import geom3d.Point;

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
import de.jreality.math.Pn;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedge.plugin.AnnotationAdapter.EdgeAnnotation;
import de.jtem.halfedge.plugin.AnnotationAdapter.VertexAnnotation;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;

public class UniformizationUtility {

	private static int log = 3;
//	private static AnnotationAdapter<?>[] adapters = {new AngleSumAdapter(), new AlphaAdapter(), new FaceIndexAnnotation<CoFace>()};
	
	private static void log(int log, String msg) {
		if (UniformizationUtility.log >= log) {
			System.out.println(msg); 
		}
	}
	
	
	
	public static Matrix makeHyperbolicMotion(CoEdge s, CoEdge t) {
		Point s1 = s.getStartVertex().getTextureCoord();
		Point s2 = s.getTargetVertex().getTextureCoord();
		Point t1 = t.getStartVertex().getTextureCoord();
		Point t2 = t.getTargetVertex().getTextureCoord();
		return makeHyperbolicMotion(s1, s2, t1, t2);
	}
		
	public static Matrix makeHyperbolicMotion(Point s1, Point s2, Point t1, Point t2) {
		geom3d.Vector ns = new geom3d.Vector(s1).cross(s2);
		geom3d.Vector nt = new geom3d.Vector(t1).cross(t2);
		
		double s1s1 = Pn.innerProduct(s1.get(), s1.get(), HYPERBOLIC);
		double t1t1 = Pn.innerProduct(t1.get(), t1.get(), HYPERBOLIC);
		double s1s2 = Pn.innerProduct(s1.get(), s2.get(), HYPERBOLIC);
		double t1t2 = Pn.innerProduct(t1.get(), t2.get(), HYPERBOLIC);
		geom3d.Vector ws = new geom3d.Vector(s2).add(new Point(s1).times(s1s2 / s1s1).times(-1));
		geom3d.Vector wt = new geom3d.Vector(t2).add(new Point(t1).times(t1t2 / t1t1).times(-1));
		
		double nss1 = Pn.innerProduct(ns.get(), s1.get(), HYPERBOLIC);
		double ntt1 = Pn.innerProduct(nt.get(), t1.get(), HYPERBOLIC);
		ns.add(new Point(s1).times(nss1 / s1s1).times(-1));
		nt.add(new Point(t1).times(ntt1 / t1t1).times(-1));
		
		double nsws = Pn.innerProduct(ns.get(), ws.get(), HYPERBOLIC);
		double ntwt = Pn.innerProduct(nt.get(), wt.get(), HYPERBOLIC);	
		double wsws = Pn.innerProduct(ws.get(), ws.get(), HYPERBOLIC);
		double wtwt = Pn.innerProduct(wt.get(), wt.get(), HYPERBOLIC);
		ns.add(new Point(ws).times(nsws / wsws).times(-1));
		nt.add(new Point(wt).times(ntwt / wtwt).times(-1));
		Pn.normalize(ns.get(), ns.get(), HYPERBOLIC);
		Pn.normalize(nt.get(), nt.get(), HYPERBOLIC);
		
		double[] sa1 = s1.get();
		double[] sa2 = s2.get();
		double[] nsa = ns.get();
		double[] ta1 = t1.get();
		double[] ta2 = t2.get();
		double[] nta = nt.get();
		Matrix S = new Matrix(
			sa1[0], sa2[0], nsa[0], 0,
			sa1[1], sa2[1], nsa[1], 0,
			0,		0, 		0, 		1,
			sa1[2], sa2[2], nsa[2], 0
		);
		Matrix T = new Matrix(
			ta1[0], ta2[0], nta[0], 0,
			ta1[1], ta2[1], nta[1], 0,
			0,		0, 		0, 		1,
			ta1[2], ta2[2], nta[2], 0
		);
		return Matrix.times(T, S.getInverse());
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
			motion = new Matrix();
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
			Matrix motion
		) {
			this(index);
			this.start = start;
			this.end = end;
			this.motion = motion;
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
		
		
		public List<double[]> getOrbit(double[] root) {
			double[] pos = root.clone();
			List<double[]> result = new LinkedList<double[]>();
			Matrix A = new Matrix();
			for (FundamentalEdge e : edgeList) {
				Matrix T = new Matrix();
				FundamentalEdge active = e;
				do {
					T.multiplyOnLeft(active.partner.motion);
					active = active.partner.nextEdge;
				} while (active.partner != e);
				T.multiplyOnLeft(active.partner.motion);
				result.add(pos.clone());
				T.transformVector(pos);
				A.multiplyOnLeft(T);
			}
			System.out.println("Polygon Cycle check: \n" + A);
			return result;
		}
		
		
		
		public List<double[]> getDualOrbit(double[] root) {
			List<double[]> result = new LinkedList<double[]>();
			FundamentalEdge start = edgeList.get(0);
			FundamentalEdge active = start;
			Matrix T = new Matrix();
			do {
				double[] pos = T.multiplyVector(root);
				result.add(pos);
				T.multiplyOnRight(active.motion);
				active = active.partner.nextEdge;
			} while (active != start);
			return result;
		}
		
		
		@Override
		public String toString() {
			StringBuffer sb = new StringBuffer();
			sb.append("Fundamental Polygon edges: " + edgeList.size());
			sb.append("\n");
			for (FundamentalEdge fe : edgeList) {
				sb.append(fe + ": ");
				//sb.append(fe.start.index + " <-> " + fe.end.index);
				sb.append(" -> " + fe.partner);
				sb.append("\n");
			}
			sb.append("---------------------------------");
			return sb.toString();
		}
		
	}
	
	
	public static FundamentalPolygon constructFundamentalPolygon(
		CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo
	) {
		FundamentalPolygon poly = new FundamentalPolygon();
		CoVertex root = cutInfo.cutRoot;
		CoEdge rootEdge = null;
		for (CoEdge e : incomingEdges(root)) {
			if (e.getOppositeEdge().getLeftFace() == null) {
				rootEdge = e.getOppositeEdge();
				break;
			}
		}
		assert rootEdge != null;
		
		
		// create fundamental vertices
		Set<CoVertex> branchSet = cutInfo.getBranchSet();
		Map<CoVertex, FundamentalVertex> branchMap = new HashMap<CoVertex, FundamentalVertex>();
		
		Set<CoVertex> tmpSet = new HashSet<CoVertex>(branchSet);
		int index = 0;
		System.out.println("brachSet: " + branchSet);
		for (CoVertex v : branchSet) {
			if (!tmpSet.contains(v)) {
				continue;
			}
			Set<CoVertex> copies = cutInfo.getCopies(v);
			System.out.println("copies for " + index + ": " + copies);
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
		while (!visited.contains(eActive)) {
			
			CoVertex vTarget = eActive.getTargetVertex();
			if (branchSet.contains(vTarget)) {
				FundamentalVertex start = lastFunV;
				FundamentalVertex end = branchMap.get(vTarget);
				// hyperbolic motion
				CoEdge coEdge = cutInfo.edgeCutMap.get(eActive.getOppositeEdge());
				Matrix A = makeHyperbolicMotion(coEdge, eActive);
				
				FundamentalEdge fEdge = new FundamentalEdge(index++, start, end, A);
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
		
		System.out.println("Cutted Polygon:\n" + poly);
		
		// canonical polygon construction --------------
		FundamentalVertex fRoot = poly.edgeList.get(0).start;
		List<FundamentalEdge> newPoly = new LinkedList<FundamentalEdge>();
		Set<FundamentalEdge> deleted = new TreeSet<FundamentalEdge>();
		for (FundamentalEdge e : poly.edgeList) {
			if (deleted.contains(e)) {
				continue;
			}
			FundamentalVertex badVertex = e.end;
			if (badVertex != fRoot) {
				deleted.add(e);
				deleted.add(e.partner);
				for (FundamentalEdge fe : poly.edgeList) {
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
		poly.edgeList = newPoly;
		
		// linkage
		for (int i = 0; i < newPoly.size(); i++) {
			FundamentalEdge prev = newPoly.get(i - 1 < 0 ? newPoly.size() - 1 : i - 1);
			FundamentalEdge next = newPoly.get((i + 1) % newPoly.size());
			FundamentalEdge act = newPoly.get(i);
			prev.nextEdge = act;
			act.prevEdge = prev;
			act.nextEdge = next;
			next.prevEdge = act;
		}
		System.out.println("Canonical Polygon:\n" + poly);
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
	
	
	
	
	public static class AngleSumAdapter extends VertexAnnotation<CoVertex> {

//		private DecimalFormat
//			df = new DecimalFormat("#.####");
		
		@Override
		public String getText(CoVertex n) {
			return "" + n.getIndex();// + " Σ" + df.format(getAngleSum(n));
		}
		
	}
	
	
	public static class AlphaAdapter extends EdgeAnnotation<CoEdge> {

		private DecimalFormat
			df = new DecimalFormat("#.####");
			
		
		@Override
		public String getText(CoEdge n) {
			return n.getIndex() + " α:" + df.format(n.getAlpha());
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
