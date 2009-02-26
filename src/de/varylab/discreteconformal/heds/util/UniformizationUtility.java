package de.varylab.discreteconformal.heds.util;

import static de.jtem.halfedge.util.HalfEdgeUtils.getGenus;
import static de.jtem.halfedge.util.HalfEdgeUtils.incomingEdges;
import static de.varylab.discreteconformal.plugin.DiscreteConformalPlugin.halfedgeDebugger;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.acos;
import static java.lang.Math.cos;
import static java.lang.Math.cosh;
import static java.lang.Math.exp;
import static java.lang.Math.sinh;
import static java.lang.Math.sqrt;
import static java.util.Collections.sort;

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

import javax.swing.text.DefaultEditorKit.CutAction;

import no.uib.cipr.matrix.Vector;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedge.functional.circlepattern.hds.CPEdge;
import de.jtem.halfedge.plugin.AnnotationAdapter;
import de.jtem.halfedge.plugin.AnnotationAdapter.EdgeAnnotation;
import de.jtem.halfedge.plugin.AnnotationAdapter.FaceIndexAnnotation;
import de.jtem.halfedge.plugin.AnnotationAdapter.VertexAnnotation;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;

public class UniformizationUtility {

	private static int log = 3;
	private static AnnotationAdapter<?>[] adapters = {new AngleSumAdapter(), new AlphaAdapter(), new FaceIndexAnnotation<CoFace>()};
	
	private static void log(int log, String msg) {
		if (UniformizationUtility.log >= log) {
			System.out.println(msg); 
		}
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
	

	public static void reduceToFundamentalPolygon(
		CoHDS hds, 
		CoVertex root,
		UAdapter u
	) {  
		Map<CoEdge, Double> lMap = getLengthMap(hds, u);
		int tries = 0;
		while (hds.numVertices() > 1 && tries++ < 5) {
			int flipSuccess = 0;
			int flipFailed = 0;
			
			List<CoVertex> vList = new ArrayList<CoVertex>(hds.getVertices());
//			Collections.sort(vList, new ValenceComparator());
			Collections.shuffle(vList);
			
			for (CoVertex v : vList) {
				if (hds.numVertices() == 1) {
					break;
				}
				Set<CoEdge> starUnique = new HashSet<CoEdge>(incomingEdges(v));
				List<CoEdge> star = new LinkedList<CoEdge>(starUnique);
				sort(star, new EdgeAngleDivideComparator());
				
				int degree = star.size();
				for (int flips = 0; flips < degree * 10; flips++) {
					CoEdge flipEdge = star.get(flips % degree);
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
		log(3, hds.toString());
		if (hds.numVertices() == 1) {
			cutToDisk(hds);
		}
	}
	
	
	
	
	public static void cutToDisk(CoHDS hds) {
		int g = getGenus(hds);
		System.out.println("Genus before cutting: " + g);
//		Map<CoVertex, CoVertex> vMap = CuttingUtility.cutLoopEdge(hds.getEdge(0));
//		for (CoVertex v : vMap.keySet()) {
//			CoVertex newV = vMap.get(v);
//			newV.setPosition(v.getPosition());
//			newV.setSolverIndex(v.getSolverIndex());
//		}
//		
//		for (int i = 0; i < 2 * g; i++) {
//			CoEdge cutEdge = null;  
//			for (CoEdge e : hds.getPositiveEdges()) {
//				CoVertex s = e.getStartVertex();
//				CoVertex t = e.getTargetVertex();
//				if (s != t && HalfEdgeUtils.isInteriorEdge(e)) {
//					cutEdge = e;
//					break;
//				}
//			}
//			if (cutEdge == null) {
//				System.out.println("alarm!");
//			}
//			assert cutEdge != null;
//			vMap = CuttingUtility.cutAtEdge(cutEdge);
//			for (CoVertex v : vMap.keySet()) {
//				CoVertex newV = vMap.get(v);
//				newV.setPosition(v.getPosition());
//				newV.setSolverIndex(v.getSolverIndex());
//			}
//		}
		Set<CoEdge> edges = new HashSet<CoEdge>(hds.getEdges());
		Set<CoEdge> tree = SpanningTreeUtility.getDualSpanningTree(edges, hds.getEdge(0));
		edges.removeAll(tree);
		Set<CoEdge> cutCycles = new HashSet<CoEdge>();
		for (CoEdge e : edges) {
			if (!cutCycles.contains(e.getOppositeEdge())) {
				cutCycles.add(e);
			}
		}
		System.out.println("Cuts will be made at: " + cutCycles);
		for (CoEdge e : cutCycles) {
			Map<CoVertex, CoVertex> vMap = CuttingUtility.cutAtEdge(e);
			for (CoVertex v : vMap.keySet()) {
				CoVertex newV = vMap.get(v);
				newV.setPosition(v.getPosition());
				newV.setSolverIndex(v.getSolverIndex());
			}
		}
		System.out.println("Genus after cutting: " + getGenus(hds));
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
	
	
	
	public static class AngleSumAdapter extends VertexAnnotation<CoVertex> {

		private DecimalFormat
			df = new DecimalFormat("#.####");
		
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
		double alpha = aEdge.getAlpha();
		double alphaP = aPEdge.getAlpha();
		double beta = bEdge.getAlpha();
		double betaP = bPEdge.getAlpha();
		double gamma = cEdge.getAlpha();
		double gammaP = cPEdge.getAlpha();
		double betaSum = beta + betaP;
		double alphaSum = alpha + alphaP;
		double cP = arcosh(cosh(a)*cosh(aP) - sinh(a)*sinh(aP)*cos(betaSum));
		double cPCheck = arcosh(cosh(b)*cosh(bP) - sinh(b)*sinh(bP)*cos(alphaSum));
		
		// check compatibility of both ways of calculation
		if (abs(cP - cPCheck) > 1E-5) {
			log(4, "lengths differ too much");
			return false;
		}  
		
		// take the mean
		double l = (cP + cPCheck) / 2;
		if (Double.isNaN(l) || Double.isInfinite(l)) {
			log(4, "We have a NaN/Inf length here -> will not flip");
			return false;
		}
		
		// new angles
		double alphaNew = calculateAlpha(a, aP, l);
		double betaNew = calculateAlpha(b, bP, l);
		double alphaPNew = calculateAlpha(aP, a, l);
		double betaPNew = calculateAlpha(bP, b, l);
		if (Double.isNaN(alphaNew) || Double.isNaN(alphaPNew)
				|| Double.isNaN(betaNew) || Double.isNaN(betaPNew)) {
			log(4, "We have a new NaN/0.0 angle -> no flip possible");
			return false;
		}

		// check angle calculations by comparison with the known sums
		double checkGammaP = alphaNew + betaNew;
		double checkGamma = alphaPNew + betaPNew;
		if (abs(checkGamma - gamma) > 1E-5 || abs(checkGammaP - gammaP) > 1E-5) {
			log(4, "Error calculating new angles -> no flip");
			return false;
		}
		
		// check if the triangles are hyperbolic 
		double checkTop = betaSum + alphaNew + alphaPNew;
		double checkBottom = alphaSum + betaNew + betaPNew;
		if (checkTop >= PI || checkBottom >= PI) {
			log(4, "no hyperbolic triangle after flip -> no flip");
			return false;
		}


		// flip the edge
		flip(e);
		
		// insert new angles and lengths
		aEdge.setAlpha(alphaNew);
		bEdge.setAlpha(betaNew);
		aPEdge.setAlpha(alphaPNew);
		bPEdge.setAlpha(betaPNew);
		cEdge.setAlpha(betaSum);
		cPEdge.setAlpha(alphaSum); 
		lMap.put(e, l);
		lMap.put(e.getOppositeEdge(), l);

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
		Double u1 = u.getU(v1); 
		Double u2 = u.getU(v2);
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
	
	
}
