package de.varylab.discreteconformal.plugin;

import static de.varylab.discreteconformal.util.CuttingUtility.cutManifoldToDisk;
import geom3d.Point;
import geom3d.Triangle;

import java.io.FileWriter;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import java.util.Set;

import javax.swing.JOptionPane;

import no.uib.cipr.matrix.Vector;
import de.jreality.util.NativePathUtility;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Node;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AbstractAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.CalculatorException;
import de.jtem.halfedgetools.adapter.CalculatorSet;
import de.jtem.halfedgetools.adapter.type.Color;
import de.jtem.halfedgetools.algorithm.computationalgeometry.ConvexHull;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmCategory;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmPlugin;
import de.jtem.mfc.field.Complex;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.adapter.PositionAdapter;
import de.varylab.discreteconformal.heds.calculator.EdgeLengthCalculator;
import de.varylab.discreteconformal.heds.calculator.SubdivisionCalculator;
import de.varylab.discreteconformal.unwrapper.EuclideanLayout;
import de.varylab.discreteconformal.unwrapper.EuclideanUnwrapperPETSc;
import de.varylab.discreteconformal.unwrapper.Unwrapper;
import de.varylab.discreteconformal.util.AlgebraicCurveUtility;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;
import de.varylab.discreteconformal.util.Delaunay;
import de.varylab.discreteconformal.util.PathUtility;
import de.varylab.discreteconformal.util.Search;
import de.varylab.discreteconformal.util.Search.DefaultWeightAdapter;
import de.varylab.discreteconformal.util.SurgeryUtility;

public class EllipticModulusEngine extends AlgorithmPlugin {

	private static Random 
		rnd = new Random();
	
	@Override
	public AlgorithmCategory getAlgorithmCategory() {
		return AlgorithmCategory.Generator;
	}

	@Override
	public String getAlgorithmName() {
		return "Elliptic Curve";
	}
	
	
	public static void setRandomSeeed(long seeed) {
		rnd.setSeed(seeed);
	}
	
	
	public < 
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void execute(HDS h, CalculatorSet c, HalfedgeInterface hif) throws CalculatorException {
		String numString = JOptionPane.showInputDialog("Number of extra points", 0);
		if (numString == null) return;
		int extraPoints = Integer.parseInt(numString);
		CoHDS hds = hif.get(new CoHDS());
		Set<CoEdge> glueSet = new HashSet<CoEdge>();
		Set<CoEdge> cutSet = new HashSet<CoEdge>();
		generateEllipticCurve(hds, extraPoints, glueSet, cutSet);
		PathVisualizer pathVisualizer = new PathVisualizer();
		for (CoEdge e : glueSet) {
			pathVisualizer.add(e);
			pathVisualizer.add(e.getOppositeEdge());
		}
		// show the result
		AdapterSet set = new AdapterSet(pathVisualizer);
		set.add(new PositionAdapter());
		hif.set(hds, set);
	}
	
	
	public static void generateEllipticCurve(CoHDS hds, int numExtraPoints, Set<CoEdge> glueEdges, Set<CoEdge> cutEdges) {
		if (hds.numVertices() < 4) { // create regulat tetrahedron
			throw new RuntimeException("No branch point set in generateEllipticCurve()");
		} 
		CoVertex v1 = hds.getVertex(0);
		CoVertex v2 = hds.getVertex(1);
		CoVertex v3 = hds.getVertex(2);
		CoVertex v4 = hds.getVertex(3);
		for (CoEdge e : new HashSet<CoEdge>(hds.getEdges())) {
			hds.removeEdge(e);
		}
		for (CoFace f : new HashSet<CoFace>(hds.getFaces())) {
			hds.removeFace(f);
		}

		// additional points
		for (int i = 0; i < numExtraPoints; i++) {
			double[] pos = {rnd.nextGaussian(), rnd.nextGaussian(), rnd.nextGaussian()};
			CoVertex v = hds.addNewVertex();
			v.setPosition(new Point(pos));
		}
		
		// on the sphere
		for (CoVertex v : hds.getVertices()) {
			Point vec = v.getPosition();
			vec.normalize();
		}
		
		// convex hull
		ConvexHull.convexHull(hds, new SubdivisionCalculator(), 1E-8);
		int vOffset = hds.numVertices();
		int eOffset = hds.numEdges();
		HalfEdgeUtils.copy(hds, hds);
		for (int i = 0; i < vOffset; i++) {
			CoVertex v = hds.getVertex(i);
			CoVertex vc = hds.getVertex(vOffset + i); 
			Point p = v.getPosition();
			Point p2 = new Point(p);
			vc.setPosition(p2);
		}
		
		
		List<CoEdge> path1 = Search.bFS(v1, v2, new HashSet<CoVertex>());
		Set<CoVertex> path1Vertices = PathUtility.getVerticesOnPath(path1);
		List<CoEdge> path2 = Search.bFS(v3, v4, path1Vertices);
		List<CoEdge> cutPath2 = Search.bFS(v1, v3, path1Vertices);
		
		
		List<CoEdge> path1c = new LinkedList<CoEdge>();
		List<CoEdge> path2c = new LinkedList<CoEdge>();
		List<CoEdge> cutPath2c = new LinkedList<CoEdge>();
		for (CoEdge e : path1) {
			path1c.add(hds.getEdge(eOffset + e.getIndex()));
		}
		for (CoEdge e : path2) {
			path2c.add(hds.getEdge(eOffset + e.getIndex()));
		}
		for (CoEdge e : cutPath2) {
			cutPath2c.add(hds.getEdge(eOffset + e.getIndex()));
		}
		
		// generators
		cutEdges.addAll(path1);
		cutEdges.addAll(path1c);
		cutEdges.addAll(cutPath2);
		cutEdges.addAll(cutPath2c);
		
		SurgeryUtility.cutAndGluePaths(path1, path1c);
		SurgeryUtility.cutAndGluePaths(path2, path2c);
		
		glueEdges.addAll(path1);
		glueEdges.addAll(path2);
		glueEdges.addAll(path1c);
		glueEdges.addAll(path2c);
	}
	
	
	
	@Color
	private static class PathVisualizer extends AbstractAdapter<double[]> {

		private Set<CoEdge>
			edges = new HashSet<CoEdge>();
		private double[]
		    pathColor = {1, 0, 0},
		    defaultColor = {1, 1, 1};
		
		public PathVisualizer() {
			super(double[].class, true, false);
		}
		
		@Override
		public <N extends Node<?, ?, ?>> boolean canAccept(Class<N> nodeClass) {
			return nodeClass == CoEdge.class;
		}

		@Override
		public double getPriority() {
			return 10;
		}
		
		@Override
		public <
			V extends de.jtem.halfedge.Vertex<V,E,F>, 
			E extends de.jtem.halfedge.Edge<V,E,F>, 
			F extends de.jtem.halfedge.Face<V,E,F>
		> double[] getE(E e, de.jtem.halfedgetools.adapter.AdapterSet a) {
			if (edges.contains(e)) {
				return pathColor;
			} else {
				return defaultColor;
			}
		};
		
		public void add(CoEdge edge) {
			edges.add(edge);
		}

	}

	
	public static Complex calculateModulus(CoHDS hds) {

		Unwrapper unwrapper = new EuclideanUnwrapperPETSc();
		unwrapper.setGradientTolerance(1E-7);
		unwrapper.setMaxIterations(500);
		Vector u = null;
		try {
			u = unwrapper.unwrap(hds);
		} catch (Exception e) {
			e.printStackTrace();
			return new Complex();
		}
		DefaultWeightAdapter<CoEdge> w = new DefaultWeightAdapter<CoEdge>();
		CoVertex cutRoot = hds.getVertex(0);
		CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo = cutManifoldToDisk(hds, cutRoot, w);
		EuclideanLayout.doLayout(hds, u);
		
		return AlgebraicCurveUtility.calculateCutModulus(cutInfo);
	}
	
	
	public static double getMaxCircumRadius(CoHDS hds) {
		double max = 0;
		for (CoFace f : hds.getFaces()) {
			Triangle t = f.toTriangle();
			double rad = t.getCircumRadius();
			if (rad > max) {
				max = rad;
			}
		}
		return max;
	}
	
	
	
	
	public static void main(String[] args) throws Exception {
		NativePathUtility.set("native");
		
		FileWriter fwRe = new FileWriter("data/tauReErrEquator01.dat");
		FileWriter fwIm = new FileWriter("data/tauImErrEquator01.dat");
		FileWriter fwAbs = new FileWriter("data/tauAbsErrEquator01.dat");
		for (int i = 1; i <= 5000; i++) {
			System.out.println(i + " extra points --------------------");
			Set<CoEdge> glueSet = new HashSet<CoEdge>();
			Set<CoEdge> cutSet = new HashSet<CoEdge>();
			Complex tau = null;
//			double maxCircRad = 0.0;
			
			CoHDS hds = new CoHDS();
			CoVertex v1 = hds.addNewVertex();
			CoVertex v2 = hds.addNewVertex();
			CoVertex v3 = hds.addNewVertex();
			CoVertex v4 = hds.addNewVertex();
			v1.getPosition().set(1, 0, 0);
			v2.getPosition().set(0, 1, 0);
			v3.getPosition().set(-1, 0, 0);
			v4.getPosition().set(0, -1, 0);
//			double a = 2 * sin(PI / 4);
//			v1.getPosition().set(a, 1, 0);
//			v2.getPosition().set(-a, 1, 0);
//			v3.getPosition().set(0, -1, -a);
//			v4.getPosition().set(0, -1, a);
			
			try {
				generateEllipticCurve(hds, i, glueSet, cutSet);
				Delaunay.constructDelaunay(hds, new EdgeLengthCalculator());
				tau = calculateModulus(hds);
//				maxCircRad = getMaxCircumRadius(hds);
			} catch (Exception e) {
				continue;
			}
			double reErr = Math.abs(Math.abs(tau.re) - 0);
			double imErr = Math.abs(Math.abs(tau.im) - 1);
			double absErr = Math.abs(tau.abs() - 1);
			System.out.println("Re err " + reErr);
			System.out.println("Im err " + imErr);
			System.out.println("Abs err " + absErr);
			fwAbs.write(i + "\t" + absErr + "\n");
			fwRe.write(i + "\t" + reErr + "\n");
			fwIm.write(i + "\t" + imErr + "\n");
			fwAbs.flush();
			fwRe.flush();
			fwIm.flush();
		}
		fwAbs.close();
		fwRe.close();
		fwIm.close();
	}

}

