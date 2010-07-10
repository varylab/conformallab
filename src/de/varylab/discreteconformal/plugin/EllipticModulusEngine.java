package de.varylab.discreteconformal.plugin;

import geom3d.Point;
import geom3d.Vector;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import java.util.Set;

import de.jreality.plugin.JRViewer;
import de.jreality.plugin.JRViewer.ContentType;
import de.jreality.plugin.basic.Content;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.Node;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AbstractAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Color;
import de.jtem.halfedgetools.algorithm.computationalgeometry.ConvexHull;
import de.jtem.halfedgetools.plugin.GeneratorPlugin;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.adapter.PositionAdapter;
import de.varylab.discreteconformal.heds.calculator.SubdivisionCalculator;
import de.varylab.discreteconformal.util.PathUtility;
import de.varylab.discreteconformal.util.Search;

public class EllipticModulusEngine extends GeneratorPlugin {

	private Random 
		rnd = new Random();
	
	@Override
	protected void generate(Content content, HalfedgeInterface hif) {
		CoHDS hds = new CoHDS();
		
		// branch points
		CoVertex v1 = hds.addNewVertex();
		CoVertex v2 = hds.addNewVertex();
		CoVertex v3 = hds.addNewVertex();
		CoVertex v4 = hds.addNewVertex();
		double[] pos1 = {rnd.nextGaussian(), rnd.nextGaussian(), rnd.nextGaussian()};
		double[] pos2 = {rnd.nextGaussian(), rnd.nextGaussian(), rnd.nextGaussian()};
		double[] pos3 = {rnd.nextGaussian(), rnd.nextGaussian(), rnd.nextGaussian()};
		double[] pos4 = {rnd.nextGaussian(), rnd.nextGaussian(), rnd.nextGaussian()};
		v1.setPosition(new Point(pos1));
		v2.setPosition(new Point(pos2));
		v3.setPosition(new Point(pos3));
		v4.setPosition(new Point(pos4));
		
		// additional points
		for (int i = 0; i < 100; i++) {
			double[] pos = {rnd.nextGaussian(), rnd.nextGaussian(), rnd.nextGaussian()};
			CoVertex v = hds.addNewVertex();
			v.setPosition(new Point(pos));
		}
		
		// on the sphere
		for (CoVertex v : hds.getVertices()) {
			Vector vec = v.getPosition();
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
			p2.get()[0] += 2;
			vc.setPosition(p2);
		}
		
		
		List<CoEdge> path1 = Search.bFS(v1, v2, new HashSet<CoVertex>());
		Set<CoVertex> path1Vertices = PathUtility.getVerticesOnPath(path1);
		List<CoEdge> path2 = Search.bFS(v3, v4, path1Vertices);
		
		List<CoEdge> path1c = new LinkedList<CoEdge>();
		List<CoEdge> path2c = new LinkedList<CoEdge>();		
		for (CoEdge e : path1) {
			path1c.add(hds.getEdge(eOffset + e.getIndex()));
		}
		for (CoEdge e : path2) {
			path2c.add(hds.getEdge(eOffset + e.getIndex()));
		}
		
		
		
		
		
		PathVisualizer pathVisualizer = new PathVisualizer();
		for (CoEdge e : path1) {
			pathVisualizer.add(e);
			pathVisualizer.add(e.getOppositeEdge());
		}
		for (CoEdge e : path2) {
			pathVisualizer.add(e);
			pathVisualizer.add(e.getOppositeEdge());
		}
		
		for (CoEdge e : path1c) {
			pathVisualizer.add(e);
			pathVisualizer.add(e.getOppositeEdge());
		}
		for (CoEdge e : path2c) {
			pathVisualizer.add(e);
			pathVisualizer.add(e.getOppositeEdge());
		}
		
		// show the result
		AdapterSet set = new AdapterSet(pathVisualizer);
		set.add(new PositionAdapter());
		hif.set(hds, set);
	}
	
	
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> void cutAndGluePaths(List<E> p1, List<E> p2) {
		if (p1.size() != p2.size()) {
			throw new IllegalArgumentException("Paths of different lengths in cutAndGluePaths()");
		}
		V start = p1.get(0).getStartVertex();
		V end = p1.get(p1.size() - 1).getTargetVertex();

		V startC = p2.get(0).getStartVertex();
		V endC = p2.get(p2.size() - 1).getTargetVertex();
		
		for (int i = 0; i < p1.size(); i++) {
			E e = p1.get(i);
			
		}
	}
	
	
	
	
	@Color
	private class PathVisualizer extends AbstractAdapter<double[]> {

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

	@Override
	protected String[] getMenuPath() {
		return new String[] {"Algebraic Curves"};
	}
	
	public static void main(String[] args) {
		JRViewer v = new JRViewer();
		v.addBasicUI();
		v.addContentUI();
		v.addContentSupport(ContentType.CenteredAndScaled);
		v.registerPlugin(new EllipticModulusEngine());
		v.startup();
	}

}

