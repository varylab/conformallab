package de.varylab.discreteconformal.plugin;

import static javax.swing.JOptionPane.showInputDialog;

import java.awt.EventQueue;
import java.util.LinkedList;
import java.util.Set;
import java.util.TreeSet;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.algorithm.computationalgeometry.ConvexHull;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmCategory;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmPlugin;
import de.jtem.halfedgetools.selection.Selection;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.SphereUtility;

public class SphericalEqualizerPlugin extends AlgorithmPlugin {

	private String
		iterationsString = null;
	
	@Override
	public AlgorithmCategory getAlgorithmCategory() {
		return AlgorithmCategory.Geometry;
	}

	@Override
	public String getAlgorithmName() {
		return "Spherical Equalizer";
	}

	@Override
	public <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>, 
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void execute(HDS hds, AdapterSet a, HalfedgeInterface hi) throws Exception {
		Runnable r = new Runnable() {
			@Override
			public void run() {
				iterationsString = showInputDialog(getOptionParent(), "Number of iterations", 20);				
			}
		};
		EventQueue.invokeAndWait(r);
		if (iterationsString == null) return;
		int iterations = Integer.parseInt(iterationsString);
		CoHDS chds = hi.get(new CoHDS());
		Selection sel = hi.getSelection();
		Set<CoVertex> fivedVertices = new TreeSet<CoVertex>(sel.getVertices(chds));
		SphereUtility.equalizeSphereVertices(chds, fivedVertices, iterations, 1E-6);
		for (CoEdge e : new LinkedList<CoEdge>(chds.getEdges())) {
			chds.removeEdge(e);
		}
		for (CoFace f : new LinkedList<CoFace>(chds.getFaces())) {
			chds.removeFace(f);
		}
		ConvexHull.convexHull(chds, a, 1E-6);
		hi.update();
	}

}
