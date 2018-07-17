package de.varylab.discreteconformal.util;

import static de.varylab.discreteconformal.util.PathUtility.getVerticesOnPath;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import de.jreality.math.Pn;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.algorithm.computationalgeometry.ConvexHull;
import de.varylab.discreteconformal.ConformalAdapterSet;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;


public class HyperellipticUtility {

	
	/**
	 * Generates a two-sheeted version of hds branched at the vertices with 
	 * indices specified by branchVertices.
	 * @param hds the base data structure for one sheet
	 * @param useConvexHull use a convex hull algorithm to create new edges and faces
	 * @param numExtraPoints add a number of extra points on S^2
	 * @param glueEdges the set of edges that has been used to cut and re-glue the two sheets
	 * @param branchVertices set indices of the vertices that will be used as branch vertices
	 * @return the set of branch vertices in the result
	 */
	public static Set<CoVertex> generateHyperellipticImage(
		CoHDS hds, 
		boolean useConvexHull,
		int numExtraPoints, 
		boolean projectToSphere,
		Set<CoEdge> glueEdges, 
		int... branchVertices
	) {
		if (hds.numVertices() < 4) {
			throw new RuntimeException("Not enough vertices in generateHyperellipticImage()");
		} 
		if (branchVertices.length < 3) {
			throw new RuntimeException("Not enough branch points specified in generateHyperellipticImage()");
		}
		if (!useConvexHull && numExtraPoints > 0) {
			throw new RuntimeException("connot have extra vertices without convex hull generation in generateHyperellipticImage()");
		}
		List<CoVertex> branchList = new LinkedList<CoVertex>();
		for (int i : branchVertices) {
			CoVertex v = hds.getVertex(i);
			branchList.add(v);
		}
		// add point at infinity if necessary
		if (branchVertices.length % 2 != 0) {
			CoVertex v = hds.addNewVertex();
			v.P = new double[] {0.0, 0.0, 1.0, 1.0};
			branchList.add(v);
		}
		if (numExtraPoints > 0) {
			// additional points
			Random rnd = new Random();
			for (int i = 0; i < numExtraPoints; i++) {
				CoVertex v = hds.addNewVertex();
				v.P = new double[] {rnd.nextGaussian(), rnd.nextGaussian(), rnd.nextGaussian(), 1.0};
			}
		}
		if (projectToSphere) {
			// on the sphere
			for (CoVertex v : hds.getVertices()) {
				Pn.setToLength(v.P, v.P, 1.0, Pn.EUCLIDEAN);
			}
		}
		if (useConvexHull) {
			// clean up and create new edges and faces of the convex hull
			for (CoEdge e : new HashSet<CoEdge>(hds.getEdges())) {
				hds.removeEdge(e);
			}
			for (CoFace f : new HashSet<CoFace>(hds.getFaces())) {
				hds.removeFace(f);
			}
			ConformalAdapterSet a = new ConformalAdapterSet();
			ConvexHull.convexHull(hds, a, 1E-8);
		}

		int vOffset = hds.numVertices();
		int eOffset = hds.numEdges();
		HalfEdgeUtils.copy(hds, hds);
		Map<CoVertex, CoVertex> sheetVertexMap = new HashMap<>(); 
		for (int i = 0; i < vOffset; i++) {
			CoVertex v = hds.getVertex(i);
			CoVertex vc = hds.getVertex(vOffset + i); 
			double[] p = v.P;
			vc.P = p.clone();
			sheetVertexMap.put(v, vc);
			sheetVertexMap.put(vc, v);
		}
		
		glueEdges.clear();
		Set<CoVertex> avoid = new HashSet<CoVertex>(branchList);
		Set<CoVertex> originalBranchVertices = new HashSet<>(branchList);
		while (!branchList.isEmpty()) {
			CoVertex start = branchList.remove(0);
			CoVertex end = branchList.remove(0);
			avoid.remove(start);
			avoid.remove(end);
			List<CoEdge> path1 = Search.bFS(start, end, avoid);
			List<CoEdge> path1c = new LinkedList<CoEdge>();
			for (CoEdge e : path1) {
				int sheetOffset = e.getIndex() >= eOffset ? -eOffset : eOffset;
				path1c.add(hds.getEdge(sheetOffset + e.getIndex()));
			}
			SurgeryUtility.cutAndGluePaths(path1, path1c);
			Set<CoVertex> path1Vertices = getVerticesOnPath(path1);
			avoid.addAll(path1Vertices);
			glueEdges.addAll(path1);
			glueEdges.addAll(path1c);
		}
		
		// create result branch vertices
		Set<CoVertex> resultBranchVertices = new HashSet<>();
		for (CoVertex v : originalBranchVertices) {
			CoVertex bv = v;
			if (!v.isValid()) {
				bv = sheetVertexMap.get(v);
			}
			resultBranchVertices.add(bv);
		}
		return resultBranchVertices;
	}
	
}
