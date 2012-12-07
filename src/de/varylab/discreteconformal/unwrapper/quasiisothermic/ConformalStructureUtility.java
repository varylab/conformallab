package de.varylab.discreteconformal.unwrapper.quasiisothermic;

import static de.jtem.jpetsc.InsertMode.INSERT_VALUES;
import static java.lang.Math.abs;
import static java.lang.Math.log;
import static java.lang.Math.sqrt;

import java.util.HashMap;
import java.util.Map;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.jpetsc.KSP;
import de.jtem.jpetsc.Mat;
import de.jtem.jpetsc.MatStructure;
import de.jtem.jpetsc.PETSc;
import de.jtem.jpetsc.Vec;

public class ConformalStructureUtility {

	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> Map<E, Double> calculatePseudoConformalStructure(HDS hds, Map<E, Double> alphaMap) {
		Map<E, Double> lcrPseudo = new HashMap<E, Double>();
		for (E e : hds.getPositiveEdges()) {
			if (e.getLeftFace() == null) {
				e = e.getOppositeEdge();
			}
			E eim = e.getNextEdge();
			E ejk = e.getOppositeEdge().getNextEdge();
			E emj = e.getPreviousEdge();
			E eki = e.getOppositeEdge().getPreviousEdge();
			double lim = QuasiisothermicUtility.getEdgeLength(eim, 1.0, alphaMap);
			double ljk = QuasiisothermicUtility.getEdgeLength(ejk, 1.0, alphaMap);
			double lmj = QuasiisothermicUtility.getEdgeLength(emj, lim, alphaMap);
			double lki = QuasiisothermicUtility.getEdgeLength(eki, ljk, alphaMap);
			if (e.getRightFace() == null) {
				ljk = 1.0;
				lki = 1.0;
			}
			double elcrPseudo = (lim*ljk) / (lki*lmj);			
			lcrPseudo.put(e, elcrPseudo);
			lcrPseudo.put(e.getOppositeEdge(), elcrPseudo);
		}
		return lcrPseudo;
	}
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> double calculateVertexCrossRatioProduct(V v, Map<E, Double> lcrMap) {
		double p = 1;
		for (E e : HalfEdgeUtils.incomingEdges(v)) {
			p *= lcrMap.get(e);
		}
		return p;
	}

	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> Map<E, Double> calculateConformalStructure(HDS hds, Map<E, Double> lcrPseudoMap) {
		int numBoundaryEdge = HalfEdgeUtils.boundaryEdges(hds).size();
		int dim = hds.numEdges() / 2 - numBoundaryEdge;
		Map<E, Integer> domainIndices = new HashMap<E, Integer>();
		Map<V, Integer> imageIndices = new HashMap<V, Integer>();
		
		// create indices
		int i = 0;
		int[] nz = new int[dim];
		for (V v : hds.getVertices()) {
			if (HalfEdgeUtils.isBoundaryVertex(v)) continue;
			nz[i] = HalfEdgeUtils.incomingEdges(v).size();
			imageIndices.put(v, i++);
		}
		int j = 0;
		for (E e : hds.getPositiveEdges()) {
			if (HalfEdgeUtils.isBoundaryEdge(e)) continue;
			domainIndices.put(e, j);
			domainIndices.put(e.getOppositeEdge(), j++);
		}		
		
		int dimDomain = j;

		// initialize matrices
		Vec r = new Vec(dimDomain);
		Mat A = Mat.createSeqAIJ(dimDomain, dimDomain, -1, nz);
		A.zeroEntries();
		r.zeroEntries();
		for (V v : hds.getVertices()) {
			if (HalfEdgeUtils.isBoundaryVertex(v)) continue;
			double p = calculateVertexCrossRatioProduct(v, lcrPseudoMap);
			i = imageIndices.get(v);
			r.setValue(i, -log(p), INSERT_VALUES);
			for (E e : HalfEdgeUtils.incomingEdges(v)) {
				j = domainIndices.get(e);
				A.setValue(i, j, 1.0, INSERT_VALUES);
			}
		}
		A.assemble();

		// solve
		Vec g = new Vec(dimDomain);
		g.zeroEntries();
		KSP ksp = KSP.create();
		ksp.setOptionsPrefix("cs_");
		PETSc.optionsSetValue("-cs_ksp_type", "lsqr");
		ksp.setFromOptions();
		ksp.setOperators(A, A, MatStructure.SAME_NONZERO_PATTERN);
		ksp.solve(r, g);
		System.out.println("conformal structure calculation ------------");
		System.out.println("reason: " + ksp.getConvergedReason());
		System.out.println("iteration number: " + ksp.getIterationNumber());
		System.out.println("ksp residual: " + ksp.getResidualNorm());	
		System.out.println("solution: " + g);
		Map<E, Double> lcrResult = new HashMap<E, Double>();
		for (E e : hds.getPositiveEdges()) {
			if (HalfEdgeUtils.isBoundaryEdge(e)) {
				double lcr = lcrPseudoMap.get(e);
				lcrResult.put(e, lcr);
				lcrResult.put(e.getOppositeEdge(), lcr);
				continue;
			}
			i = domainIndices.get(e);
			double lcrPseudo = lcrPseudoMap.get(e);
			double ge = g.getValue(i);
			double lcr = lcrPseudo * Math.exp(ge);
			lcrResult.put(e, lcr);
			lcrResult.put(e.getOppositeEdge(), lcr);
		}
		
		return lcrResult;
	}
	
	
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> Map<E, Double> lengthsFromCrossRatios(HDS hds, Map<E, Double> lcrMap) {
		Map<E, Double> lMap = new HashMap<E, Double>();
		// define lengths from cross-ratios
		Map<E, Double> aMap = new HashMap<E, Double>();
		for (V v : hds.getVertices()) {
			double ai = 100.0;
			for (E e : HalfEdgeUtils.incomingEdges(v)) {
				double q = lcrMap.get(e);
				ai *= q;
				aMap.put(e, ai);
			}
		}
		for (E e : hds.getEdges()) {
			if (e.getLeftFace() == null) continue;
			double ai = aMap.get(e);
			double aj = aMap.get(e.getPreviousEdge());
			double l = sqrt(1 / (ai * aj));
			lMap.put(e, l);
			if (e.getRightFace() == null) {
				lMap.put(e.getOppositeEdge(), l);	
			}
		}
		// check opposite lengths
		for (E e : hds.getEdges()) {
			double l = lMap.get(e);
			double lo = lMap.get(e.getOppositeEdge());
			assert abs(l - lo) < 1E-8 : "lengths";
		}
		return lMap;
	}
	
	
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> Map<V, Double> boundaryAnglesFromAlphas(HDS hds, Map<E, Double> alphaMap) {
		Map<V, Double> thetaMap = new HashMap<V, Double>();
		for (V v : HalfEdgeUtils.boundaryVertices(hds)) {
			if (!HalfEdgeUtils.isBoundaryVertex(v)) continue;
			double theta = 0.0;
			for (E e : HalfEdgeUtils.incomingEdges(v)) {
				if (e.getLeftFace() == null) continue;
				double a0 = alphaMap.get(e);
				double a1 = alphaMap.get(e.getNextEdge());
				double a2 = alphaMap.get(e.getPreviousEdge());
				double beta = QuasiisothermicUtility.calculateBeta(a0, a1, a2);
				theta += beta;
			}
			thetaMap.put(v, theta);
		}
		return thetaMap;
	}
	
}
