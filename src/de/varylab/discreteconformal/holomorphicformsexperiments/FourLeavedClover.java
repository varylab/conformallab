package de.varylab.discreteconformal.holomorphicformsexperiments;

import java.util.List;

import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.algorithm.triangulation.MappedLengthAdapter;
import de.jtem.mfc.field.Complex;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;

/**
 * Construction taken from arXiv:1502.05381v2 p. 38
 * 
 * @author sechel
 *
 */
public class FourLeavedClover {

	private static Complex TAU = new Complex(5./4., 2./3.);
	private static Complex I = new Complex(0, 1);
	private static double l0 = 1.0;
	private static double l1 = I.minus(I.times(TAU)).abs();
	private static double l2 = TAU.abs();
	
	private static List<CoEdge> createFace(CoHDS S, MappedLengthAdapter l) {
		CoFace f = S.addNewFace();
		CoVertex v0 = S.getVertex(0);
		CoVertex v1 = S.getVertex(1);
		CoVertex v2 = S.getVertex(2);
		CoVertex v3 = S.getVertex(3);
		List<CoEdge> E = S.addNewEdges(6);
		CoEdge e0 = E.get(0);
		CoEdge e1 = E.get(1);
		CoEdge e2 = E.get(2);
		CoEdge e3 = E.get(3);
		CoEdge e4 = E.get(4);
		CoEdge e5 = E.get(5);
		e0.linkNextEdge(e1);
		e0.setTargetVertex(v1);
		e0.setLeftFace(f);
		l.setE(e0, FourLeavedClover.l0, null);
		e1.linkNextEdge(e2);
		e1.setTargetVertex(v2);
		e1.setLeftFace(f);
		l.setE(e1, FourLeavedClover.l1, null);
		e2.linkNextEdge(e3);
		e2.setTargetVertex(v1);
		e2.setLeftFace(f);
		l.setE(e2, FourLeavedClover.l1, null);
		e3.linkNextEdge(e4);
		e3.setTargetVertex(v3);
		e3.setLeftFace(f);
		l.setE(e3, FourLeavedClover.l2, null);
		e4.linkNextEdge(e5);
		e4.setTargetVertex(v1);
		e4.setLeftFace(f);
		l.setE(e4, FourLeavedClover.l2, null);
		e5.linkNextEdge(e0);
		e5.setTargetVertex(v0);
		e5.setLeftFace(f);
		l.setE(e5, FourLeavedClover.l0, null);
		return E;
	}
	
	public static void main(String[] args) {
		CoHDS S = new CoHDS();
		S.addNewVertices(4);
		MappedLengthAdapter l = new MappedLengthAdapter();
		List<CoEdge> f0 = FourLeavedClover.createFace(S, l);
		List<CoEdge> f1 = FourLeavedClover.createFace(S, l);
		List<CoEdge> f2 = FourLeavedClover.createFace(S, l);
		List<CoEdge> f3 = FourLeavedClover.createFace(S, l);
		f3.get(3).linkOppositeEdge(f0.get(4)); // a
		f2.get(3).linkOppositeEdge(f3.get(4)); // b
		f1.get(3).linkOppositeEdge(f2.get(4)); // c
		f0.get(3).linkOppositeEdge(f1.get(4)); // d
		f3.get(1).linkOppositeEdge(f0.get(2)); // e
		f2.get(1).linkOppositeEdge(f3.get(2)); // f
		f1.get(1).linkOppositeEdge(f2.get(2)); // g
		f0.get(1).linkOppositeEdge(f1.get(2)); // h
		f3.get(5).linkOppositeEdge(f0.get(0));
		f2.get(5).linkOppositeEdge(f3.get(0));
		f1.get(5).linkOppositeEdge(f2.get(0));
		f0.get(5).linkOppositeEdge(f1.get(0));
		int g = HalfEdgeUtils.getGenus(S);
		System.out.println(S);
		System.out.println("Genus " + g);
		
	}

}
