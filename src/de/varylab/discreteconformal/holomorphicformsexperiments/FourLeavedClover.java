package de.varylab.discreteconformal.holomorphicformsexperiments;

import static de.jtem.halfedgetools.algorithm.triangulation.Delaunay.constructDelaunay;
import static de.varylab.discreteconformal.util.DiscreteRiemannUtility.getHolomorphicForms;
import static de.varylab.discreteconformal.util.LaplaceUtility.calculateCotanWeights;

import java.util.List;

import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.generic.UndirectedEdgeIndex;
import de.jtem.halfedgetools.adapter.type.Length;
import de.jtem.halfedgetools.algorithm.triangulation.MappedLengthAdapter;
import de.jtem.mfc.field.Complex;
import de.varylab.discreteconformal.adapter.MappedWeightAdapter;
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

	private static Complex TAU = new Complex(5. / 4., 2. / 3.);
	private static Complex I = new Complex(0, 1);
	private static double l0 = 1.0;
	private static double l1 = I.minus(I.times(TAU)).abs();
	private static double l2 = TAU.abs();
	private static double l3 = Math.sqrt(2);
	private static double l4 = Math.sqrt(2) * FourLeavedClover.l1;
	private static double l5 = Math.sqrt(2) * FourLeavedClover.l2;
	private static Loop loop = new Loop();

	@Length
	private static class HighPriorityLengthAdapter extends MappedLengthAdapter {
		
		@Override
		public double getPriority() {
			return 100;
		}
		
	}
	
	private static List<CoEdge> createFace(CoHDS S, MappedLengthAdapter l) {
		CoFace f0 = S.addNewFace();
		CoFace f1 = S.addNewFace();
		CoFace f2 = S.addNewFace();
		CoFace f3 = S.addNewFace();
		CoVertex v0 = S.getVertex(0);
		CoVertex v1 = S.getVertex(1);
		CoVertex v2 = S.getVertex(2);
		CoVertex v3 = S.getVertex(3);
		List<CoEdge> E = S.addNewEdges(12);
		CoEdge e0 = E.get(0);
		CoEdge e1 = E.get(1);
		CoEdge e2 = E.get(2);
		CoEdge e3 = E.get(3);
		CoEdge e4 = E.get(4);
		CoEdge e5 = E.get(5);
		CoEdge e6 = E.get(6);
		CoEdge e7 = E.get(7);
		CoEdge e8 = E.get(8);
		CoEdge e9 = E.get(9);
		CoEdge e10 = E.get(10);
		CoEdge e11 = E.get(11);
		// outer edges
		e0.linkNextEdge(e9);
		e0.setTargetVertex(v1);
		e0.setLeftFace(f0);
		l.setE(e0, FourLeavedClover.l0, null);
		e1.linkNextEdge(e2);
		e1.setTargetVertex(v2);
		e1.setLeftFace(f1);
		l.setE(e1, FourLeavedClover.l1, null);
		e2.linkNextEdge(e10);
		e2.setTargetVertex(v1);
		e2.setLeftFace(f1);
		l.setE(e2, FourLeavedClover.l1, null);
		e3.linkNextEdge(e4);
		e3.setTargetVertex(v3);
		e3.setLeftFace(f2);
		l.setE(e3, FourLeavedClover.l2, null);
		e4.linkNextEdge(e11);
		e4.setTargetVertex(v1);
		e4.setLeftFace(f2);
		l.setE(e4, FourLeavedClover.l2, null);
		e5.linkNextEdge(e0);
		e5.setTargetVertex(v0);
		e5.setLeftFace(f0);
		l.setE(e5, FourLeavedClover.l0, null);
		// inner edges
		e6.linkNextEdge(e7);
		e6.setTargetVertex(v1);
		e6.setLeftFace(f3);
		e6.linkOppositeEdge(e9);
		l.setE(e6, FourLeavedClover.l3, null);
		e7.linkNextEdge(e8);
		e7.setTargetVertex(v1);
		e7.setLeftFace(f3);
		e7.linkOppositeEdge(e10);
		l.setE(e7, FourLeavedClover.l4, null);
		e8.linkNextEdge(e6);
		e8.setTargetVertex(v1);
		e8.setLeftFace(f3);
		e8.linkOppositeEdge(e11);
		l.setE(e8, FourLeavedClover.l5, null);
		e9.linkNextEdge(e5);
		e9.setTargetVertex(v1);
		e9.setLeftFace(f0);
		l.setE(e9, FourLeavedClover.l3, null);
		e10.linkNextEdge(e1);
		e10.setTargetVertex(v1);
		e10.setLeftFace(f1);
		l.setE(e10, FourLeavedClover.l4, null);
		e11.linkNextEdge(e3);
		e11.setTargetVertex(v1);
		e11.setLeftFace(f2);
		l.setE(e11, FourLeavedClover.l5, null);
		
		return E;
	}

	public static void main(String[] args) {
		CoHDS Sstart = new CoHDS();
		Sstart.addNewVertices(4);
		MappedLengthAdapter l = new HighPriorityLengthAdapter();
		List<CoEdge> f0 = FourLeavedClover.createFace(Sstart, l);
		List<CoEdge> f1 = FourLeavedClover.createFace(Sstart, l);
		List<CoEdge> f2 = FourLeavedClover.createFace(Sstart, l);
		List<CoEdge> f3 = FourLeavedClover.createFace(Sstart, l);
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
		
		AdapterSet a = new AdapterSet(l, new UndirectedEdgeIndex());
		CoHDS S2 = new CoHDS();
		FourLeavedClover.loop.subdivide(Sstart, S2, a);
		CoHDS S3 = new CoHDS();
		FourLeavedClover.loop.subdivide(S2, S3, a);
		CoHDS S4 = new CoHDS();
		FourLeavedClover.loop.subdivide(S3, S4, a);
		CoHDS S5 = new CoHDS();
		FourLeavedClover.loop.subdivide(S4, S5, a);
		CoHDS S = new CoHDS();
		FourLeavedClover.loop.subdivide(S5, S, a);
		
		int g = HalfEdgeUtils.getGenus(S);
		boolean valid = HalfEdgeUtils.isValidSurface(S, true);
		System.out.println(S);
		System.out.println("Valid: " + valid);
		System.out.println("Genus " + g);		
		
		constructDelaunay(S, a);
		MappedWeightAdapter cotanWeights = calculateCotanWeights(S, a);
		System.out.println(cotanWeights.getWeights());
		a.add(cotanWeights);
		Complex[][] dhs = getHolomorphicForms(S, a, null);
		System.out.println(dhs.length);
		System.out.println(dhs[0].length);
	}

}
