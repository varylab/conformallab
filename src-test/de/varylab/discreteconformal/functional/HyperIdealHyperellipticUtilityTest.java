package de.varylab.discreteconformal.functional;

import static java.lang.Math.PI;
import static java.lang.Math.toDegrees;

import java.util.List;

import org.junit.Assert;
import org.junit.Test;

import de.jreality.math.Pn;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.mfc.field.Complex;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.math.ComplexUtility;

public class HyperIdealHyperellipticUtilityTest extends Assert {

	@Test
	public void testCalculateCircleIntersections() throws Exception {
		CoHDS hds = new CoHDS();
		List<CoVertex> verts = hds.addNewVertices(4);
		CoVertex a = verts.get(0);
		CoVertex b = verts.get(1);
		CoVertex c = verts.get(2);
		CoVertex d = verts.get(3);
		a.P = Pn.homogenize(null, ComplexUtility.inverseStereographic(new Complex(2, 0)));
		b.P = Pn.homogenize(null, ComplexUtility.inverseStereographic(new Complex(2, 1)));
		c.P = Pn.homogenize(null, ComplexUtility.inverseStereographic(new Complex(1, 0)));
		d.P = Pn.homogenize(null, ComplexUtility.inverseStereographic(new Complex(3, 0)));
		CoFace fl = HalfEdgeUtils.constructFaceByVertices(hds, a, b, c).getLeftFace();
		CoFace fr = HalfEdgeUtils.constructFaceByVertices(hds, a, d, b).getLeftFace();
		CoEdge e = HalfEdgeUtils.findEdgeBetweenFaces(fl, fr);
		HyperIdealHyperellipticUtility.calculateCircleIntersections(hds);
		assertEquals(PI/2, e.getTheta(), 1e-8);
	}
	
	@Test
	public void testCalculateCircleIntersectionsInfinite() throws Exception {
		CoHDS hds = new CoHDS();
		List<CoVertex> verts = hds.addNewVertices(4);
		CoVertex a = verts.get(0);
		CoVertex b = verts.get(1);
		CoVertex c = verts.get(2);
		CoVertex d = verts.get(3);
		a.P = new double[]{0,1,0,1};
		b.P = new double[]{0,0,1,1};
		c.P = new double[]{-1,0,0,1};
		d.P = new double[]{1,0,0,1};
		CoFace fl = HalfEdgeUtils.constructFaceByVertices(hds, a, b, c).getLeftFace();
		CoFace fr = HalfEdgeUtils.constructFaceByVertices(hds, a, d, b).getLeftFace();
		CoEdge e = HalfEdgeUtils.findEdgeBetweenFaces(fl, fr);
		HyperIdealHyperellipticUtility.calculateCircleIntersections(hds);
		assertEquals(PI/2, e.getTheta(), 1e-8);
	}
	
	@Test
	public void testLawsonHyperellipticAngles() throws Exception {
		CoHDS hds = HyperIdealGenerator.createLawsonHyperelliptic();
		for (CoEdge e : hds.getPositiveEdges()) {
				int i = e.getStartVertex().getIndex();
				int j = e.getTargetVertex().getIndex();
				// an edge not connected to a branch point
				if (i != 0 && i != 1 && i != 2 && i != 3 && i != 6 && i != 7 &&
					j != 0 && j != 1 && j != 2 && j != 3 && j != 6 && j != 7) {
					Assert.assertEquals(150, toDegrees(e.getTheta()), 1E-12);
				} else
				// edge from a branch point to the north or south pole
				if (i == 4 || i == 5 || i == 8 || i == 9 || j == 4 || j == 5 || j == 8 || j == 9) {
					Assert.assertEquals(150, toDegrees(e.getTheta()), 1E-12);					
				} 
				// edge from a branch point to the mid point on the equator
				else {
					Assert.assertEquals(30, toDegrees(e.getTheta()), 1E-12);					
				}
		}
	}
	
}
