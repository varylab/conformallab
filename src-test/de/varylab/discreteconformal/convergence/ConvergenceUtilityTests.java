package de.varylab.discreteconformal.convergence;

import static de.jtem.halfedge.util.HalfEdgeUtils.constructFaceByVertices;
import static de.varylab.discreteconformal.convergence.ConvergenceUtility.getMaxMeanSumScaleInvariantCircumRadius;
import static de.varylab.discreteconformal.convergence.ConvergenceUtility.getTextureTriangleArea;
import static java.lang.Math.sin;

import org.junit.Assert;
import org.junit.Test;

import de.jtem.halfedge.util.HalfEdgeUtils;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;

public class ConvergenceUtilityTests {

	@Test
	public void testGetTextureCircumRadius() throws Exception {
		CoHDS hds = new CoHDS();
		HalfEdgeUtils.addTetrahedron(hds);
		CoFace f = hds.getFace(0);
		
		CoVertex A = f.getBoundaryEdge().getStartVertex();
		CoVertex B = f.getBoundaryEdge().getNextEdge().getStartVertex();
		CoVertex C = f.getBoundaryEdge().getPreviousEdge().getStartVertex();
		
		A.T = new double[] {0,0,0,1};
		B.T = new double[] {1,0,0,1};
		C.T = new double[] {0,1,0,1};
		
		double r = ConvergenceUtility.getTextureCircumCircleRadius(f);
		Assert.assertEquals(Math.sqrt(2)/2, r, 1E-10);
		
		B.T = new double[] {0.5,0.5,0,1};
		r = ConvergenceUtility.getTextureCircumCircleRadius(f);
		Assert.assertEquals(0.5, r, 1E-10);
	}
	
	@Test
	public void testGetTextureTriangleArea() throws Exception {
		CoHDS hds = new CoHDS();
		HalfEdgeUtils.addTetrahedron(hds);
		CoFace f = hds.getFace(0);
		
		CoVertex A = f.getBoundaryEdge().getStartVertex();
		CoVertex B = f.getBoundaryEdge().getNextEdge().getStartVertex();
		CoVertex C = f.getBoundaryEdge().getPreviousEdge().getStartVertex();
		
		A.T = new double[] {0,0,0,1};
		B.T = new double[] {1,0,0,1};
		C.T = new double[] {0,1,0,1};
		
		double r = ConvergenceUtility.getTextureTriangleArea(f);
		Assert.assertEquals(0.5, r, 1E-10);
		
		B.T = new double[] {0.5,0.5,0,1};
		r = ConvergenceUtility.getTextureTriangleArea(f);
		Assert.assertEquals(0.25, r, 1E-10);
	}
	
	
	@Test
	public void testScaleInvariantCircumCircleRadius() throws Exception {
		CoHDS hds = new CoHDS();
		CoVertex v1 = hds.addNewVertex();
		CoVertex v2 = hds.addNewVertex();
		CoVertex v3 = hds.addNewVertex();
		CoVertex v4 = hds.addNewVertex();
		CoFace f1 = constructFaceByVertices(hds, v1, v2, v3).getLeftFace();
		CoFace f2 = constructFaceByVertices(hds, v1, v3, v4).getLeftFace();
		v1.T = new double[]{0,0,0,1};
		v2.T = new double[]{1,0,0,1};
		v3.T = new double[]{0,1,0,1};
		v4.T = new double[]{-1,0,0,1};
		Assert.assertEquals(0.5, getTextureTriangleArea(f1), 1E-10);
		Assert.assertEquals(0.5, getTextureTriangleArea(f2), 1E-10);
		double[] data1 = getMaxMeanSumScaleInvariantCircumRadius(hds);
		Assert.assertEquals(sin(Math.PI/4), data1[0], 1E-10);
		Assert.assertEquals(sin(Math.PI/4), data1[1], 1E-10);
		Assert.assertEquals(2.0*sin(Math.PI/4), data1[2], 1E-10);
		// scale down by 2.0
		v1.T[3] = 2.0;
		v2.T[3] = 2.0;
		v3.T[3] = 2.0;
		v4.T[3] = 2.0;
		Assert.assertEquals(0.125, getTextureTriangleArea(f1), 1E-10);
		Assert.assertEquals(0.125, getTextureTriangleArea(f2), 1E-10);
		double[] data2 = getMaxMeanSumScaleInvariantCircumRadius(hds);
		Assert.assertEquals(data1[0], data2[0], 1E-10);
		Assert.assertEquals(data1[1], data2[1], 1E-10);
		Assert.assertEquals(data1[2], data2[2], 1E-10);
	}
	
}
