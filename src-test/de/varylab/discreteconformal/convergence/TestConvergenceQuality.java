package de.varylab.discreteconformal.convergence;

import org.junit.Test;

import de.jreality.junitutils.Assert;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;

public class TestConvergenceQuality {

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
		
		double r = ConvergenceQuality.getTextureCircumCircleRadius(f);
		Assert.assertEquals(Math.sqrt(2)/2, r, 1E-10);
		
		B.T = new double[] {0.5,0.5,0,1};
		r = ConvergenceQuality.getTextureCircumCircleRadius(f);
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
		
		double r = ConvergenceQuality.getTextureTriangleArea(f);
		Assert.assertEquals(0.5, r, 1E-10);
		
		B.T = new double[] {0.5,0.5,0,1};
		r = ConvergenceQuality.getTextureTriangleArea(f);
		Assert.assertEquals(0.25, r, 1E-10);
	}
	
}
