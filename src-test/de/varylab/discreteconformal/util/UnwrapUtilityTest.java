package de.varylab.discreteconformal.util;

import static java.lang.Math.PI;

import org.junit.Assert;
import org.junit.Test;

import de.jtem.halfedge.util.HalfEdgeUtils;
import de.varylab.discreteconformal.ConformalAdapterSet;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;

public class UnwrapUtilityTest {

	@Test
	public void testGetAngleReturnsPI() throws Exception {
		CoHDS hds = new CoHDS();
		CoVertex v0 = hds.addNewVertex();
		CoVertex v1 = hds.addNewVertex();
		CoVertex v2 = hds.addNewVertex();
		v0.P[0] = -1;
		v2.P[0] = 1;
		HalfEdgeUtils.constructFaceByVertices(hds, v0, v1, v2);
		CoEdge e = HalfEdgeUtils.findEdgeBetweenVertices(v2, v1);
		double angle = UnwrapUtility.getAngle(e, new ConformalAdapterSet());
		Assert.assertEquals(PI, angle, 1E-15);
	}
	
}
