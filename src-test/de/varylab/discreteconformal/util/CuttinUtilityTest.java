package de.varylab.discreteconformal.util;

import org.junit.Assert;
import org.junit.Test;

import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.adapter.CoTexturePositionAdapter;

public class CuttinUtilityTest {

	@Test
	public void testIsInConvexTextureFace_False() throws Exception {
		double[] p1 = {0.7488102998904661, 0.06293998610761144, 0.0, 1.0};
		double[] p2 = {0.7487811940754379, 0.06289451051246124, 0.0, 1.0};
		double[] p3 = {0.7487254625255592, 0.06291429499873116, 0.0, 1.0};
		double[] pp = {0.44661534423161037, 2.2808373704822393E-4, 0.0, 1.0};
		CoHDS hds = new CoHDS();
		CoVertex v1 = hds.addNewVertex();
		CoVertex v2 = hds.addNewVertex();
		CoVertex v3 = hds.addNewVertex();
		v1.T = p1;
		v2.T = p2;
		v3.T = p3;
		CoFace f = HalfEdgeUtils.constructFaceByVertices(hds, v1, v2, v3).getLeftFace();
		AdapterSet a = AdapterSet.createGenericAdapters();
		a.add(new CoTexturePositionAdapter());
		boolean isIn = CuttingUtility.isInConvexTextureFace(pp, f, a);
		Assert.assertFalse(isIn);
	}
	
	@Test
	public void testIsInConvexTextureFace_True() throws Exception {
		double[] p1 = {0.0, 1.0E-8, 0.0, 1.0};
		double[] p2 = {0.0, 1.0E-8, 1E-8, 0.0, 1.0};
		double[] p3 = {0.0, 0.0, 0.0, 1.0};
		double[] pp = {0.0, 0.5E-8, 0.25E-8, 1.0};
		CoHDS hds = new CoHDS();
		CoVertex v1 = hds.addNewVertex();
		CoVertex v2 = hds.addNewVertex();
		CoVertex v3 = hds.addNewVertex();
		v1.T = p1;
		v2.T = p2;
		v3.T = p3;
		CoFace f = HalfEdgeUtils.constructFaceByVertices(hds, v1, v2, v3).getLeftFace();
		AdapterSet a = AdapterSet.createGenericAdapters();
		a.add(new CoTexturePositionAdapter());
		boolean isIn = CuttingUtility.isInConvexTextureFace(pp, f, a);
		Assert.assertTrue(isIn);
	}

	
}
