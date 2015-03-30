package de.varylab.discreteconformal.functional;

import java.io.InputStream;

import org.junit.Test;

import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Position;
import de.varylab.conformallab.data.DataIO;
import de.varylab.conformallab.data.DataUtility;
import de.varylab.conformallab.data.types.HalfedgeEmbedding;
import de.varylab.discreteconformal.ConformalAdapterSet;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;

public class HyperIdealHyperEllipticCurveTests {

	@Test
	public void testCreateDataOnHyperEllipticCurveLawson() throws Exception {
		InputStream in = getClass().getResourceAsStream("lawson_curve_source.xml");
		HalfedgeEmbedding he = DataIO.readConformalData(HalfedgeEmbedding.class, in);
		CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo = new CuttingInfo<>();
		CoHDS hds = new CoHDS();
		AdapterSet a = new ConformalAdapterSet();
		DataUtility.toHalfedge(he, a, Position.class, hds, cutInfo);
		
	}
	
	
}
