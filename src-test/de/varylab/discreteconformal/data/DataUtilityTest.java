package de.varylab.discreteconformal.data;

import java.io.ByteArrayOutputStream;
import java.util.LinkedList;
import java.util.List;

import org.custommonkey.xmlunit.XMLAssert;
import org.junit.Test;

import de.jreality.math.Matrix;
import de.varylab.conformallab.data.DataIO;
import de.varylab.conformallab.data.DataUtility;
import de.varylab.conformallab.data.types.UniformizationData;
import de.varylab.discreteconformal.math.RnBig;
import de.varylab.discreteconformal.uniformization.FundamentalEdge;
import de.varylab.discreteconformal.uniformization.FundamentalPolygon;
import de.varylab.discreteconformal.uniformization.FundamentalVertex;


public class DataUtilityTest {

	@Test
	public void testToFuchsianData() throws Exception {
		List<FundamentalEdge> edges = new LinkedList<>();
		FundamentalVertex v = new FundamentalVertex(0);
		FundamentalEdge e0 = new FundamentalEdge(0);
		FundamentalEdge e1 = new FundamentalEdge(1);
		FundamentalEdge e2 = new FundamentalEdge(2);
		FundamentalEdge e3 = new FundamentalEdge(3);
		edges.add(e0);
		edges.add(e1);
		edges.add(e2);
		edges.add(e3);
		e0.end = v;
		e1.end = v;
		e2.end = v;
		e3.end = v;
		e0.motion = new Matrix();
		e1.motion = new Matrix();
		e2.motion = new Matrix();
		e3.motion = new Matrix();
		e0.motionBig = RnBig.toBig(null, e0.motion.getArray());
		e1.motionBig = RnBig.toBig(null, e1.motion.getArray());
		e2.motionBig = RnBig.toBig(null, e2.motion.getArray());
		e3.motionBig = RnBig.toBig(null, e3.motion.getArray());
		e0.nextEdge = e1; 
		e1.nextEdge = e2; 
		e2.nextEdge = e3; 
		e3.nextEdge = e0;
		e0.partner = e2;
		e1.partner = e3;
		e2.partner = e0;
		e3.partner = e1;
		e0.prevEdge = e3;
		e1.prevEdge = e0;
		e2.prevEdge = e1;
		e3.prevEdge = e2;
		e0.sourceEdgeCount = 1;
		e1.sourceEdgeCount = 1;
		e2.sourceEdgeCount = 1;
		e3.sourceEdgeCount = 1;
		e0.start = v;
		e1.start = v;
		e2.start = v;
		e3.start = v;
		e0.startPosition = RnBig.toBig(null, new double[]{1.0,0,0,1}); 
		e1.startPosition = RnBig.toBig(null, new double[]{0,2.0,0,1}); 
		e2.startPosition = RnBig.toBig(null, new double[]{0,0,3.0,1}); 
		e3.startPosition = RnBig.toBig(null, new double[]{1.0,0,0,4.0});
		
		FundamentalPolygon P = new FundamentalPolygon(edges);
		UniformizationData data = DataUtility.toUniformizationData("Fuchsian Group", P);
		ByteArrayOutputStream out = new ByteArrayOutputStream();
		DataIO.writeConformalData(data, out);
		String xml = new String(out.toByteArray());
		XMLAssert.assertXMLEqual(
			"<UniformizationData xmlns='http://www.varylab.com/conformallab/types' name='Fuchsian Group'>\n" + 
			"    <UniformizingGroup>\n" + 
			"        <IsometryPSL3R m11='1.0' m12='0.0' m13='0.0' m14='0.0' m21='0.0' m22='1.0' m23='0.0' m24='0.0' m31='0.0' m32='0.0' m33='1.0' m34='0.0' m41='0.0' m42='0.0' m43='0.0' m44='1.0'/>\n" + 
			"        <IsometryPSL3R m11='1.0' m12='0.0' m13='0.0' m14='0.0' m21='0.0' m22='1.0' m23='0.0' m24='0.0' m31='0.0' m32='0.0' m33='1.0' m34='0.0' m41='0.0' m42='0.0' m43='0.0' m44='1.0'/>\n" + 
			"        <IsometryPSL3R m11='1.0' m12='0.0' m13='0.0' m14='0.0' m21='0.0' m22='1.0' m23='0.0' m24='0.0' m31='0.0' m32='0.0' m33='1.0' m34='0.0' m41='0.0' m42='0.0' m43='0.0' m44='1.0'/>\n" + 
			"        <IsometryPSL3R m11='1.0' m12='0.0' m13='0.0' m14='0.0' m21='0.0' m22='1.0' m23='0.0' m24='0.0' m31='0.0' m32='0.0' m33='1.0' m34='0.0' m41='0.0' m42='0.0' m43='0.0' m44='1.0'/>\n" + 
			"    </UniformizingGroup>\n" + 
			"    <FundamentalPolygon>\n" + 
			"        <FundamentalEdge identifiedEdge='2'>\n" + 
			"            <StartPosition re='1.0' im='0.0'/>\n" + 
			"        </FundamentalEdge>\n" + 
			"        <FundamentalEdge identifiedEdge='3'>\n" + 
			"            <StartPosition re='0.0' im='2.0'/>\n" + 
			"        </FundamentalEdge>\n" + 
			"        <FundamentalEdge identifiedEdge='0'>\n" + 
			"            <StartPosition re='0.0' im='0.0'/>\n" + 
			"        </FundamentalEdge>\n" + 
			"        <FundamentalEdge identifiedEdge='1'>\n" + 
			"            <StartPosition re='0.25' im='0.0'/>\n" + 
			"        </FundamentalEdge>\n" + 
			"    </FundamentalPolygon>\n" + 
			"</UniformizationData>" , xml
		);
	}

}
