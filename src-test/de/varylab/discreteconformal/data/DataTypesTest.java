package de.varylab.discreteconformal.data;

import java.io.ByteArrayOutputStream;
import java.io.InputStream;

import org.custommonkey.xmlunit.XMLAssert;
import org.junit.Assert;
import org.junit.Test;

import de.varylab.conformallab.data.DataIO;
import de.varylab.conformallab.data.types.Complex;
import de.varylab.conformallab.data.types.DiscreteEmbedding;
import de.varylab.conformallab.data.types.DiscreteMetric;
import de.varylab.conformallab.data.types.EmbeddedTriangle;
import de.varylab.conformallab.data.types.EmbeddedVertex;
import de.varylab.conformallab.data.types.FundamentalEdge;
import de.varylab.conformallab.data.types.FundamentalPolygon;
import de.varylab.conformallab.data.types.FundamentalVertex;
import de.varylab.conformallab.data.types.IsometryPSL2R;
import de.varylab.conformallab.data.types.MetricEdge;
import de.varylab.conformallab.data.types.MetricTriangle;
import de.varylab.conformallab.data.types.ObjectFactory;
import de.varylab.conformallab.data.types.UniformizationData;
import de.varylab.conformallab.data.types.UniformizingGroup;
import de.varylab.conformallab.data.types.VertexIdentification;

public class DataTypesTest {

	@Test
	public void testReadFuchsianData() throws Exception {
		InputStream in = getClass().getResourceAsStream("wente_uniformization.xml");
		UniformizationData data = DataIO.readConformalData(UniformizationData.class, in);
		UniformizingGroup group = data.getUniformizingGroup();
		IsometryPSL2R T0 = group.getGenerators().get(0);
		Assert.assertEquals(1.0, T0.getM11(), 0.0);
		Assert.assertEquals(1.0, T0.getM22(), 0.0);
		Assert.assertEquals(1.0, T0.getM33(), 0.0);
		Assert.assertEquals(0.0, T0.getM12(), 0.0);
		Assert.assertEquals(1.414213562373094, T0.getM13(), 1E-15);
		FundamentalEdge e0 = data.getFundamentalPolygon().getEdges().get(0);
		FundamentalEdge e1 = data.getFundamentalPolygon().getEdges().get(1);
		FundamentalEdge e2 = data.getFundamentalPolygon().getEdges().get(2);
		FundamentalEdge e3 = data.getFundamentalPolygon().getEdges().get(3);
		Assert.assertEquals(2, e0.getIdentifiedEdge());
		Assert.assertEquals(3, e1.getIdentifiedEdge());
		Assert.assertEquals(0, e2.getIdentifiedEdge());
		Assert.assertEquals(1, e3.getIdentifiedEdge());
		Assert.assertEquals(2.4494897427831788, e0.getStartPosition().getIm(), 1E-15);
		Assert.assertEquals(2.8284271247461885, e1.getStartPosition().getRe(), 1E-15);
		Assert.assertEquals(2.4494897427831797, e1.getStartPosition().getIm(), 1E-15);
		Assert.assertEquals(1.4142135623730951, e2.getStartPosition().getRe(), 1E-15);
		Assert.assertEquals(-1.4142135623730951, e3.getStartPosition().getRe(), 1E-15);
		Assert.assertEquals(0, e0.getStartVertex());
		FundamentalVertex v = data.getFundamentalPolygon().getVertices().get(0);
		Assert.assertEquals(0, v.getIndex());
	}
	
	@Test
	public void testWriteFuchsianData() throws Exception {
		ObjectFactory of = new ObjectFactory();
		UniformizationData fd = of.createUniformizationData();
		UniformizingGroup G = of.createUniformizingGroup();
		IsometryPSL2R T0 = of.createIsometryPSL2R();
		IsometryPSL2R T1 = of.createIsometryPSL2R();
		IsometryPSL2R T2 = of.createIsometryPSL2R();
		IsometryPSL2R T3 = of.createIsometryPSL2R();
		T0.setM11(1.0);
		T0.setM11(2.0);
		T0.setM11(3.0);
		T0.setM11(4.0);
		G.getGenerators().add(T0);
		G.getGenerators().add(T1);
		G.getGenerators().add(T2);
		G.getGenerators().add(T3);
		FundamentalPolygon P = of.createFundamentalPolygon();
		for (int i = 0; i < 8; i++) {
			FundamentalEdge e = of.createFundamentalEdge();
			e.setIndex(i);
			e.setIdentifiedEdge(0);
			e.setNextEdge((i + 1) % 8);
			e.setPreviousEdge((i + 7) % 8);
			Complex s = new Complex();
			s.setRe(1.0);
			s.setIm(2.0);
			e.setStartPosition(s);
			e.setStartVertex(0);
			P.getEdges().add(e);
		}
		FundamentalVertex v = of.createFundamentalVertex();
		v.setIndex(0);
		P.getVertices().add(v);
		fd.setUniformizingGroup(G);
		fd.setFundamentalPolygon(P);
		fd.setName("Test Name");
		ByteArrayOutputStream out = new ByteArrayOutputStream();
		DataIO.writeConformalData(fd, out);
		String xml = new String(out.toByteArray());
		XMLAssert.assertXMLEqual(
			"<?xml version='1.0' encoding='UTF-8' standalone='yes'?>\n" + 
			"<UniformizationData xmlns='http://www.varylab.com/conformallab/types' name='Test Name'>\n" + 
			"    <UniformizingGroup>\n" + 
			"        <IsometryPSL2R m11='4.0'/>\n" + 
			"        <IsometryPSL2R/>\n" + 
			"        <IsometryPSL2R/>\n" + 
			"        <IsometryPSL2R/>\n" + 
			"    </UniformizingGroup>\n" + 
			"    <FundamentalPolygon>\n" + 
			"        <FundamentalVertex index='0'/>\n" + 
			"        <FundamentalEdge identifiedEdge='0' startVertex='0' index='0' nextEdge='1' previousEdge='7'>\n" + 
			"            <StartPosition re='1.0' im='2.0'/>\n" + 
			"        </FundamentalEdge>\n" + 
			"        <FundamentalEdge identifiedEdge='0' startVertex='0' index='1' nextEdge='2' previousEdge='0'>\n" + 
			"            <StartPosition re='1.0' im='2.0'/>\n" + 
			"        </FundamentalEdge>\n" + 
			"        <FundamentalEdge identifiedEdge='0' startVertex='0' index='2' nextEdge='3' previousEdge='1'>\n" + 
			"            <StartPosition re='1.0' im='2.0'/>\n" + 
			"        </FundamentalEdge>\n" + 
			"        <FundamentalEdge identifiedEdge='0' startVertex='0' index='3' nextEdge='4' previousEdge='2'>\n" + 
			"            <StartPosition re='1.0' im='2.0'/>\n" + 
			"        </FundamentalEdge>\n" + 
			"        <FundamentalEdge identifiedEdge='0' startVertex='0' index='4' nextEdge='5' previousEdge='3'>\n" + 
			"            <StartPosition re='1.0' im='2.0'/>\n" + 
			"        </FundamentalEdge>\n" + 
			"        <FundamentalEdge identifiedEdge='0' startVertex='0' index='5' nextEdge='6' previousEdge='4'>\n" + 
			"            <StartPosition re='1.0' im='2.0'/>\n" + 
			"        </FundamentalEdge>\n" + 
			"        <FundamentalEdge identifiedEdge='0' startVertex='0' index='6' nextEdge='7' previousEdge='5'>\n" + 
			"            <StartPosition re='1.0' im='2.0'/>\n" + 
			"        </FundamentalEdge>\n" + 
			"        <FundamentalEdge identifiedEdge='0' startVertex='0' index='7' nextEdge='0' previousEdge='6'>\n" + 
			"            <StartPosition re='1.0' im='2.0'/>\n" + 
			"        </FundamentalEdge>\n" + 
			"    </FundamentalPolygon>\n" + 
			"</UniformizationData>\n", xml
		);
	}
	
	
	@Test
	public void testReadDiscreteEmbedding() throws Exception {
		InputStream in = getClass().getResourceAsStream("discreteembedding.xml");
		DiscreteEmbedding de = DataIO.readConformalData(DiscreteEmbedding.class, in);
		Assert.assertEquals(1, de.getTriangles().size());
		EmbeddedTriangle t = de.getTriangles().get(0);
		Assert.assertEquals(0, t.getVertex1());
		Assert.assertEquals(1, t.getVertex2());
		Assert.assertEquals(2, t.getVertex3());
		Assert.assertEquals(3, de.getVertices().size());
		
		Assert.assertEquals(0.0, de.getVertices().get(0).getX(), 1E-20);
		Assert.assertEquals(0.0, de.getVertices().get(0).getY(), 1E-20);
		Assert.assertEquals(0.0, de.getVertices().get(0).getZ(), 1E-20);
		Assert.assertEquals(1.0, de.getVertices().get(0).getW(), 1E-20);
		
		Assert.assertEquals(1.0, de.getVertices().get(1).getX(), 1E-20);
		Assert.assertEquals(0.0, de.getVertices().get(1).getY(), 1E-20);
		Assert.assertEquals(0.0, de.getVertices().get(1).getZ(), 1E-20);
		Assert.assertEquals(1.0, de.getVertices().get(1).getW(), 1E-20);
		
		Assert.assertEquals(1.0, de.getVertices().get(2).getX(), 1E-20);
		Assert.assertEquals(1.0, de.getVertices().get(2).getY(), 1E-20);
		Assert.assertEquals(0.0, de.getVertices().get(2).getZ(), 1E-20);
		Assert.assertEquals(1.0, de.getVertices().get(2).getW(), 1E-20);
		
		Assert.assertEquals(1, de.getIdentifications().size());
		VertexIdentification id = de.getIdentifications().get(0);
		Assert.assertEquals(2, id.getVertices().size());
		Assert.assertEquals(0, (int)id.getVertices().get(0));
		Assert.assertEquals(1, (int)id.getVertices().get(1));
	}
	
	
	@Test
	public void testWriteDiscreteEmbedding() throws Exception {
		DiscreteEmbedding de = new DiscreteEmbedding();
		EmbeddedVertex v0 = new EmbeddedVertex();
		EmbeddedVertex v1 = new EmbeddedVertex();
		EmbeddedVertex v2 = new EmbeddedVertex();
		v0.setIndex(0);
		v1.setIndex(1);
		v2.setIndex(2);
		v0.setX(Math.PI);
		v1.setY(Math.E);
		v2.setW(Math.PI);
		de.getVertices().add(v0);
		de.getVertices().add(v1);
		de.getVertices().add(v2);
		EmbeddedTriangle t = new EmbeddedTriangle();
		t.setVertex1(0);
		t.setVertex2(1);
		t.setVertex3(2);
		de.getTriangles().add(t);
		ByteArrayOutputStream out = new ByteArrayOutputStream();
		DataIO.writeConformalData(de, out);
		String xml = new String(out.toByteArray());
		XMLAssert.assertXMLEqual(
			"<?xml version='1.0' encoding='UTF-8'?>" + 
			"<DiscreteEmbedding xmlns='http://www.varylab.com/conformallab/types'>\n" + 
			"    <EmbeddedVertex index='0' w='0.0' z='0.0' y='0.0' x='3.141592653589793'/>\n" + 
			"    <EmbeddedVertex index='1' w='0.0' z='0.0' y='2.718281828459045' x='0.0'/>\n" + 
			"    <EmbeddedVertex index='2' w='3.141592653589793' z='0.0' y='0.0' x='0.0'/>\n" + 
			"    <EmbeddedTriangle vertex3='2' vertex2='1' vertex1='0'/>\n" + 
			"</DiscreteEmbedding>", xml
		);
	}
	
	@Test
	public void testReadDiscreteMetric() throws Exception {
		InputStream in = getClass().getResourceAsStream("discretemetric.xml");
		DiscreteMetric dm = DataIO.readConformalData(DiscreteMetric.class, in);
		Assert.assertEquals(1, dm.getTriangles().size());
		Assert.assertEquals(3, dm.getEdges().size());
		MetricTriangle t = dm.getTriangles().get(0);
		Assert.assertEquals(0, t.getEdge1());
		Assert.assertEquals(1, t.getEdge2());
		Assert.assertEquals(2, t.getEdge3());
		Assert.assertEquals(Math.PI, dm.getEdges().get(0).getLength(), 1E-20);
		Assert.assertEquals(Math.E, dm.getEdges().get(1).getLength(), 1E-20);
		Assert.assertEquals(1.0, dm.getEdges().get(2).getLength(), 1E-20);
	}
	
	
	@Test
	public void testWriteDiscreteMetric() throws Exception {
		DiscreteMetric dm = new DiscreteMetric();
		MetricEdge e1 = new MetricEdge();
		MetricEdge e2 = new MetricEdge();
		MetricEdge e3 = new MetricEdge();
		e1.setIndex(0);
		e2.setIndex(1);
		e3.setIndex(2);
		e1.setLength(Math.PI);
		e2.setLength(Math.E);
		e3.setLength(2.0);
		dm.getEdges().add(e1);
		dm.getEdges().add(e2);
		dm.getEdges().add(e3);
		MetricTriangle t = new MetricTriangle();
		t.setEdge1(0);
		t.setEdge2(1);
		t.setEdge3(2);
		dm.getTriangles().add(t);
		
		ByteArrayOutputStream out = new ByteArrayOutputStream();
		DataIO.writeConformalData(dm, out);
		String xml = new String(out.toByteArray());
		XMLAssert.assertXMLEqual(
			"<?xml version='1.0' encoding='UTF-8'?>" + 
			"<DiscreteMetric xmlns='http://www.varylab.com/conformallab/types'>\n" + 
			"    <MetricEdge index='0' length='3.141592653589793'/>\n" + 
			"    <MetricEdge index='1' length='2.718281828459045'/>\n" + 
			"    <MetricEdge index='2' length='2.0'/>\n" + 
			"    <MetricTriangle edge3='2' edge2='1' edge1='0'/>\n" + 
			"</DiscreteMetric>", xml
		);
	}
	
}
