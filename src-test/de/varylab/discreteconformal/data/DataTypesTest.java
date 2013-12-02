package de.varylab.discreteconformal.data;

import java.io.ByteArrayOutputStream;
import java.io.InputStream;

import org.junit.Assert;
import org.junit.Test;

import de.varylab.conformallab.data.DataFactory;
import de.varylab.conformallab.data.types.DiscreteEmbedding;
import de.varylab.conformallab.data.types.DiscreteMetric;
import de.varylab.conformallab.data.types.EmbeddedTriangle;
import de.varylab.conformallab.data.types.EmbeddedVertex;
import de.varylab.conformallab.data.types.MetricEdge;
import de.varylab.conformallab.data.types.MetricTriangle;

public class DataTypesTest {

	@Test
	public void testReadDiscreteEmbedding() throws Exception {
		InputStream in = getClass().getResourceAsStream("discreteembedding.xml");
		DiscreteEmbedding de = DataFactory.readDiscreteEmbedding(in);
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
		DataFactory.writeDiscreteEmbedding(de, out);
		String xml = new String(out.toByteArray());
		Assert.assertEquals(
			"<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n" + 
			"<DiscreteEmbedding xmlns=\"http://www.varylab.com/conformallab/types\">\n" + 
			"    <EmbeddedVertex index=\"0\" w=\"0.0\" z=\"0.0\" y=\"0.0\" x=\"3.141592653589793\"/>\n" + 
			"    <EmbeddedVertex index=\"1\" w=\"0.0\" z=\"0.0\" y=\"2.718281828459045\" x=\"0.0\"/>\n" + 
			"    <EmbeddedVertex index=\"2\" w=\"3.141592653589793\" z=\"0.0\" y=\"0.0\" x=\"0.0\"/>\n" + 
			"    <EmbeddedTriangle vertex3=\"2\" vertex2=\"1\" vertex1=\"0\"/>\n" + 
			"</DiscreteEmbedding>\n", xml
		);
	}
	
	@Test
	public void testReadDiscreteMetric() throws Exception {
		InputStream in = getClass().getResourceAsStream("discretemetric.xml");
		DiscreteMetric dm = DataFactory.readDiscreteMetric(in);
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
		DataFactory.writeDiscreteMetric(dm, out);
		String xml = new String(out.toByteArray());
		Assert.assertEquals(
			"<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n" + 
			"<DiscreteMetric xmlns=\"http://www.varylab.com/conformallab/types\">\n" + 
			"    <MetricEdge index=\"0\" length=\"3.141592653589793\"/>\n" + 
			"    <MetricEdge index=\"1\" length=\"2.718281828459045\"/>\n" + 
			"    <MetricEdge index=\"2\" length=\"2.0\"/>\n" + 
			"    <MetricTriangle edge3=\"2\" edge2=\"1\" edge1=\"0\"/>\n" + 
			"</DiscreteMetric>\n", xml
		);
	}
	
}
