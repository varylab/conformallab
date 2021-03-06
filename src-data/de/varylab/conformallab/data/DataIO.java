package de.varylab.conformallab.data;

import static javax.xml.XMLConstants.W3C_XML_SCHEMA_NS_URI;
import static javax.xml.bind.Marshaller.JAXB_FORMATTED_OUTPUT;

import java.io.InputStream;
import java.io.OutputStream;
import java.util.logging.Logger;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBElement;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;
import javax.xml.transform.Source;
import javax.xml.transform.stream.StreamSource;
import javax.xml.validation.Schema;
import javax.xml.validation.SchemaFactory;

import de.varylab.conformallab.data.types.ConformalData;
import de.varylab.conformallab.data.types.ConformalDataList;
import de.varylab.conformallab.data.types.DiscreteEmbedding;
import de.varylab.conformallab.data.types.DiscreteMap;
import de.varylab.conformallab.data.types.DiscreteMetric;
import de.varylab.conformallab.data.types.HalfedgeEmbedding;
import de.varylab.conformallab.data.types.HalfedgeMap;
import de.varylab.conformallab.data.types.HyperEllipticAlgebraicCurve;
import de.varylab.conformallab.data.types.ObjectFactory;
import de.varylab.conformallab.data.types.SchottkyData;
import de.varylab.conformallab.data.types.UniformizationData;

public class DataIO {

	private static Logger
		log = Logger.getLogger(DataIO.class.getName());
	private static Schema
		typesSchema = null;
	private static JAXBContext
		typesContex = null;
	private static Unmarshaller
		typesUnmarshaller = null;
	private static Marshaller
		typesMarshaller = null;
	
	static {
		try {
			SchemaFactory factory = SchemaFactory.newInstance(W3C_XML_SCHEMA_NS_URI);
			Source typesSource = new StreamSource(DataIO.class.getResourceAsStream("types.xsd"));
			typesSchema = factory.newSchema(typesSource);
			typesContex = JAXBContext.newInstance("de.varylab.conformallab.data.types");
			typesUnmarshaller = typesContex.createUnmarshaller();
			typesUnmarshaller.setSchema(typesSchema);
			typesMarshaller = typesContex.createMarshaller();
			typesMarshaller.setSchema(typesSchema);
			typesMarshaller.setProperty(JAXB_FORMATTED_OUTPUT, true);
		} catch (Exception e) {
			log.severe(e.getMessage());
		}
	}

	public static <T extends ConformalData> T readConformalData(Class<T> dataClass, InputStream in) throws JAXBException {
		StreamSource inSource = new StreamSource(in);
		JAXBElement<T> e = typesUnmarshaller.unmarshal(inSource, dataClass);
		return e.getValue();
	}
	@SuppressWarnings("unchecked")
	public static ConformalData readConformalData(InputStream in) throws JAXBException {
		StreamSource inSource = new StreamSource(in);
		JAXBElement<? extends ConformalData> e = (JAXBElement<? extends ConformalData>)typesUnmarshaller.unmarshal(inSource);
		return e.getValue();
	}
	public static Object read(InputStream in) throws JAXBException {
		StreamSource inSource = new StreamSource(in);
		JAXBElement<?> e = (JAXBElement<?>)typesUnmarshaller.unmarshal(inSource);
		return e.getValue();
	}	
	public static <T extends ConformalData> void writeConformalData(T data, OutputStream out) throws JAXBException {
		ObjectFactory of = new ObjectFactory();
		JAXBElement<?> e = null;
		if (data instanceof UniformizationData) {
			e = of.createUniformizationData((UniformizationData)data);
		} else
		if (data instanceof HyperEllipticAlgebraicCurve) {
			e = of.createHyperEllipticAlgebraicCurve((HyperEllipticAlgebraicCurve)data);
		} else 
		if (data instanceof SchottkyData) {
			e = of.createSchottkyData((SchottkyData)data);
		} else
		if (data instanceof DiscreteMetric) {
			e = of.createDiscreteMetric((DiscreteMetric)data);
		} else
		if (data instanceof DiscreteEmbedding) {
			// convert to HalfedgeEmbedding to avoid deprecated DiscreteEmbedding
			DiscreteEmbedding de = (DiscreteEmbedding)data;
			HalfedgeEmbedding he = DataUtility.toHalfedgeEmbedding(de);
			e = of.createHalfedgeEmbedding(he);
		} else
		if (data instanceof DiscreteMap) {
			// convert to HalfedgeMap to avoid deprecated DiscreteMap
			DiscreteMap dm = (DiscreteMap)data;
			HalfedgeMap hm = DataUtility.toHalfedgeMap(dm);
			e = of.createHalfedgeMap(hm);
		} else
		if (data instanceof HalfedgeEmbedding) {
			e = of.createHalfedgeEmbedding((HalfedgeEmbedding)data);
		} else
		if (data instanceof HalfedgeMap) {
			e = of.createHalfedgeMap((HalfedgeMap)data);
		}
		typesMarshaller.marshal(e, out);
	}
	
	public static ConformalDataList readConformalDataList(InputStream in) throws JAXBException {
		StreamSource inSource = new StreamSource(in);
		JAXBElement<ConformalDataList> e = typesUnmarshaller.unmarshal(inSource, ConformalDataList.class);
		return e.getValue();
	}
	public static void writeConformalDataList(ConformalDataList data, OutputStream out) throws JAXBException {
		ObjectFactory of = new ObjectFactory();
		ConformalDataList writeData = of.createConformalDataList();
		// convert to avoid writing DiscreteEmbedding or DiscreteMap
		for (ConformalData d : data.getData()) {
			if (d instanceof DiscreteEmbedding) {
				DiscreteEmbedding de = (DiscreteEmbedding)d;
				HalfedgeEmbedding he = DataUtility.toHalfedgeEmbedding(de);
				writeData.getData().add(he);
			} else
			if (d instanceof DiscreteMap) {
				DiscreteMap dm = (DiscreteMap)d;
				HalfedgeMap hm = DataUtility.toHalfedgeMap(dm);
				writeData.getData().add(hm);
			} else {
				writeData.getData().add(d);
			}
		}
		JAXBElement<ConformalDataList> e = of.createConformalDataList(writeData);
		typesMarshaller.marshal(e, out);
	}
	
}
