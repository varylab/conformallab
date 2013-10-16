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

import de.varylab.conformallab.data.types.ObjectFactory;
import de.varylab.conformallab.data.types.SchottkyData;

public class DataFactory {

	private static Logger
		log = Logger.getLogger(DataFactory.class.getName());
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
			Source typesSource = new StreamSource(DataFactory.class.getResourceAsStream("types.xsd"));
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
	
	public static SchottkyData loadSchottkyData(InputStream in) throws JAXBException {
		StreamSource inSource = new StreamSource(in);
		JAXBElement<SchottkyData> e = typesUnmarshaller.unmarshal(inSource, SchottkyData.class);
		return e.getValue();
	}
	
	public static void writeSchottkyData(SchottkyData data, OutputStream out) throws JAXBException {
		ObjectFactory of = new ObjectFactory();
		JAXBElement<SchottkyData> e = of.createSchottkyData(data);
		typesMarshaller.marshal(e, out);
	}
	
	
}
