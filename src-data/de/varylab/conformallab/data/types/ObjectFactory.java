//
// Diese Datei wurde mit der JavaTM Architecture for XML Binding(JAXB) Reference Implementation, v2.2.6 generiert 
// Siehe <a href="http://java.sun.com/xml/jaxb">http://java.sun.com/xml/jaxb</a> 
// Änderungen an dieser Datei gehen bei einer Neukompilierung des Quellschemas verloren. 
// Generiert: 2013.10.16 um 05:25:04 PM CEST 
//


package de.varylab.conformallab.data.types;

import javax.xml.bind.JAXBElement;
import javax.xml.bind.annotation.XmlElementDecl;
import javax.xml.bind.annotation.XmlRegistry;
import javax.xml.namespace.QName;


/**
 * This object contains factory methods for each 
 * Java content interface and Java element interface 
 * generated in the de.varylab.conformallab.data.types package. 
 * <p>An ObjectFactory allows you to programatically 
 * construct new instances of the Java representation 
 * for XML content. The Java representation of XML 
 * content can consist of schema derived interfaces 
 * and classes representing the binding of schema 
 * type definitions, element declarations and model 
 * groups.  Factory methods for each of these are 
 * provided in this class.
 * 
 */
@XmlRegistry
public class ObjectFactory {

    private final static QName _SchottkyData_QNAME = new QName("http://www.varylab.com/conformallab/types", "SchottkyData");

    /**
     * Create a new ObjectFactory that can be used to create new instances of schema derived classes for package: de.varylab.conformallab.data.types
     * 
     */
    public ObjectFactory() {
    }

    /**
     * Create an instance of {@link SchottkyData }
     * 
     */
    public SchottkyData createSchottkyData() {
        return new SchottkyData();
    }

    /**
     * Create an instance of {@link Moebius }
     * 
     */
    public Moebius createMoebius() {
        return new Moebius();
    }

    /**
     * Create an instance of {@link Circle }
     * 
     */
    public Circle createCircle() {
        return new Circle();
    }

    /**
     * Create an instance of {@link SchottkyGenerator }
     * 
     */
    public SchottkyGenerator createSchottkyGenerator() {
        return new SchottkyGenerator();
    }

    /**
     * Create an instance of {@link Complex }
     * 
     */
    public Complex createComplex() {
        return new Complex();
    }

    /**
     * Create an instance of {@link JAXBElement }{@code <}{@link SchottkyData }{@code >}}
     * 
     */
    @XmlElementDecl(namespace = "http://www.varylab.com/conformallab/types", name = "SchottkyData")
    public JAXBElement<SchottkyData> createSchottkyData(SchottkyData value) {
        return new JAXBElement<SchottkyData>(_SchottkyData_QNAME, SchottkyData.class, null, value);
    }

}
