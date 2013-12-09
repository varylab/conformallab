//
// Diese Datei wurde mit der JavaTM Architecture for XML Binding(JAXB) Reference Implementation, v2.2.6 generiert 
// Siehe <a href="http://java.sun.com/xml/jaxb">http://java.sun.com/xml/jaxb</a> 
// Änderungen an dieser Datei gehen bei einer Neukompilierung des Quellschemas verloren. 
// Generiert: 2013.12.09 um 07:42:44 PM CET 
//


package de.varylab.conformallab.data.surface;

import javax.xml.bind.JAXBElement;
import javax.xml.bind.annotation.XmlElementDecl;
import javax.xml.bind.annotation.XmlRegistry;
import javax.xml.namespace.QName;


/**
 * This object contains factory methods for each 
 * Java content interface and Java element interface 
 * generated in the de.varylab.conformallab.data.surface package. 
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

    private final static QName _Surface_QNAME = new QName("http://www.varylab.com/conformallab/surface", "Surface");

    /**
     * Create a new ObjectFactory that can be used to create new instances of schema derived classes for package: de.varylab.conformallab.data.surface
     * 
     */
    public ObjectFactory() {
    }

    /**
     * Create an instance of {@link Surface }
     * 
     */
    public Surface createSurface() {
        return new Surface();
    }

    /**
     * Create an instance of {@link JAXBElement }{@code <}{@link Surface }{@code >}}
     * 
     */
    @XmlElementDecl(namespace = "http://www.varylab.com/conformallab/surface", name = "Surface")
    public JAXBElement<Surface> createSurface(Surface value) {
        return new JAXBElement<Surface>(_Surface_QNAME, Surface.class, null, value);
    }

}
