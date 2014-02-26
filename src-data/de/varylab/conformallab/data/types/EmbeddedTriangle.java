//
// Diese Datei wurde mit der JavaTM Architecture for XML Binding(JAXB) Reference Implementation, v2.2.6 generiert 
// Siehe <a href="http://java.sun.com/xml/jaxb">http://java.sun.com/xml/jaxb</a> 
// Änderungen an dieser Datei gehen bei einer Neukompilierung des Quellschemas verloren. 
// Generiert: 2014.02.26 um 05:27:05 PM CET 
//


package de.varylab.conformallab.data.types;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlType;


/**
 * <p>Java-Klasse für EmbeddedTriangle complex type.
 * 
 * <p>Das folgende Schemafragment gibt den erwarteten Content an, der in dieser Klasse enthalten ist.
 * 
 * <pre>
 * &lt;complexType name="EmbeddedTriangle">
 *   &lt;complexContent>
 *     &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
 *       &lt;attribute name="vertex1" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
 *       &lt;attribute name="vertex2" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
 *       &lt;attribute name="vertex3" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
 *     &lt;/restriction>
 *   &lt;/complexContent>
 * &lt;/complexType>
 * </pre>
 * 
 * 
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "EmbeddedTriangle")
public class EmbeddedTriangle {

    @XmlAttribute(name = "vertex1", required = true)
    protected int vertex1;
    @XmlAttribute(name = "vertex2", required = true)
    protected int vertex2;
    @XmlAttribute(name = "vertex3", required = true)
    protected int vertex3;

    /**
     * Ruft den Wert der vertex1-Eigenschaft ab.
     * 
     */
    public int getVertex1() {
        return vertex1;
    }

    /**
     * Legt den Wert der vertex1-Eigenschaft fest.
     * 
     */
    public void setVertex1(int value) {
        this.vertex1 = value;
    }

    /**
     * Ruft den Wert der vertex2-Eigenschaft ab.
     * 
     */
    public int getVertex2() {
        return vertex2;
    }

    /**
     * Legt den Wert der vertex2-Eigenschaft fest.
     * 
     */
    public void setVertex2(int value) {
        this.vertex2 = value;
    }

    /**
     * Ruft den Wert der vertex3-Eigenschaft ab.
     * 
     */
    public int getVertex3() {
        return vertex3;
    }

    /**
     * Legt den Wert der vertex3-Eigenschaft fest.
     * 
     */
    public void setVertex3(int value) {
        this.vertex3 = value;
    }

}
