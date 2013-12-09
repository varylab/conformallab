//
// Diese Datei wurde mit der JavaTM Architecture for XML Binding(JAXB) Reference Implementation, v2.2.6 generiert 
// Siehe <a href="http://java.sun.com/xml/jaxb">http://java.sun.com/xml/jaxb</a> 
// Änderungen an dieser Datei gehen bei einer Neukompilierung des Quellschemas verloren. 
// Generiert: 2013.12.09 um 07:42:44 PM CET 
//


package de.varylab.conformallab.data.types;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlType;


/**
 * <p>Java-Klasse für FundamentalEdge complex type.
 * 
 * <p>Das folgende Schemafragment gibt den erwarteten Content an, der in dieser Klasse enthalten ist.
 * 
 * <pre>
 * &lt;complexType name="FundamentalEdge">
 *   &lt;complexContent>
 *     &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
 *       &lt;attribute name="vertex1" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
 *       &lt;attribute name="vertex2" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
 *       &lt;attribute name="identificationIndex" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
 *     &lt;/restriction>
 *   &lt;/complexContent>
 * &lt;/complexType>
 * </pre>
 * 
 * 
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "FundamentalEdge")
public class FundamentalEdge {

    @XmlAttribute(name = "vertex1", required = true)
    protected int vertex1;
    @XmlAttribute(name = "vertex2", required = true)
    protected int vertex2;
    @XmlAttribute(name = "identificationIndex", required = true)
    protected int identificationIndex;

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
     * Ruft den Wert der identificationIndex-Eigenschaft ab.
     * 
     */
    public int getIdentificationIndex() {
        return identificationIndex;
    }

    /**
     * Legt den Wert der identificationIndex-Eigenschaft fest.
     * 
     */
    public void setIdentificationIndex(int value) {
        this.identificationIndex = value;
    }

}
