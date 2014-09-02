//
// Diese Datei wurde mit der JavaTM Architecture for XML Binding(JAXB) Reference Implementation, v2.2.6 generiert 
// Siehe <a href="http://java.sun.com/xml/jaxb">http://java.sun.com/xml/jaxb</a> 
// Änderungen an dieser Datei gehen bei einer Neukompilierung des Quellschemas verloren. 
// Generiert: 2014.09.02 um 04:53:14 PM CEST 
//


package de.varylab.conformallab.data.types;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlType;


/**
 * <p>Java-Klasse für EdgeIdentification complex type.
 * 
 * <p>Das folgende Schemafragment gibt den erwarteten Content an, der in dieser Klasse enthalten ist.
 * 
 * <pre>
 * &lt;complexType name="EdgeIdentification">
 *   &lt;complexContent>
 *     &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
 *       &lt;attribute name="edge1" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
 *       &lt;attribute name="edge2" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
 *       &lt;attribute name="edge3" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
 *       &lt;attribute name="edge4" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
 *     &lt;/restriction>
 *   &lt;/complexContent>
 * &lt;/complexType>
 * </pre>
 * 
 * 
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "EdgeIdentification")
public class EdgeIdentification {

    @XmlAttribute(name = "edge1", required = true)
    protected int edge1;
    @XmlAttribute(name = "edge2", required = true)
    protected int edge2;
    @XmlAttribute(name = "edge3", required = true)
    protected int edge3;
    @XmlAttribute(name = "edge4", required = true)
    protected int edge4;

    /**
     * Ruft den Wert der edge1-Eigenschaft ab.
     * 
     */
    public int getEdge1() {
        return edge1;
    }

    /**
     * Legt den Wert der edge1-Eigenschaft fest.
     * 
     */
    public void setEdge1(int value) {
        this.edge1 = value;
    }

    /**
     * Ruft den Wert der edge2-Eigenschaft ab.
     * 
     */
    public int getEdge2() {
        return edge2;
    }

    /**
     * Legt den Wert der edge2-Eigenschaft fest.
     * 
     */
    public void setEdge2(int value) {
        this.edge2 = value;
    }

    /**
     * Ruft den Wert der edge3-Eigenschaft ab.
     * 
     */
    public int getEdge3() {
        return edge3;
    }

    /**
     * Legt den Wert der edge3-Eigenschaft fest.
     * 
     */
    public void setEdge3(int value) {
        this.edge3 = value;
    }

    /**
     * Ruft den Wert der edge4-Eigenschaft ab.
     * 
     */
    public int getEdge4() {
        return edge4;
    }

    /**
     * Legt den Wert der edge4-Eigenschaft fest.
     * 
     */
    public void setEdge4(int value) {
        this.edge4 = value;
    }

}
