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
 * <p>Java-Klasse für HalfedgeEdge complex type.
 * 
 * <p>Das folgende Schemafragment gibt den erwarteten Content an, der in dieser Klasse enthalten ist.
 * 
 * <pre>
 * &lt;complexType name="HalfedgeEdge">
 *   &lt;complexContent>
 *     &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
 *       &lt;attribute name="left" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
 *       &lt;attribute name="target" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
 *       &lt;attribute name="next" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
 *       &lt;attribute name="opposite" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
 *       &lt;attribute name="index" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
 *     &lt;/restriction>
 *   &lt;/complexContent>
 * &lt;/complexType>
 * </pre>
 * 
 * 
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "HalfedgeEdge")
public class HalfedgeEdge {

    @XmlAttribute(name = "left", required = true)
    protected int left;
    @XmlAttribute(name = "target", required = true)
    protected int target;
    @XmlAttribute(name = "next", required = true)
    protected int next;
    @XmlAttribute(name = "opposite", required = true)
    protected int opposite;
    @XmlAttribute(name = "index", required = true)
    protected int index;

    /**
     * Ruft den Wert der left-Eigenschaft ab.
     * 
     */
    public int getLeft() {
        return left;
    }

    /**
     * Legt den Wert der left-Eigenschaft fest.
     * 
     */
    public void setLeft(int value) {
        this.left = value;
    }

    /**
     * Ruft den Wert der target-Eigenschaft ab.
     * 
     */
    public int getTarget() {
        return target;
    }

    /**
     * Legt den Wert der target-Eigenschaft fest.
     * 
     */
    public void setTarget(int value) {
        this.target = value;
    }

    /**
     * Ruft den Wert der next-Eigenschaft ab.
     * 
     */
    public int getNext() {
        return next;
    }

    /**
     * Legt den Wert der next-Eigenschaft fest.
     * 
     */
    public void setNext(int value) {
        this.next = value;
    }

    /**
     * Ruft den Wert der opposite-Eigenschaft ab.
     * 
     */
    public int getOpposite() {
        return opposite;
    }

    /**
     * Legt den Wert der opposite-Eigenschaft fest.
     * 
     */
    public void setOpposite(int value) {
        this.opposite = value;
    }

    /**
     * Ruft den Wert der index-Eigenschaft ab.
     * 
     */
    public int getIndex() {
        return index;
    }

    /**
     * Legt den Wert der index-Eigenschaft fest.
     * 
     */
    public void setIndex(int value) {
        this.index = value;
    }

}
