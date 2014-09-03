//
// Diese Datei wurde mit der JavaTM Architecture for XML Binding(JAXB) Reference Implementation, v2.2.6 generiert 
// Siehe <a href="http://java.sun.com/xml/jaxb">http://java.sun.com/xml/jaxb</a> 
// Änderungen an dieser Datei gehen bei einer Neukompilierung des Quellschemas verloren. 
// Generiert: 2014.09.03 um 04:35:53 PM CEST 
//


package de.varylab.conformallab.data.types;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
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
 *       &lt;sequence>
 *         &lt;element name="StartPosition" type="{http://www.varylab.com/conformallab/types}Complex"/>
 *       &lt;/sequence>
 *       &lt;attribute name="index" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
 *       &lt;attribute name="nextEdge" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
 *       &lt;attribute name="previousEdge" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
 *       &lt;attribute name="identifiedEdge" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
 *       &lt;attribute name="startVertex" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
 *     &lt;/restriction>
 *   &lt;/complexContent>
 * &lt;/complexType>
 * </pre>
 * 
 * 
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "FundamentalEdge", propOrder = {
    "startPosition"
})
public class FundamentalEdge {

    @XmlElement(name = "StartPosition", required = true)
    protected Complex startPosition;
    @XmlAttribute(name = "index", required = true)
    protected int index;
    @XmlAttribute(name = "nextEdge", required = true)
    protected int nextEdge;
    @XmlAttribute(name = "previousEdge", required = true)
    protected int previousEdge;
    @XmlAttribute(name = "identifiedEdge", required = true)
    protected int identifiedEdge;
    @XmlAttribute(name = "startVertex", required = true)
    protected int startVertex;

    /**
     * Ruft den Wert der startPosition-Eigenschaft ab.
     * 
     * @return
     *     possible object is
     *     {@link Complex }
     *     
     */
    public Complex getStartPosition() {
        return startPosition;
    }

    /**
     * Legt den Wert der startPosition-Eigenschaft fest.
     * 
     * @param value
     *     allowed object is
     *     {@link Complex }
     *     
     */
    public void setStartPosition(Complex value) {
        this.startPosition = value;
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

    /**
     * Ruft den Wert der nextEdge-Eigenschaft ab.
     * 
     */
    public int getNextEdge() {
        return nextEdge;
    }

    /**
     * Legt den Wert der nextEdge-Eigenschaft fest.
     * 
     */
    public void setNextEdge(int value) {
        this.nextEdge = value;
    }

    /**
     * Ruft den Wert der previousEdge-Eigenschaft ab.
     * 
     */
    public int getPreviousEdge() {
        return previousEdge;
    }

    /**
     * Legt den Wert der previousEdge-Eigenschaft fest.
     * 
     */
    public void setPreviousEdge(int value) {
        this.previousEdge = value;
    }

    /**
     * Ruft den Wert der identifiedEdge-Eigenschaft ab.
     * 
     */
    public int getIdentifiedEdge() {
        return identifiedEdge;
    }

    /**
     * Legt den Wert der identifiedEdge-Eigenschaft fest.
     * 
     */
    public void setIdentifiedEdge(int value) {
        this.identifiedEdge = value;
    }

    /**
     * Ruft den Wert der startVertex-Eigenschaft ab.
     * 
     */
    public int getStartVertex() {
        return startVertex;
    }

    /**
     * Legt den Wert der startVertex-Eigenschaft fest.
     * 
     */
    public void setStartVertex(int value) {
        this.startVertex = value;
    }

}
