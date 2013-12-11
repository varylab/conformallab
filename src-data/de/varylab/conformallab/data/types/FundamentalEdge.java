//
// Diese Datei wurde mit der JavaTM Architecture for XML Binding(JAXB) Reference Implementation, v2.2.6 generiert 
// Siehe <a href="http://java.sun.com/xml/jaxb">http://java.sun.com/xml/jaxb</a> 
// Änderungen an dieser Datei gehen bei einer Neukompilierung des Quellschemas verloren. 
// Generiert: 2013.12.10 um 06:00:53 PM CET 
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
 *       &lt;attribute name="identifiedEdge" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
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
    @XmlAttribute(name = "identifiedEdge", required = true)
    protected int identifiedEdge;

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

}
