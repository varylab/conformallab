//
// Diese Datei wurde mit der JavaTM Architecture for XML Binding(JAXB) Reference Implementation, v2.2.6 generiert 
// Siehe <a href="http://java.sun.com/xml/jaxb">http://java.sun.com/xml/jaxb</a> 
// Änderungen an dieser Datei gehen bei einer Neukompilierung des Quellschemas verloren. 
// Generiert: 2014.02.26 um 05:27:05 PM CET 
//


package de.varylab.conformallab.data.types;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlType;


/**
 * <p>Java-Klasse für SchottkyGenerator complex type.
 * 
 * <p>Das folgende Schemafragment gibt den erwarteten Content an, der in dieser Klasse enthalten ist.
 * 
 * <pre>
 * &lt;complexType name="SchottkyGenerator">
 *   &lt;complexContent>
 *     &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
 *       &lt;sequence>
 *         &lt;element name="A" type="{http://www.varylab.com/conformallab/types}Complex"/>
 *         &lt;element name="B" type="{http://www.varylab.com/conformallab/types}Complex"/>
 *         &lt;element name="Mu" type="{http://www.varylab.com/conformallab/types}Complex"/>
 *         &lt;element name="Circle" type="{http://www.varylab.com/conformallab/types}Circle"/>
 *       &lt;/sequence>
 *     &lt;/restriction>
 *   &lt;/complexContent>
 * &lt;/complexType>
 * </pre>
 * 
 * 
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "SchottkyGenerator", propOrder = {
    "a",
    "b",
    "mu",
    "circle"
})
public class SchottkyGenerator {

    @XmlElement(name = "A", required = true)
    protected Complex a;
    @XmlElement(name = "B", required = true)
    protected Complex b;
    @XmlElement(name = "Mu", required = true)
    protected Complex mu;
    @XmlElement(name = "Circle", required = true)
    protected Circle circle;

    /**
     * Ruft den Wert der a-Eigenschaft ab.
     * 
     * @return
     *     possible object is
     *     {@link Complex }
     *     
     */
    public Complex getA() {
        return a;
    }

    /**
     * Legt den Wert der a-Eigenschaft fest.
     * 
     * @param value
     *     allowed object is
     *     {@link Complex }
     *     
     */
    public void setA(Complex value) {
        this.a = value;
    }

    /**
     * Ruft den Wert der b-Eigenschaft ab.
     * 
     * @return
     *     possible object is
     *     {@link Complex }
     *     
     */
    public Complex getB() {
        return b;
    }

    /**
     * Legt den Wert der b-Eigenschaft fest.
     * 
     * @param value
     *     allowed object is
     *     {@link Complex }
     *     
     */
    public void setB(Complex value) {
        this.b = value;
    }

    /**
     * Ruft den Wert der mu-Eigenschaft ab.
     * 
     * @return
     *     possible object is
     *     {@link Complex }
     *     
     */
    public Complex getMu() {
        return mu;
    }

    /**
     * Legt den Wert der mu-Eigenschaft fest.
     * 
     * @param value
     *     allowed object is
     *     {@link Complex }
     *     
     */
    public void setMu(Complex value) {
        this.mu = value;
    }

    /**
     * Ruft den Wert der circle-Eigenschaft ab.
     * 
     * @return
     *     possible object is
     *     {@link Circle }
     *     
     */
    public Circle getCircle() {
        return circle;
    }

    /**
     * Legt den Wert der circle-Eigenschaft fest.
     * 
     * @param value
     *     allowed object is
     *     {@link Circle }
     *     
     */
    public void setCircle(Circle value) {
        this.circle = value;
    }

}
