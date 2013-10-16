//
// Diese Datei wurde mit der JavaTM Architecture for XML Binding(JAXB) Reference Implementation, v2.2.6 generiert 
// Siehe <a href="http://java.sun.com/xml/jaxb">http://java.sun.com/xml/jaxb</a> 
// \u00c4nderungen an dieser Datei gehen bei einer Neukompilierung des Quellschemas verloren. 
// Generiert: 2013.10.10 um 12:59:01 PM CEST 
//


package de.varylab.conformallab.data.types;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlType;


/**
 * <p>Java-Klasse f\u00fcr SchottkyGenerator complex type.
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
 *         &lt;element name="mu" type="{http://www.varylab.com/conformallab/types}Complex"/>
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
    "mu"
})
public class SchottkyGenerator {

    @XmlElement(name = "A", required = true)
    protected Complex a;
    @XmlElement(name = "B", required = true)
    protected Complex b;
    @XmlElement(required = true)
    protected Complex mu;

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

}
