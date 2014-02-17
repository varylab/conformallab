//
// Diese Datei wurde mit der JavaTM Architecture for XML Binding(JAXB) Reference Implementation, v2.2.6 generiert 
// Siehe <a href="http://java.sun.com/xml/jaxb">http://java.sun.com/xml/jaxb</a> 
// Änderungen an dieser Datei gehen bei einer Neukompilierung des Quellschemas verloren. 
// Generiert: 2014.02.17 um 03:17:36 PM CET 
//


package de.varylab.conformallab.data.types;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlType;


/**
 * <p>Java-Klasse für Complex complex type.
 * 
 * <p>Das folgende Schemafragment gibt den erwarteten Content an, der in dieser Klasse enthalten ist.
 * 
 * <pre>
 * &lt;complexType name="Complex">
 *   &lt;complexContent>
 *     &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
 *       &lt;attribute name="re" use="required" type="{http://www.w3.org/2001/XMLSchema}double" />
 *       &lt;attribute name="im" use="required" type="{http://www.w3.org/2001/XMLSchema}double" />
 *     &lt;/restriction>
 *   &lt;/complexContent>
 * &lt;/complexType>
 * </pre>
 * 
 * 
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "Complex")
public class Complex {

    @XmlAttribute(name = "re", required = true)
    protected double re;
    @XmlAttribute(name = "im", required = true)
    protected double im;

    /**
     * Ruft den Wert der re-Eigenschaft ab.
     * 
     */
    public double getRe() {
        return re;
    }

    /**
     * Legt den Wert der re-Eigenschaft fest.
     * 
     */
    public void setRe(double value) {
        this.re = value;
    }

    /**
     * Ruft den Wert der im-Eigenschaft ab.
     * 
     */
    public double getIm() {
        return im;
    }

    /**
     * Legt den Wert der im-Eigenschaft fest.
     * 
     */
    public void setIm(double value) {
        this.im = value;
    }

}
