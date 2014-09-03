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
import javax.xml.bind.annotation.XmlType;


/**
 * <p>Java-Klasse für IsometryPSL2R complex type.
 * 
 * <p>Das folgende Schemafragment gibt den erwarteten Content an, der in dieser Klasse enthalten ist.
 * 
 * <pre>
 * &lt;complexType name="IsometryPSL2R">
 *   &lt;complexContent>
 *     &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
 *       &lt;attribute name="m11" type="{http://www.w3.org/2001/XMLSchema}double" default="0.0" />
 *       &lt;attribute name="m12" type="{http://www.w3.org/2001/XMLSchema}double" default="0.0" />
 *       &lt;attribute name="m13" type="{http://www.w3.org/2001/XMLSchema}double" default="0.0" />
 *       &lt;attribute name="m21" type="{http://www.w3.org/2001/XMLSchema}double" default="0.0" />
 *       &lt;attribute name="m22" type="{http://www.w3.org/2001/XMLSchema}double" default="0.0" />
 *       &lt;attribute name="m23" type="{http://www.w3.org/2001/XMLSchema}double" default="0.0" />
 *       &lt;attribute name="m31" type="{http://www.w3.org/2001/XMLSchema}double" default="0.0" />
 *       &lt;attribute name="m32" type="{http://www.w3.org/2001/XMLSchema}double" default="0.0" />
 *       &lt;attribute name="m33" type="{http://www.w3.org/2001/XMLSchema}double" default="0.0" />
 *     &lt;/restriction>
 *   &lt;/complexContent>
 * &lt;/complexType>
 * </pre>
 * 
 * 
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "IsometryPSL2R")
public class IsometryPSL2R {

    @XmlAttribute(name = "m11")
    protected Double m11;
    @XmlAttribute(name = "m12")
    protected Double m12;
    @XmlAttribute(name = "m13")
    protected Double m13;
    @XmlAttribute(name = "m21")
    protected Double m21;
    @XmlAttribute(name = "m22")
    protected Double m22;
    @XmlAttribute(name = "m23")
    protected Double m23;
    @XmlAttribute(name = "m31")
    protected Double m31;
    @XmlAttribute(name = "m32")
    protected Double m32;
    @XmlAttribute(name = "m33")
    protected Double m33;

    /**
     * Ruft den Wert der m11-Eigenschaft ab.
     * 
     * @return
     *     possible object is
     *     {@link Double }
     *     
     */
    public double getM11() {
        if (m11 == null) {
            return  0.0D;
        } else {
            return m11;
        }
    }

    /**
     * Legt den Wert der m11-Eigenschaft fest.
     * 
     * @param value
     *     allowed object is
     *     {@link Double }
     *     
     */
    public void setM11(Double value) {
        this.m11 = value;
    }

    /**
     * Ruft den Wert der m12-Eigenschaft ab.
     * 
     * @return
     *     possible object is
     *     {@link Double }
     *     
     */
    public double getM12() {
        if (m12 == null) {
            return  0.0D;
        } else {
            return m12;
        }
    }

    /**
     * Legt den Wert der m12-Eigenschaft fest.
     * 
     * @param value
     *     allowed object is
     *     {@link Double }
     *     
     */
    public void setM12(Double value) {
        this.m12 = value;
    }

    /**
     * Ruft den Wert der m13-Eigenschaft ab.
     * 
     * @return
     *     possible object is
     *     {@link Double }
     *     
     */
    public double getM13() {
        if (m13 == null) {
            return  0.0D;
        } else {
            return m13;
        }
    }

    /**
     * Legt den Wert der m13-Eigenschaft fest.
     * 
     * @param value
     *     allowed object is
     *     {@link Double }
     *     
     */
    public void setM13(Double value) {
        this.m13 = value;
    }

    /**
     * Ruft den Wert der m21-Eigenschaft ab.
     * 
     * @return
     *     possible object is
     *     {@link Double }
     *     
     */
    public double getM21() {
        if (m21 == null) {
            return  0.0D;
        } else {
            return m21;
        }
    }

    /**
     * Legt den Wert der m21-Eigenschaft fest.
     * 
     * @param value
     *     allowed object is
     *     {@link Double }
     *     
     */
    public void setM21(Double value) {
        this.m21 = value;
    }

    /**
     * Ruft den Wert der m22-Eigenschaft ab.
     * 
     * @return
     *     possible object is
     *     {@link Double }
     *     
     */
    public double getM22() {
        if (m22 == null) {
            return  0.0D;
        } else {
            return m22;
        }
    }

    /**
     * Legt den Wert der m22-Eigenschaft fest.
     * 
     * @param value
     *     allowed object is
     *     {@link Double }
     *     
     */
    public void setM22(Double value) {
        this.m22 = value;
    }

    /**
     * Ruft den Wert der m23-Eigenschaft ab.
     * 
     * @return
     *     possible object is
     *     {@link Double }
     *     
     */
    public double getM23() {
        if (m23 == null) {
            return  0.0D;
        } else {
            return m23;
        }
    }

    /**
     * Legt den Wert der m23-Eigenschaft fest.
     * 
     * @param value
     *     allowed object is
     *     {@link Double }
     *     
     */
    public void setM23(Double value) {
        this.m23 = value;
    }

    /**
     * Ruft den Wert der m31-Eigenschaft ab.
     * 
     * @return
     *     possible object is
     *     {@link Double }
     *     
     */
    public double getM31() {
        if (m31 == null) {
            return  0.0D;
        } else {
            return m31;
        }
    }

    /**
     * Legt den Wert der m31-Eigenschaft fest.
     * 
     * @param value
     *     allowed object is
     *     {@link Double }
     *     
     */
    public void setM31(Double value) {
        this.m31 = value;
    }

    /**
     * Ruft den Wert der m32-Eigenschaft ab.
     * 
     * @return
     *     possible object is
     *     {@link Double }
     *     
     */
    public double getM32() {
        if (m32 == null) {
            return  0.0D;
        } else {
            return m32;
        }
    }

    /**
     * Legt den Wert der m32-Eigenschaft fest.
     * 
     * @param value
     *     allowed object is
     *     {@link Double }
     *     
     */
    public void setM32(Double value) {
        this.m32 = value;
    }

    /**
     * Ruft den Wert der m33-Eigenschaft ab.
     * 
     * @return
     *     possible object is
     *     {@link Double }
     *     
     */
    public double getM33() {
        if (m33 == null) {
            return  0.0D;
        } else {
            return m33;
        }
    }

    /**
     * Legt den Wert der m33-Eigenschaft fest.
     * 
     * @param value
     *     allowed object is
     *     {@link Double }
     *     
     */
    public void setM33(Double value) {
        this.m33 = value;
    }

}
