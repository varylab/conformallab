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
 * <p>Java-Klasse für HyperbolicMotion complex type.
 * 
 * <p>Das folgende Schemafragment gibt den erwarteten Content an, der in dieser Klasse enthalten ist.
 * 
 * <pre>
 * &lt;complexType name="HyperbolicMotion">
 *   &lt;complexContent>
 *     &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
 *       &lt;attribute name="m11" use="required" type="{http://www.w3.org/2001/XMLSchema}double" />
 *       &lt;attribute name="m12" use="required" type="{http://www.w3.org/2001/XMLSchema}double" />
 *       &lt;attribute name="m13" use="required" type="{http://www.w3.org/2001/XMLSchema}double" />
 *       &lt;attribute name="m14" use="required" type="{http://www.w3.org/2001/XMLSchema}double" />
 *       &lt;attribute name="m21" use="required" type="{http://www.w3.org/2001/XMLSchema}double" />
 *       &lt;attribute name="m22" use="required" type="{http://www.w3.org/2001/XMLSchema}double" />
 *       &lt;attribute name="m23" use="required" type="{http://www.w3.org/2001/XMLSchema}double" />
 *       &lt;attribute name="m24" use="required" type="{http://www.w3.org/2001/XMLSchema}double" />
 *       &lt;attribute name="m31" use="required" type="{http://www.w3.org/2001/XMLSchema}double" />
 *       &lt;attribute name="m32" use="required" type="{http://www.w3.org/2001/XMLSchema}double" />
 *       &lt;attribute name="m33" use="required" type="{http://www.w3.org/2001/XMLSchema}double" />
 *       &lt;attribute name="m34" use="required" type="{http://www.w3.org/2001/XMLSchema}double" />
 *       &lt;attribute name="m41" use="required" type="{http://www.w3.org/2001/XMLSchema}double" />
 *       &lt;attribute name="m42" use="required" type="{http://www.w3.org/2001/XMLSchema}double" />
 *       &lt;attribute name="m43" use="required" type="{http://www.w3.org/2001/XMLSchema}double" />
 *       &lt;attribute name="m44" use="required" type="{http://www.w3.org/2001/XMLSchema}double" />
 *     &lt;/restriction>
 *   &lt;/complexContent>
 * &lt;/complexType>
 * </pre>
 * 
 * 
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "HyperbolicMotion")
public class HyperbolicMotion {

    @XmlAttribute(name = "m11", required = true)
    protected double m11;
    @XmlAttribute(name = "m12", required = true)
    protected double m12;
    @XmlAttribute(name = "m13", required = true)
    protected double m13;
    @XmlAttribute(name = "m14", required = true)
    protected double m14;
    @XmlAttribute(name = "m21", required = true)
    protected double m21;
    @XmlAttribute(name = "m22", required = true)
    protected double m22;
    @XmlAttribute(name = "m23", required = true)
    protected double m23;
    @XmlAttribute(name = "m24", required = true)
    protected double m24;
    @XmlAttribute(name = "m31", required = true)
    protected double m31;
    @XmlAttribute(name = "m32", required = true)
    protected double m32;
    @XmlAttribute(name = "m33", required = true)
    protected double m33;
    @XmlAttribute(name = "m34", required = true)
    protected double m34;
    @XmlAttribute(name = "m41", required = true)
    protected double m41;
    @XmlAttribute(name = "m42", required = true)
    protected double m42;
    @XmlAttribute(name = "m43", required = true)
    protected double m43;
    @XmlAttribute(name = "m44", required = true)
    protected double m44;

    /**
     * Ruft den Wert der m11-Eigenschaft ab.
     * 
     */
    public double getM11() {
        return m11;
    }

    /**
     * Legt den Wert der m11-Eigenschaft fest.
     * 
     */
    public void setM11(double value) {
        this.m11 = value;
    }

    /**
     * Ruft den Wert der m12-Eigenschaft ab.
     * 
     */
    public double getM12() {
        return m12;
    }

    /**
     * Legt den Wert der m12-Eigenschaft fest.
     * 
     */
    public void setM12(double value) {
        this.m12 = value;
    }

    /**
     * Ruft den Wert der m13-Eigenschaft ab.
     * 
     */
    public double getM13() {
        return m13;
    }

    /**
     * Legt den Wert der m13-Eigenschaft fest.
     * 
     */
    public void setM13(double value) {
        this.m13 = value;
    }

    /**
     * Ruft den Wert der m14-Eigenschaft ab.
     * 
     */
    public double getM14() {
        return m14;
    }

    /**
     * Legt den Wert der m14-Eigenschaft fest.
     * 
     */
    public void setM14(double value) {
        this.m14 = value;
    }

    /**
     * Ruft den Wert der m21-Eigenschaft ab.
     * 
     */
    public double getM21() {
        return m21;
    }

    /**
     * Legt den Wert der m21-Eigenschaft fest.
     * 
     */
    public void setM21(double value) {
        this.m21 = value;
    }

    /**
     * Ruft den Wert der m22-Eigenschaft ab.
     * 
     */
    public double getM22() {
        return m22;
    }

    /**
     * Legt den Wert der m22-Eigenschaft fest.
     * 
     */
    public void setM22(double value) {
        this.m22 = value;
    }

    /**
     * Ruft den Wert der m23-Eigenschaft ab.
     * 
     */
    public double getM23() {
        return m23;
    }

    /**
     * Legt den Wert der m23-Eigenschaft fest.
     * 
     */
    public void setM23(double value) {
        this.m23 = value;
    }

    /**
     * Ruft den Wert der m24-Eigenschaft ab.
     * 
     */
    public double getM24() {
        return m24;
    }

    /**
     * Legt den Wert der m24-Eigenschaft fest.
     * 
     */
    public void setM24(double value) {
        this.m24 = value;
    }

    /**
     * Ruft den Wert der m31-Eigenschaft ab.
     * 
     */
    public double getM31() {
        return m31;
    }

    /**
     * Legt den Wert der m31-Eigenschaft fest.
     * 
     */
    public void setM31(double value) {
        this.m31 = value;
    }

    /**
     * Ruft den Wert der m32-Eigenschaft ab.
     * 
     */
    public double getM32() {
        return m32;
    }

    /**
     * Legt den Wert der m32-Eigenschaft fest.
     * 
     */
    public void setM32(double value) {
        this.m32 = value;
    }

    /**
     * Ruft den Wert der m33-Eigenschaft ab.
     * 
     */
    public double getM33() {
        return m33;
    }

    /**
     * Legt den Wert der m33-Eigenschaft fest.
     * 
     */
    public void setM33(double value) {
        this.m33 = value;
    }

    /**
     * Ruft den Wert der m34-Eigenschaft ab.
     * 
     */
    public double getM34() {
        return m34;
    }

    /**
     * Legt den Wert der m34-Eigenschaft fest.
     * 
     */
    public void setM34(double value) {
        this.m34 = value;
    }

    /**
     * Ruft den Wert der m41-Eigenschaft ab.
     * 
     */
    public double getM41() {
        return m41;
    }

    /**
     * Legt den Wert der m41-Eigenschaft fest.
     * 
     */
    public void setM41(double value) {
        this.m41 = value;
    }

    /**
     * Ruft den Wert der m42-Eigenschaft ab.
     * 
     */
    public double getM42() {
        return m42;
    }

    /**
     * Legt den Wert der m42-Eigenschaft fest.
     * 
     */
    public void setM42(double value) {
        this.m42 = value;
    }

    /**
     * Ruft den Wert der m43-Eigenschaft ab.
     * 
     */
    public double getM43() {
        return m43;
    }

    /**
     * Legt den Wert der m43-Eigenschaft fest.
     * 
     */
    public void setM43(double value) {
        this.m43 = value;
    }

    /**
     * Ruft den Wert der m44-Eigenschaft ab.
     * 
     */
    public double getM44() {
        return m44;
    }

    /**
     * Legt den Wert der m44-Eigenschaft fest.
     * 
     */
    public void setM44(double value) {
        this.m44 = value;
    }

}
