//
// Diese Datei wurde mit der JavaTM Architecture for XML Binding(JAXB) Reference Implementation, v2.2.6 generiert 
// Siehe <a href="http://java.sun.com/xml/jaxb">http://java.sun.com/xml/jaxb</a> 
// Änderungen an dieser Datei gehen bei einer Neukompilierung des Quellschemas verloren. 
// Generiert: 2014.09.03 um 12:51:40 PM CEST 
//


package de.varylab.conformallab.data.types;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlType;


/**
 * <p>Java-Klasse für UniformizationData complex type.
 * 
 * <p>Das folgende Schemafragment gibt den erwarteten Content an, der in dieser Klasse enthalten ist.
 * 
 * <pre>
 * &lt;complexType name="UniformizationData">
 *   &lt;complexContent>
 *     &lt;extension base="{http://www.varylab.com/conformallab/types}ConformalData">
 *       &lt;sequence>
 *         &lt;element name="UniformizingGroup" type="{http://www.varylab.com/conformallab/types}UniformizingGroup"/>
 *         &lt;element name="FundamentalPolygon" type="{http://www.varylab.com/conformallab/types}FundamentalPolygon"/>
 *       &lt;/sequence>
 *     &lt;/extension>
 *   &lt;/complexContent>
 * &lt;/complexType>
 * </pre>
 * 
 * 
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "UniformizationData", propOrder = {
    "uniformizingGroup",
    "fundamentalPolygon"
})
public class UniformizationData
    extends ConformalData
{

    @XmlElement(name = "UniformizingGroup", required = true)
    protected UniformizingGroup uniformizingGroup;
    @XmlElement(name = "FundamentalPolygon", required = true)
    protected FundamentalPolygon fundamentalPolygon;

    /**
     * Ruft den Wert der uniformizingGroup-Eigenschaft ab.
     * 
     * @return
     *     possible object is
     *     {@link UniformizingGroup }
     *     
     */
    public UniformizingGroup getUniformizingGroup() {
        return uniformizingGroup;
    }

    /**
     * Legt den Wert der uniformizingGroup-Eigenschaft fest.
     * 
     * @param value
     *     allowed object is
     *     {@link UniformizingGroup }
     *     
     */
    public void setUniformizingGroup(UniformizingGroup value) {
        this.uniformizingGroup = value;
    }

    /**
     * Ruft den Wert der fundamentalPolygon-Eigenschaft ab.
     * 
     * @return
     *     possible object is
     *     {@link FundamentalPolygon }
     *     
     */
    public FundamentalPolygon getFundamentalPolygon() {
        return fundamentalPolygon;
    }

    /**
     * Legt den Wert der fundamentalPolygon-Eigenschaft fest.
     * 
     * @param value
     *     allowed object is
     *     {@link FundamentalPolygon }
     *     
     */
    public void setFundamentalPolygon(FundamentalPolygon value) {
        this.fundamentalPolygon = value;
    }

}
