//
// Diese Datei wurde mit der JavaTM Architecture for XML Binding(JAXB) Reference Implementation, v2.2.6 generiert 
// Siehe <a href="http://java.sun.com/xml/jaxb">http://java.sun.com/xml/jaxb</a> 
// Änderungen an dieser Datei gehen bei einer Neukompilierung des Quellschemas verloren. 
// Generiert: 2013.12.10 um 06:00:53 PM CET 
//


package de.varylab.conformallab.data.types;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlType;


/**
 * <p>Java-Klasse für FuchsianData complex type.
 * 
 * <p>Das folgende Schemafragment gibt den erwarteten Content an, der in dieser Klasse enthalten ist.
 * 
 * <pre>
 * &lt;complexType name="FuchsianData">
 *   &lt;complexContent>
 *     &lt;extension base="{http://www.varylab.com/conformallab/types}ConformalData">
 *       &lt;sequence>
 *         &lt;element name="FuchsianGroup" type="{http://www.varylab.com/conformallab/types}FuchsianGroup"/>
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
@XmlType(name = "FuchsianData", propOrder = {
    "fuchsianGroup",
    "fundamentalPolygon"
})
public class FuchsianData
    extends ConformalData
{

    @XmlElement(name = "FuchsianGroup", required = true)
    protected FuchsianGroup fuchsianGroup;
    @XmlElement(name = "FundamentalPolygon", required = true)
    protected FundamentalPolygon fundamentalPolygon;

    /**
     * Ruft den Wert der fuchsianGroup-Eigenschaft ab.
     * 
     * @return
     *     possible object is
     *     {@link FuchsianGroup }
     *     
     */
    public FuchsianGroup getFuchsianGroup() {
        return fuchsianGroup;
    }

    /**
     * Legt den Wert der fuchsianGroup-Eigenschaft fest.
     * 
     * @param value
     *     allowed object is
     *     {@link FuchsianGroup }
     *     
     */
    public void setFuchsianGroup(FuchsianGroup value) {
        this.fuchsianGroup = value;
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
