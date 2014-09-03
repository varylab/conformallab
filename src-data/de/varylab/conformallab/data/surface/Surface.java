//
// Diese Datei wurde mit der JavaTM Architecture for XML Binding(JAXB) Reference Implementation, v2.2.6 generiert 
// Siehe <a href="http://java.sun.com/xml/jaxb">http://java.sun.com/xml/jaxb</a> 
// Änderungen an dieser Datei gehen bei einer Neukompilierung des Quellschemas verloren. 
// Generiert: 2014.09.03 um 04:35:53 PM CEST 
//


package de.varylab.conformallab.data.surface;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlType;
import de.varylab.conformallab.data.types.SchottkyData;


/**
 * <p>Java-Klasse für Surface complex type.
 * 
 * <p>Das folgende Schemafragment gibt den erwarteten Content an, der in dieser Klasse enthalten ist.
 * 
 * <pre>
 * &lt;complexType name="Surface">
 *   &lt;complexContent>
 *     &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
 *       &lt;choice>
 *         &lt;element name="SchottkyData" type="{http://www.varylab.com/conformallab/types}SchottkyData"/>
 *       &lt;/choice>
 *     &lt;/restriction>
 *   &lt;/complexContent>
 * &lt;/complexType>
 * </pre>
 * 
 * 
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "Surface", propOrder = {
    "schottkyData"
})
public class Surface {

    @XmlElement(name = "SchottkyData")
    protected SchottkyData schottkyData;

    /**
     * Ruft den Wert der schottkyData-Eigenschaft ab.
     * 
     * @return
     *     possible object is
     *     {@link SchottkyData }
     *     
     */
    public SchottkyData getSchottkyData() {
        return schottkyData;
    }

    /**
     * Legt den Wert der schottkyData-Eigenschaft fest.
     * 
     * @param value
     *     allowed object is
     *     {@link SchottkyData }
     *     
     */
    public void setSchottkyData(SchottkyData value) {
        this.schottkyData = value;
    }

}
