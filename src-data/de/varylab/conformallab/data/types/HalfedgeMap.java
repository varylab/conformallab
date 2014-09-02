//
// Diese Datei wurde mit der JavaTM Architecture for XML Binding(JAXB) Reference Implementation, v2.2.6 generiert 
// Siehe <a href="http://java.sun.com/xml/jaxb">http://java.sun.com/xml/jaxb</a> 
// Änderungen an dieser Datei gehen bei einer Neukompilierung des Quellschemas verloren. 
// Generiert: 2014.09.02 um 04:53:14 PM CEST 
//


package de.varylab.conformallab.data.types;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlType;


/**
 * <p>Java-Klasse für HalfedgeMap complex type.
 * 
 * <p>Das folgende Schemafragment gibt den erwarteten Content an, der in dieser Klasse enthalten ist.
 * 
 * <pre>
 * &lt;complexType name="HalfedgeMap">
 *   &lt;complexContent>
 *     &lt;extension base="{http://www.varylab.com/conformallab/types}ConformalData">
 *       &lt;sequence>
 *         &lt;element name="Domain" type="{http://www.varylab.com/conformallab/types}HalfedgeEmbedding"/>
 *         &lt;element name="Image" type="{http://www.varylab.com/conformallab/types}HalfedgeEmbedding"/>
 *       &lt;/sequence>
 *     &lt;/extension>
 *   &lt;/complexContent>
 * &lt;/complexType>
 * </pre>
 * 
 * 
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "HalfedgeMap", propOrder = {
    "domain",
    "image"
})
public class HalfedgeMap
    extends ConformalData
{

    @XmlElement(name = "Domain", required = true)
    protected HalfedgeEmbedding domain;
    @XmlElement(name = "Image", required = true)
    protected HalfedgeEmbedding image;

    /**
     * Ruft den Wert der domain-Eigenschaft ab.
     * 
     * @return
     *     possible object is
     *     {@link HalfedgeEmbedding }
     *     
     */
    public HalfedgeEmbedding getDomain() {
        return domain;
    }

    /**
     * Legt den Wert der domain-Eigenschaft fest.
     * 
     * @param value
     *     allowed object is
     *     {@link HalfedgeEmbedding }
     *     
     */
    public void setDomain(HalfedgeEmbedding value) {
        this.domain = value;
    }

    /**
     * Ruft den Wert der image-Eigenschaft ab.
     * 
     * @return
     *     possible object is
     *     {@link HalfedgeEmbedding }
     *     
     */
    public HalfedgeEmbedding getImage() {
        return image;
    }

    /**
     * Legt den Wert der image-Eigenschaft fest.
     * 
     * @param value
     *     allowed object is
     *     {@link HalfedgeEmbedding }
     *     
     */
    public void setImage(HalfedgeEmbedding value) {
        this.image = value;
    }

}
