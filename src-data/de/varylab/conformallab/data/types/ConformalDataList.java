//
// Diese Datei wurde mit der JavaTM Architecture for XML Binding(JAXB) Reference Implementation, v2.2.6 generiert 
// Siehe <a href="http://java.sun.com/xml/jaxb">http://java.sun.com/xml/jaxb</a> 
// Änderungen an dieser Datei gehen bei einer Neukompilierung des Quellschemas verloren. 
// Generiert: 2014.05.30 um 12:40:20 PM CEST 
//


package de.varylab.conformallab.data.types;

import java.util.ArrayList;
import java.util.List;
import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlElements;
import javax.xml.bind.annotation.XmlType;


/**
 * <p>Java-Klasse für ConformalDataList complex type.
 * 
 * <p>Das folgende Schemafragment gibt den erwarteten Content an, der in dieser Klasse enthalten ist.
 * 
 * <pre>
 * &lt;complexType name="ConformalDataList">
 *   &lt;complexContent>
 *     &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
 *       &lt;sequence>
 *         &lt;choice maxOccurs="unbounded" minOccurs="0">
 *           &lt;element name="UniformizationData" type="{http://www.varylab.com/conformallab/types}UniformizationData"/>
 *           &lt;element name="HyperEllipticAlgebraicCurve" type="{http://www.varylab.com/conformallab/types}HyperEllipticAlgebraicCurve"/>
 *           &lt;element name="SchottkyData" type="{http://www.varylab.com/conformallab/types}SchottkyData"/>
 *           &lt;element name="DiscreteMetric" type="{http://www.varylab.com/conformallab/types}DiscreteMetric"/>
 *           &lt;element name="DiscreteEmbedding" type="{http://www.varylab.com/conformallab/types}DiscreteEmbedding"/>
 *           &lt;element name="DiscreteMap" type="{http://www.varylab.com/conformallab/types}DiscreteMap"/>
 *         &lt;/choice>
 *       &lt;/sequence>
 *     &lt;/restriction>
 *   &lt;/complexContent>
 * &lt;/complexType>
 * </pre>
 * 
 * 
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "ConformalDataList", propOrder = {
    "data"
})
public class ConformalDataList {

    @XmlElements({
        @XmlElement(name = "UniformizationData", type = UniformizationData.class),
        @XmlElement(name = "HyperEllipticAlgebraicCurve", type = HyperEllipticAlgebraicCurve.class),
        @XmlElement(name = "SchottkyData", type = SchottkyData.class),
        @XmlElement(name = "DiscreteMetric", type = DiscreteMetric.class),
        @XmlElement(name = "DiscreteEmbedding", type = DiscreteEmbedding.class),
        @XmlElement(name = "DiscreteMap", type = DiscreteMap.class)
    })
    protected List<ConformalData> data;

    /**
     * Gets the value of the data property.
     * 
     * <p>
     * This accessor method returns a reference to the live list,
     * not a snapshot. Therefore any modification you make to the
     * returned list will be present inside the JAXB object.
     * This is why there is not a <CODE>set</CODE> method for the data property.
     * 
     * <p>
     * For example, to add a new item, do as follows:
     * <pre>
     *    getData().add(newItem);
     * </pre>
     * 
     * 
     * <p>
     * Objects of the following type(s) are allowed in the list
     * {@link UniformizationData }
     * {@link HyperEllipticAlgebraicCurve }
     * {@link SchottkyData }
     * {@link DiscreteMetric }
     * {@link DiscreteEmbedding }
     * {@link DiscreteMap }
     * 
     * 
     */
    public List<ConformalData> getData() {
        if (data == null) {
            data = new ArrayList<ConformalData>();
        }
        return this.data;
    }

}
