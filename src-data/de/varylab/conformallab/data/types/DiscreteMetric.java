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
import javax.xml.bind.annotation.XmlType;


/**
 * <p>Java-Klasse für DiscreteMetric complex type.
 * 
 * <p>Das folgende Schemafragment gibt den erwarteten Content an, der in dieser Klasse enthalten ist.
 * 
 * <pre>
 * &lt;complexType name="DiscreteMetric">
 *   &lt;complexContent>
 *     &lt;extension base="{http://www.varylab.com/conformallab/types}ConformalData">
 *       &lt;sequence>
 *         &lt;element name="MetricEdge" type="{http://www.varylab.com/conformallab/types}MetricEdge" maxOccurs="unbounded" minOccurs="0"/>
 *         &lt;element name="MetricTriangle" type="{http://www.varylab.com/conformallab/types}MetricTriangle" maxOccurs="unbounded" minOccurs="0"/>
 *       &lt;/sequence>
 *     &lt;/extension>
 *   &lt;/complexContent>
 * &lt;/complexType>
 * </pre>
 * 
 * 
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "DiscreteMetric", propOrder = {
    "edges",
    "triangles"
})
public class DiscreteMetric
    extends ConformalData
{

    @XmlElement(name = "MetricEdge")
    protected List<MetricEdge> edges;
    @XmlElement(name = "MetricTriangle")
    protected List<MetricTriangle> triangles;

    /**
     * Gets the value of the edges property.
     * 
     * <p>
     * This accessor method returns a reference to the live list,
     * not a snapshot. Therefore any modification you make to the
     * returned list will be present inside the JAXB object.
     * This is why there is not a <CODE>set</CODE> method for the edges property.
     * 
     * <p>
     * For example, to add a new item, do as follows:
     * <pre>
     *    getEdges().add(newItem);
     * </pre>
     * 
     * 
     * <p>
     * Objects of the following type(s) are allowed in the list
     * {@link MetricEdge }
     * 
     * 
     */
    public List<MetricEdge> getEdges() {
        if (edges == null) {
            edges = new ArrayList<MetricEdge>();
        }
        return this.edges;
    }

    /**
     * Gets the value of the triangles property.
     * 
     * <p>
     * This accessor method returns a reference to the live list,
     * not a snapshot. Therefore any modification you make to the
     * returned list will be present inside the JAXB object.
     * This is why there is not a <CODE>set</CODE> method for the triangles property.
     * 
     * <p>
     * For example, to add a new item, do as follows:
     * <pre>
     *    getTriangles().add(newItem);
     * </pre>
     * 
     * 
     * <p>
     * Objects of the following type(s) are allowed in the list
     * {@link MetricTriangle }
     * 
     * 
     */
    public List<MetricTriangle> getTriangles() {
        if (triangles == null) {
            triangles = new ArrayList<MetricTriangle>();
        }
        return this.triangles;
    }

}
