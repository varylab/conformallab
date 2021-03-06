//
// Diese Datei wurde mit der JavaTM Architecture for XML Binding(JAXB) Reference Implementation, v2.2.6 generiert 
// Siehe <a href="http://java.sun.com/xml/jaxb">http://java.sun.com/xml/jaxb</a> 
// Änderungen an dieser Datei gehen bei einer Neukompilierung des Quellschemas verloren. 
// Generiert: 2014.09.03 um 04:35:53 PM CEST 
//


package de.varylab.conformallab.data.types;

import java.util.ArrayList;
import java.util.List;
import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlType;


/**
 * <p>Java-Klasse für DiscreteEmbedding complex type.
 * 
 * <p>Das folgende Schemafragment gibt den erwarteten Content an, der in dieser Klasse enthalten ist.
 * 
 * <pre>
 * &lt;complexType name="DiscreteEmbedding">
 *   &lt;complexContent>
 *     &lt;extension base="{http://www.varylab.com/conformallab/types}ConformalData">
 *       &lt;sequence>
 *         &lt;element name="EmbeddedVertex" type="{http://www.varylab.com/conformallab/types}EmbeddedVertex" maxOccurs="unbounded" minOccurs="0"/>
 *         &lt;element name="Identification" type="{http://www.varylab.com/conformallab/types}VertexIdentification" maxOccurs="unbounded" minOccurs="0"/>
 *         &lt;element name="EmbeddedTriangle" type="{http://www.varylab.com/conformallab/types}EmbeddedTriangle" maxOccurs="unbounded" minOccurs="0"/>
 *         &lt;element name="Selection" type="{http://www.varylab.com/conformallab/types}EmbeddingSelection" minOccurs="0"/>
 *       &lt;/sequence>
 *     &lt;/extension>
 *   &lt;/complexContent>
 * &lt;/complexType>
 * </pre>
 * 
 * 
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "DiscreteEmbedding", propOrder = {
    "vertices",
    "identifications",
    "triangles",
    "selection"
})
public class DiscreteEmbedding
    extends ConformalData
{

    @XmlElement(name = "EmbeddedVertex")
    protected List<EmbeddedVertex> vertices;
    @XmlElement(name = "Identification")
    protected List<VertexIdentification> identifications;
    @XmlElement(name = "EmbeddedTriangle")
    protected List<EmbeddedTriangle> triangles;
    @XmlElement(name = "Selection")
    protected EmbeddingSelection selection;

    /**
     * Gets the value of the vertices property.
     * 
     * <p>
     * This accessor method returns a reference to the live list,
     * not a snapshot. Therefore any modification you make to the
     * returned list will be present inside the JAXB object.
     * This is why there is not a <CODE>set</CODE> method for the vertices property.
     * 
     * <p>
     * For example, to add a new item, do as follows:
     * <pre>
     *    getVertices().add(newItem);
     * </pre>
     * 
     * 
     * <p>
     * Objects of the following type(s) are allowed in the list
     * {@link EmbeddedVertex }
     * 
     * 
     */
    public List<EmbeddedVertex> getVertices() {
        if (vertices == null) {
            vertices = new ArrayList<EmbeddedVertex>();
        }
        return this.vertices;
    }

    /**
     * Gets the value of the identifications property.
     * 
     * <p>
     * This accessor method returns a reference to the live list,
     * not a snapshot. Therefore any modification you make to the
     * returned list will be present inside the JAXB object.
     * This is why there is not a <CODE>set</CODE> method for the identifications property.
     * 
     * <p>
     * For example, to add a new item, do as follows:
     * <pre>
     *    getIdentifications().add(newItem);
     * </pre>
     * 
     * 
     * <p>
     * Objects of the following type(s) are allowed in the list
     * {@link VertexIdentification }
     * 
     * 
     */
    public List<VertexIdentification> getIdentifications() {
        if (identifications == null) {
            identifications = new ArrayList<VertexIdentification>();
        }
        return this.identifications;
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
     * {@link EmbeddedTriangle }
     * 
     * 
     */
    public List<EmbeddedTriangle> getTriangles() {
        if (triangles == null) {
            triangles = new ArrayList<EmbeddedTriangle>();
        }
        return this.triangles;
    }

    /**
     * Ruft den Wert der selection-Eigenschaft ab.
     * 
     * @return
     *     possible object is
     *     {@link EmbeddingSelection }
     *     
     */
    public EmbeddingSelection getSelection() {
        return selection;
    }

    /**
     * Legt den Wert der selection-Eigenschaft fest.
     * 
     * @param value
     *     allowed object is
     *     {@link EmbeddingSelection }
     *     
     */
    public void setSelection(EmbeddingSelection value) {
        this.selection = value;
    }

}
