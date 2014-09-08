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
 * <p>Java-Klasse für HalfedgeEmbedding complex type.
 * 
 * <p>Das folgende Schemafragment gibt den erwarteten Content an, der in dieser Klasse enthalten ist.
 * 
 * <pre>
 * &lt;complexType name="HalfedgeEmbedding">
 *   &lt;complexContent>
 *     &lt;extension base="{http://www.varylab.com/conformallab/types}ConformalData">
 *       &lt;sequence>
 *         &lt;element name="Vertex" type="{http://www.varylab.com/conformallab/types}HalfedgeVertex" maxOccurs="unbounded" minOccurs="0"/>
 *         &lt;element name="Identification" type="{http://www.varylab.com/conformallab/types}VertexIdentification" maxOccurs="unbounded" minOccurs="0"/>
 *         &lt;element name="Halfedge" type="{http://www.varylab.com/conformallab/types}HalfedgeEdge" maxOccurs="unbounded" minOccurs="0"/>
 *         &lt;element name="EdgeIdentification" type="{http://www.varylab.com/conformallab/types}EdgeIdentification" maxOccurs="unbounded" minOccurs="0"/>
 *         &lt;element name="Face" type="{http://www.varylab.com/conformallab/types}HalfedgeFace" maxOccurs="unbounded" minOccurs="0"/>
 *         &lt;element name="Selection" type="{http://www.varylab.com/conformallab/types}HalfedgeSelection" minOccurs="0"/>
 *       &lt;/sequence>
 *     &lt;/extension>
 *   &lt;/complexContent>
 * &lt;/complexType>
 * </pre>
 * 
 * 
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "HalfedgeEmbedding", propOrder = {
    "vertices",
    "identifications",
    "edges",
    "edgeIdentifications",
    "faces",
    "selection"
})
public class HalfedgeEmbedding
    extends ConformalData
{

    @XmlElement(name = "Vertex")
    protected List<HalfedgeVertex> vertices;
    @XmlElement(name = "Identification")
    protected List<VertexIdentification> identifications;
    @XmlElement(name = "Halfedge")
    protected List<HalfedgeEdge> edges;
    @XmlElement(name = "EdgeIdentification")
    protected List<EdgeIdentification> edgeIdentifications;
    @XmlElement(name = "Face")
    protected List<HalfedgeFace> faces;
    @XmlElement(name = "Selection")
    protected HalfedgeSelection selection;

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
     * {@link HalfedgeVertex }
     * 
     * 
     */
    public List<HalfedgeVertex> getVertices() {
        if (vertices == null) {
            vertices = new ArrayList<HalfedgeVertex>();
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
     * {@link HalfedgeEdge }
     * 
     * 
     */
    public List<HalfedgeEdge> getEdges() {
        if (edges == null) {
            edges = new ArrayList<HalfedgeEdge>();
        }
        return this.edges;
    }

    /**
     * Gets the value of the edgeIdentifications property.
     * 
     * <p>
     * This accessor method returns a reference to the live list,
     * not a snapshot. Therefore any modification you make to the
     * returned list will be present inside the JAXB object.
     * This is why there is not a <CODE>set</CODE> method for the edgeIdentifications property.
     * 
     * <p>
     * For example, to add a new item, do as follows:
     * <pre>
     *    getEdgeIdentifications().add(newItem);
     * </pre>
     * 
     * 
     * <p>
     * Objects of the following type(s) are allowed in the list
     * {@link EdgeIdentification }
     * 
     * 
     */
    public List<EdgeIdentification> getEdgeIdentifications() {
        if (edgeIdentifications == null) {
            edgeIdentifications = new ArrayList<EdgeIdentification>();
        }
        return this.edgeIdentifications;
    }

    /**
     * Gets the value of the faces property.
     * 
     * <p>
     * This accessor method returns a reference to the live list,
     * not a snapshot. Therefore any modification you make to the
     * returned list will be present inside the JAXB object.
     * This is why there is not a <CODE>set</CODE> method for the faces property.
     * 
     * <p>
     * For example, to add a new item, do as follows:
     * <pre>
     *    getFaces().add(newItem);
     * </pre>
     * 
     * 
     * <p>
     * Objects of the following type(s) are allowed in the list
     * {@link HalfedgeFace }
     * 
     * 
     */
    public List<HalfedgeFace> getFaces() {
        if (faces == null) {
            faces = new ArrayList<HalfedgeFace>();
        }
        return this.faces;
    }

    /**
     * Ruft den Wert der selection-Eigenschaft ab.
     * 
     * @return
     *     possible object is
     *     {@link HalfedgeSelection }
     *     
     */
    public HalfedgeSelection getSelection() {
        return selection;
    }

    /**
     * Legt den Wert der selection-Eigenschaft fest.
     * 
     * @param value
     *     allowed object is
     *     {@link HalfedgeSelection }
     *     
     */
    public void setSelection(HalfedgeSelection value) {
        this.selection = value;
    }

}
