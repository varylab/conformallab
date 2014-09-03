//
// Diese Datei wurde mit der JavaTM Architecture for XML Binding(JAXB) Reference Implementation, v2.2.6 generiert 
// Siehe <a href="http://java.sun.com/xml/jaxb">http://java.sun.com/xml/jaxb</a> 
// Änderungen an dieser Datei gehen bei einer Neukompilierung des Quellschemas verloren. 
// Generiert: 2014.09.03 um 12:51:40 PM CEST 
//


package de.varylab.conformallab.data.types;

import java.util.ArrayList;
import java.util.List;
import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlType;


/**
 * <p>Java-Klasse für HalfedgeSelection complex type.
 * 
 * <p>Das folgende Schemafragment gibt den erwarteten Content an, der in dieser Klasse enthalten ist.
 * 
 * <pre>
 * &lt;complexType name="HalfedgeSelection">
 *   &lt;complexContent>
 *     &lt;extension base="{http://www.varylab.com/conformallab/types}ConformalData">
 *       &lt;sequence>
 *         &lt;element name="VertexSelection">
 *           &lt;complexType>
 *             &lt;complexContent>
 *               &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
 *                 &lt;sequence>
 *                   &lt;element name="Vertex" maxOccurs="unbounded" minOccurs="0">
 *                     &lt;complexType>
 *                       &lt;complexContent>
 *                         &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
 *                           &lt;attribute name="index" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
 *                           &lt;attribute name="channel" type="{http://www.w3.org/2001/XMLSchema}int" default="0" />
 *                         &lt;/restriction>
 *                       &lt;/complexContent>
 *                     &lt;/complexType>
 *                   &lt;/element>
 *                 &lt;/sequence>
 *               &lt;/restriction>
 *             &lt;/complexContent>
 *           &lt;/complexType>
 *         &lt;/element>
 *         &lt;element name="EdgeSelection">
 *           &lt;complexType>
 *             &lt;complexContent>
 *               &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
 *                 &lt;sequence>
 *                   &lt;element name="Edge" maxOccurs="unbounded" minOccurs="0">
 *                     &lt;complexType>
 *                       &lt;complexContent>
 *                         &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
 *                           &lt;attribute name="index" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
 *                           &lt;attribute name="channel" type="{http://www.w3.org/2001/XMLSchema}int" default="0" />
 *                         &lt;/restriction>
 *                       &lt;/complexContent>
 *                     &lt;/complexType>
 *                   &lt;/element>
 *                 &lt;/sequence>
 *               &lt;/restriction>
 *             &lt;/complexContent>
 *           &lt;/complexType>
 *         &lt;/element>
 *         &lt;element name="FaceSelection">
 *           &lt;complexType>
 *             &lt;complexContent>
 *               &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
 *                 &lt;sequence>
 *                   &lt;element name="Face" maxOccurs="unbounded" minOccurs="0">
 *                     &lt;complexType>
 *                       &lt;complexContent>
 *                         &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
 *                           &lt;attribute name="index" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
 *                           &lt;attribute name="channel" type="{http://www.w3.org/2001/XMLSchema}int" default="0" />
 *                         &lt;/restriction>
 *                       &lt;/complexContent>
 *                     &lt;/complexType>
 *                   &lt;/element>
 *                 &lt;/sequence>
 *               &lt;/restriction>
 *             &lt;/complexContent>
 *           &lt;/complexType>
 *         &lt;/element>
 *       &lt;/sequence>
 *     &lt;/extension>
 *   &lt;/complexContent>
 * &lt;/complexType>
 * </pre>
 * 
 * 
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "HalfedgeSelection", propOrder = {
    "vertexSelection",
    "edgeSelection",
    "faceSelection"
})
public class HalfedgeSelection
    extends ConformalData
{

    @XmlElement(name = "VertexSelection", required = true)
    protected HalfedgeSelection.VertexSelection vertexSelection;
    @XmlElement(name = "EdgeSelection", required = true)
    protected HalfedgeSelection.EdgeSelection edgeSelection;
    @XmlElement(name = "FaceSelection", required = true)
    protected HalfedgeSelection.FaceSelection faceSelection;

    /**
     * Ruft den Wert der vertexSelection-Eigenschaft ab.
     * 
     * @return
     *     possible object is
     *     {@link HalfedgeSelection.VertexSelection }
     *     
     */
    public HalfedgeSelection.VertexSelection getVertexSelection() {
        return vertexSelection;
    }

    /**
     * Legt den Wert der vertexSelection-Eigenschaft fest.
     * 
     * @param value
     *     allowed object is
     *     {@link HalfedgeSelection.VertexSelection }
     *     
     */
    public void setVertexSelection(HalfedgeSelection.VertexSelection value) {
        this.vertexSelection = value;
    }

    /**
     * Ruft den Wert der edgeSelection-Eigenschaft ab.
     * 
     * @return
     *     possible object is
     *     {@link HalfedgeSelection.EdgeSelection }
     *     
     */
    public HalfedgeSelection.EdgeSelection getEdgeSelection() {
        return edgeSelection;
    }

    /**
     * Legt den Wert der edgeSelection-Eigenschaft fest.
     * 
     * @param value
     *     allowed object is
     *     {@link HalfedgeSelection.EdgeSelection }
     *     
     */
    public void setEdgeSelection(HalfedgeSelection.EdgeSelection value) {
        this.edgeSelection = value;
    }

    /**
     * Ruft den Wert der faceSelection-Eigenschaft ab.
     * 
     * @return
     *     possible object is
     *     {@link HalfedgeSelection.FaceSelection }
     *     
     */
    public HalfedgeSelection.FaceSelection getFaceSelection() {
        return faceSelection;
    }

    /**
     * Legt den Wert der faceSelection-Eigenschaft fest.
     * 
     * @param value
     *     allowed object is
     *     {@link HalfedgeSelection.FaceSelection }
     *     
     */
    public void setFaceSelection(HalfedgeSelection.FaceSelection value) {
        this.faceSelection = value;
    }


    /**
     * <p>Java-Klasse für anonymous complex type.
     * 
     * <p>Das folgende Schemafragment gibt den erwarteten Content an, der in dieser Klasse enthalten ist.
     * 
     * <pre>
     * &lt;complexType>
     *   &lt;complexContent>
     *     &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
     *       &lt;sequence>
     *         &lt;element name="Edge" maxOccurs="unbounded" minOccurs="0">
     *           &lt;complexType>
     *             &lt;complexContent>
     *               &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
     *                 &lt;attribute name="index" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
     *                 &lt;attribute name="channel" type="{http://www.w3.org/2001/XMLSchema}int" default="0" />
     *               &lt;/restriction>
     *             &lt;/complexContent>
     *           &lt;/complexType>
     *         &lt;/element>
     *       &lt;/sequence>
     *     &lt;/restriction>
     *   &lt;/complexContent>
     * &lt;/complexType>
     * </pre>
     * 
     * 
     */
    @XmlAccessorType(XmlAccessType.FIELD)
    @XmlType(name = "", propOrder = {
        "edges"
    })
    public static class EdgeSelection {

        @XmlElement(name = "Edge")
        protected List<HalfedgeSelection.EdgeSelection.Edge> edges;

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
         * {@link HalfedgeSelection.EdgeSelection.Edge }
         * 
         * 
         */
        public List<HalfedgeSelection.EdgeSelection.Edge> getEdges() {
            if (edges == null) {
                edges = new ArrayList<HalfedgeSelection.EdgeSelection.Edge>();
            }
            return this.edges;
        }


        /**
         * <p>Java-Klasse für anonymous complex type.
         * 
         * <p>Das folgende Schemafragment gibt den erwarteten Content an, der in dieser Klasse enthalten ist.
         * 
         * <pre>
         * &lt;complexType>
         *   &lt;complexContent>
         *     &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
         *       &lt;attribute name="index" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
         *       &lt;attribute name="channel" type="{http://www.w3.org/2001/XMLSchema}int" default="0" />
         *     &lt;/restriction>
         *   &lt;/complexContent>
         * &lt;/complexType>
         * </pre>
         * 
         * 
         */
        @XmlAccessorType(XmlAccessType.FIELD)
        @XmlType(name = "")
        public static class Edge {

            @XmlAttribute(name = "index", required = true)
            protected int index;
            @XmlAttribute(name = "channel")
            protected Integer channel;

            /**
             * Ruft den Wert der index-Eigenschaft ab.
             * 
             */
            public int getIndex() {
                return index;
            }

            /**
             * Legt den Wert der index-Eigenschaft fest.
             * 
             */
            public void setIndex(int value) {
                this.index = value;
            }

            /**
             * Ruft den Wert der channel-Eigenschaft ab.
             * 
             * @return
             *     possible object is
             *     {@link Integer }
             *     
             */
            public int getChannel() {
                if (channel == null) {
                    return  0;
                } else {
                    return channel;
                }
            }

            /**
             * Legt den Wert der channel-Eigenschaft fest.
             * 
             * @param value
             *     allowed object is
             *     {@link Integer }
             *     
             */
            public void setChannel(Integer value) {
                this.channel = value;
            }

        }

    }


    /**
     * <p>Java-Klasse für anonymous complex type.
     * 
     * <p>Das folgende Schemafragment gibt den erwarteten Content an, der in dieser Klasse enthalten ist.
     * 
     * <pre>
     * &lt;complexType>
     *   &lt;complexContent>
     *     &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
     *       &lt;sequence>
     *         &lt;element name="Face" maxOccurs="unbounded" minOccurs="0">
     *           &lt;complexType>
     *             &lt;complexContent>
     *               &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
     *                 &lt;attribute name="index" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
     *                 &lt;attribute name="channel" type="{http://www.w3.org/2001/XMLSchema}int" default="0" />
     *               &lt;/restriction>
     *             &lt;/complexContent>
     *           &lt;/complexType>
     *         &lt;/element>
     *       &lt;/sequence>
     *     &lt;/restriction>
     *   &lt;/complexContent>
     * &lt;/complexType>
     * </pre>
     * 
     * 
     */
    @XmlAccessorType(XmlAccessType.FIELD)
    @XmlType(name = "", propOrder = {
        "faces"
    })
    public static class FaceSelection {

        @XmlElement(name = "Face")
        protected List<HalfedgeSelection.FaceSelection.Face> faces;

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
         * {@link HalfedgeSelection.FaceSelection.Face }
         * 
         * 
         */
        public List<HalfedgeSelection.FaceSelection.Face> getFaces() {
            if (faces == null) {
                faces = new ArrayList<HalfedgeSelection.FaceSelection.Face>();
            }
            return this.faces;
        }


        /**
         * <p>Java-Klasse für anonymous complex type.
         * 
         * <p>Das folgende Schemafragment gibt den erwarteten Content an, der in dieser Klasse enthalten ist.
         * 
         * <pre>
         * &lt;complexType>
         *   &lt;complexContent>
         *     &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
         *       &lt;attribute name="index" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
         *       &lt;attribute name="channel" type="{http://www.w3.org/2001/XMLSchema}int" default="0" />
         *     &lt;/restriction>
         *   &lt;/complexContent>
         * &lt;/complexType>
         * </pre>
         * 
         * 
         */
        @XmlAccessorType(XmlAccessType.FIELD)
        @XmlType(name = "")
        public static class Face {

            @XmlAttribute(name = "index", required = true)
            protected int index;
            @XmlAttribute(name = "channel")
            protected Integer channel;

            /**
             * Ruft den Wert der index-Eigenschaft ab.
             * 
             */
            public int getIndex() {
                return index;
            }

            /**
             * Legt den Wert der index-Eigenschaft fest.
             * 
             */
            public void setIndex(int value) {
                this.index = value;
            }

            /**
             * Ruft den Wert der channel-Eigenschaft ab.
             * 
             * @return
             *     possible object is
             *     {@link Integer }
             *     
             */
            public int getChannel() {
                if (channel == null) {
                    return  0;
                } else {
                    return channel;
                }
            }

            /**
             * Legt den Wert der channel-Eigenschaft fest.
             * 
             * @param value
             *     allowed object is
             *     {@link Integer }
             *     
             */
            public void setChannel(Integer value) {
                this.channel = value;
            }

        }

    }


    /**
     * <p>Java-Klasse für anonymous complex type.
     * 
     * <p>Das folgende Schemafragment gibt den erwarteten Content an, der in dieser Klasse enthalten ist.
     * 
     * <pre>
     * &lt;complexType>
     *   &lt;complexContent>
     *     &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
     *       &lt;sequence>
     *         &lt;element name="Vertex" maxOccurs="unbounded" minOccurs="0">
     *           &lt;complexType>
     *             &lt;complexContent>
     *               &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
     *                 &lt;attribute name="index" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
     *                 &lt;attribute name="channel" type="{http://www.w3.org/2001/XMLSchema}int" default="0" />
     *               &lt;/restriction>
     *             &lt;/complexContent>
     *           &lt;/complexType>
     *         &lt;/element>
     *       &lt;/sequence>
     *     &lt;/restriction>
     *   &lt;/complexContent>
     * &lt;/complexType>
     * </pre>
     * 
     * 
     */
    @XmlAccessorType(XmlAccessType.FIELD)
    @XmlType(name = "", propOrder = {
        "vertices"
    })
    public static class VertexSelection {

        @XmlElement(name = "Vertex")
        protected List<HalfedgeSelection.VertexSelection.Vertex> vertices;

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
         * {@link HalfedgeSelection.VertexSelection.Vertex }
         * 
         * 
         */
        public List<HalfedgeSelection.VertexSelection.Vertex> getVertices() {
            if (vertices == null) {
                vertices = new ArrayList<HalfedgeSelection.VertexSelection.Vertex>();
            }
            return this.vertices;
        }


        /**
         * <p>Java-Klasse für anonymous complex type.
         * 
         * <p>Das folgende Schemafragment gibt den erwarteten Content an, der in dieser Klasse enthalten ist.
         * 
         * <pre>
         * &lt;complexType>
         *   &lt;complexContent>
         *     &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
         *       &lt;attribute name="index" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
         *       &lt;attribute name="channel" type="{http://www.w3.org/2001/XMLSchema}int" default="0" />
         *     &lt;/restriction>
         *   &lt;/complexContent>
         * &lt;/complexType>
         * </pre>
         * 
         * 
         */
        @XmlAccessorType(XmlAccessType.FIELD)
        @XmlType(name = "")
        public static class Vertex {

            @XmlAttribute(name = "index", required = true)
            protected int index;
            @XmlAttribute(name = "channel")
            protected Integer channel;

            /**
             * Ruft den Wert der index-Eigenschaft ab.
             * 
             */
            public int getIndex() {
                return index;
            }

            /**
             * Legt den Wert der index-Eigenschaft fest.
             * 
             */
            public void setIndex(int value) {
                this.index = value;
            }

            /**
             * Ruft den Wert der channel-Eigenschaft ab.
             * 
             * @return
             *     possible object is
             *     {@link Integer }
             *     
             */
            public int getChannel() {
                if (channel == null) {
                    return  0;
                } else {
                    return channel;
                }
            }

            /**
             * Legt den Wert der channel-Eigenschaft fest.
             * 
             * @param value
             *     allowed object is
             *     {@link Integer }
             *     
             */
            public void setChannel(Integer value) {
                this.channel = value;
            }

        }

    }

}
