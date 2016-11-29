import java.util.*;

import Jama.Matrix;
import Jcg.geometry.*;
import Jcg.graph.arraybased.*;
import Jcg.graph.arraybased.drawing.*;

/**
 * Provides methods for drawing graphs in 2D using Tutte barycentric method.
 * The used data structure is an array based representation of the graph
 *
 * @author Luca Castelli Aleardi (Ecole Polytechnique)
 * @version 2014
 */
public class ArrayBasedTutteDrawing<X extends Point_> extends ArrayBasedGraphDrawing<X> {
	Matrix M; // linear system to solve (Tutte laplacian restricted to inner vertices)

	int[] exteriorVertices; // indices of vertices lying on the boundary cycle: the order counts here
	Point_[] exteriorPoints; // geometric coordinates of the boundary vertices
	boolean[] isOnBoundary; // says whether a vertex is on the boundary cycle

	public ArrayBasedTutteDrawing() {}

	/**
	 * Initialize the Tutte method for barycentric drawing
	 * 
	 * @param g the input graph to draw
	 * @param exteriorVertices array storing the indices of the vertices on the boundary cycle
	 */	
	public ArrayBasedTutteDrawing(ArrayBasedGraph g, int[] exteriorVertices) {
		this.g=g; // the input graph
    	this.points=new ArrayList<X>(g.sizeVertices()); // create an array for storing vertex coordinates (to be computed)
    	this.exteriorVertices=exteriorVertices;
    	int sizeOuterFace=exteriorVertices.length;
    	
		double radius=5.;
		// place exterior vertices on a regular polygon
    	this.exteriorPoints=regularPolygonVertices(sizeOuterFace, radius);
    	
    	this.isOnBoundary=new boolean[g.sizeVertices()];
    	for(int i=0;i<exteriorVertices.length;i++)
    		isOnBoundary[exteriorVertices[i]]=true;
	}
	
	/**
	 * Compute and return the vector b, describing the constraints of vertices adjacent
	 * to the outer cycle (the system to solve is Mx=b)
	 * 
	 * @return the vector b (having size n)
	 */	
	public double[] computeVectorB(ArrayBasedGraph g) {
		int n=g.sizeVertices();
		double[] b=new double[n]; // the result
		
		throw new Error("To be completed");
		// return b;
	}
	
	/**
	 * Extract a sub-matrix from M, removing k rows and columns
	 * 
	 * @param M
	 * 			the original matrix of size nxn
	 * @param indices
	 * 				an array storing the k indices of rows/columns to delete
	 * 
	 * @return the (n-k)x(n-k) sub-matrix obtained from M, removing k rows and columns
	 */	
	public Matrix subMatrix(Matrix m, int[] indices) {
		int k=indices.length; // number of rows/columns to delete
		double[][] result=new double[this.g.sizeVertices()-k][this.g.sizeVertices()-k];
		
		throw new Error("To be completed");
	   	//return new Matrix(result);
	}
	
	/**
	 * Extract a sub-vector from b, removing k entries
	 * 
	 * @param b
	 * 			the original vector of size n
	 * @param indices
	 * 				an array storing the k indices of entries to delete
	 * 
	 * @return the (n-k) sub-vector obtained from b, removing k entries
	 */	
	public Matrix subVector(double[] b, int[] indices) {
		int k=indices.length; // number of rows/columns to delete
		double[][] result=new double[this.g.sizeVertices()-k][this.g.sizeVertices()-k];
		
		throw new Error("To be completed");
	   	//return new Matrix(result);
	}


	/**
	 * compute the Tutte drawing of a planar graph iteratively
	 * using the Force-Directed paradigm
	 */	
	public void computeDrawing() {
		if(this.exteriorPoints.length<3) 
			throw new Error("error: too few exterior points");
		
		// initialize vertex locations
		int count=0; // counter for exterior  vertices
		for(int i=0;i<g.sizeVertices();i++) {
			if(isOnBoundary[i]==true) { // vertex v_i is on the outer boundary cycle
				this.points.add((X)this.exteriorPoints[count]);
				count++;
			}
			else { // vertex v_i is an inner vertex
				Point_2 p=new Point_2();
				p.setOrigin(); // initialise the point at the origin
				this.points.add((X)p);
			}
		}
		
		// computes planar coordinates solving a system of linear equations
	}
	
	/**
	 * return the vertices of a regular polygon (equally spaced on a circle of radius r)
	 * 
	 * @param r the radius of the circle
	 * @param n the number of points on the outer cycle
	 * @return Point_2[] an array of 2D points, storing the vertices of a regular polygon
	 */	
	public static Point_2[] regularPolygonVertices(int n, double r) {
		Point_2[] vertices=new Point_2[n];
		double x,y;
		
		for(int i=0;i<n;i++) {
			x=r*Math.cos((2.*Math.PI/n)*i);
			y=r*Math.sin((2.*Math.PI/n)*i);
			vertices[i]=new Point_2(x,y);
		}
		return vertices;
	}

}
