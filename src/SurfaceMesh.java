
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;

import Jcg.geometry.*;
import Jcg.mesh.SharedVertexRepresentation;
import Jcg.polyhedron.*;
import linalg.MTJ.*;
import no.uib.cipr.matrix.sparse.CompRowMatrix;


/**
 * Class for rendering a surface triangle mesh (using Processing)
 * 
 * @author Luca Castelli Aleardi (INF555, 2012)
 *
 */
public class SurfaceMesh {
	
	double scaleFactor=60; // scaling factor: useful for 3d rendering
	MeshViewer view; // Processing 3d frame (where meshes are rendered)
	public Polyhedron_3<Point_3> polyhedron3D; // triangle mesh
	
	/**
	 * Create a surface mesh from an OFF file
	 */	
	public SurfaceMesh(MeshViewer view, String filename) {
		this.view=view;

		// shared vertex representation of the mesh
    	SharedVertexRepresentation sharedVertex=new SharedVertexRepresentation(filename);
    	LoadMesh<Point_3> load3D=new LoadMesh<Point_3>();
    	
    	polyhedron3D=load3D.createTriangleMesh(sharedVertex.points,sharedVertex.faceDegrees,
				sharedVertex.faces,sharedVertex.sizeHalfedges);

    	//System.out.println(polyhedron3D.verticesToString());   	
    	//System.out.println(polyhedron3D.facesToString());
    	polyhedron3D.isValid(false);
    	    	
    	this.scaleFactor=this.computeScaleFactor();
	}

	
	
	/* this function computes the eigenvalues and eigenvectors linked to a shape and then creates a distance histogram according to the method described in the rustamov article*/
	public double[] compute() {

		ArrayList<Vertex<Point_3>> vertices = polyhedron3D.vertices;
		int count = 0;
		HashMap<Vertex<Point_3>,ArrayList<Vertex<Point_3>>> neighboursHash = new HashMap<>();

		int nbEigen = 25;
		int n = polyhedron3D.vertices.size();
		CompRowMatrix M;
		int[][] m = new int[n][]; // indices of non-zero entries for each row
		HashMap<Integer,Double> aires = new HashMap<>();
		for (Vertex<Point_3> v : vertices) {
			
			
			ArrayList<Vertex<Point_3>> neighbours = new ArrayList<>();
			Halfedge<Point_3> baseHalfedge = v.getHalfedge();
			Halfedge<Point_3> currentHalfedge = v.getHalfedge();
			while (baseHalfedge != currentHalfedge || neighbours.size() == 0) {
				neighbours.add(currentHalfedge.getVertex());
				currentHalfedge = currentHalfedge.getOpposite().getNext();
			}
			m[count] = new int[neighbours.size()+1];
			m[count][neighbours.size()]=count;
			for(int j = 0; j < neighbours.size(); j++) {
				m[count][j]=neighbours.get(j).index; //il faut bien mettre les indices pour lesquels la matrice sera non nuls ici (Les valeurs non nulles de la matrices sont toutes celles dont les indices sont (count, index))
			}
			neighboursHash.put(v, neighbours);
			
			
			count++;
		}
		M = new CompRowMatrix(n,n,m);
		
		for (Vertex<Point_3> v : vertices) { 

			ArrayList<Vertex<Point_3>> neighbours = neighboursHash.get(v);

			double aire =0;
			for (int i = 0; i < neighbours.size(); i++) {
				
				Vertex<Point_3> current = neighbours.get(i);				
				Vertex<Point_3> prev = (Vertex<Point_3>) prev(i, neighbours);
				Vertex<Point_3> next = (Vertex<Point_3>) next(i, neighbours);

				
				double d02 = (double)( v.getPoint().minus(current.getPoint()).squaredLength());
				double d12 = (double)( v.getPoint().minus(prev.getPoint()).squaredLength());
				double d22 = (double)( prev.getPoint().minus(current.getPoint()).squaredLength());
				double d32 = (double)( v.getPoint().minus(next.getPoint()).squaredLength());
				double d42 = (double)( next.getPoint().minus(current.getPoint()).squaredLength());
				
				double pA1 = 0.5*(Math.sqrt(d02)+Math.sqrt(d12)+Math.sqrt(d22));
				
				double A1  = 0.25 * Math.sqrt(pA1*(pA1-d02)*(pA1-d12)*(pA1-d22));
				
				double pA2 = 0.5*(Math.sqrt(d02)+Math.sqrt(d32)+Math.sqrt(d42));
				
				double A2  = 0.25 * Math.sqrt(pA2*(pA2-d02)*(pA2-d32)*(pA2-d42));
				
				double wij = 1/(8*A1)*(d02-d12-d22)+1/(8*A2)*(d02-d32-d42);
				
				M.set(v.index, current.index, wij);
				aire += A1;

							
			}
			aires.put(v.index, aire/3);
			
		}

		CompRowMatrix L = new CompRowMatrix(n,n, m);
		for (int i = 0; i < n; i++) {
			double si = aires.get(i);
			L.set(i,i,0);
			for (int k = 0 ; k<m[i].length;k++){
				int j = m[i][k];
				double mij = M.get(i, j);
				if(j==i){
					;
				}
				else{
					L.set(i,j,-mij/si);
					L.set(i,i,L.get(i,i)+mij/si);
				}
			}
		}
		
		
		MTJSparseMatrix Laplacian = new MTJSparseMatrix(L);
		MTJSparseEigenSolver Solver = new MTJSparseEigenSolver(Laplacian);
		Solver.computeEigenvalueDecomposition(nbEigen); 
		double[] eigenvalues = Solver.getEigenvalues();
		double[][] eigenvectors = Solver.getEigenvectors();
		
		for(int i = 0; i<eigenvalues.length;i++){
			System.out.println(eigenvalues[i]);
		}
		
		
		//normalise eigenvectors (S-inner product is norm)
		
		for (int i=0; i<eigenvectors.length;i++){
			double norm = 0;
			for(int j = 0; j < eigenvectors[i].length;j++){
				norm += eigenvectors[i][j]*eigenvectors[i][j]*aires.get(j);
			}
			for(int j = 0; j < eigenvectors[i].length;j++){
				eigenvectors[i][j]/=norm;
			}
		}

		
		//get coordinates of each points
		
		double[][] coordinates = new double[n][nbEigen]; //second dimension is nb of eigenvectors we have
		for(int i = 0; i<n; i++){
			for(int j = 0 ; j<nbEigen; j++){
				coordinates[i][j]=eigenvectors[j][i]/Math.sqrt(eigenvalues[j]); //this formula is supposed to be the discrete version of the eigenfunction(point)*1/sqrt(eigenvalue) continuous formula
			}
		}
		
		// get histogram
		// we should get subsample, no same amout of points each time
		
		double[] histogram = new double[100]; // we'll need to know the max dist between two points 
		double maxdist = 10;
		for(int i = 0;i<n;i++){
			for(int j = 0;j<n;j++){
				histogram[(int)(dist(coordinates[i],coordinates[j])*100/maxdist)]+=1;
			}
		}
		
		//compute returns histo and then MeshViewer compares several histos
		return histogram;
	}
	
	
	public Object next(int pos, ArrayList<?> neighbours) {
		if (pos == neighbours.size() - 1)
			return neighbours.get(0);
		else
			return neighbours.get(pos + 1);
	}
	
	public Object prev(int pos, ArrayList<?> neighbours) {
		if (pos == 0)
			return neighbours.get(neighbours.size() - 1);
		else
			return neighbours.get(pos - 1);
	}
	
	public double dist(double[] coordinatesA, double[] coordinatesB){
		double dist = 0;
		for(int i = 0; i < coordinatesA.length;i++){
			dist += Math.pow(coordinatesA[i]-coordinatesB[i],2);
		}
		dist = Math.sqrt(dist);
		return dist;
	}
	
	
	/**
	 * Draw a segment between two points
	 */	
	public void drawSegment(Point_3 p, Point_3 q) {
		float s=(float)this.scaleFactor;
		float x1=(float)p.getX().doubleValue()*s;
		float y1=(float)p.getY().doubleValue()*s;
		float z1=(float)p.getZ().doubleValue()*s;
		float x2=(float)q.getX().doubleValue()*s;
		float y2=(float)q.getY().doubleValue()*s;
		float z2=(float)q.getZ().doubleValue()*s;
		this.view.line(	x1, y1, z1, x2, y2, z2 );		
	}
	
	/**
	 * Draw a vertex (as a small sphere)
	 */	
	public void drawVertex(Point_3 p) {
		float s=(float)this.scaleFactor;
		float x1=(float)p.getX().doubleValue()*s;
		float y1=(float)p.getY().doubleValue()*s;
		float z1=(float)p.getZ().doubleValue()*s;
		
		view.translate(x1, y1, z1);
		view.sphere(s/25f);
		view.translate(-x1, -y1, -z1);
	}


	/**
	 * Draw a triangle
	 */	
	public void drawTriangle(Point_3 p, Point_3 q, Point_3 r) {
		float s=(float)this.scaleFactor;
		view.vertex( (float)(p.getX().doubleValue()*s), (float)(p.getY().doubleValue()*s), (float)(p.getZ().doubleValue()*s));
		view.vertex( (float)(q.getX().doubleValue()*s), (float)(q.getY().doubleValue()*s), (float)(q.getZ().doubleValue()*s));
		view.vertex( (float)(r.getX().doubleValue()*s), (float)(r.getY().doubleValue()*s), (float)(r.getZ().doubleValue()*s));
	}

	/**
	 * Draw a (triangle or polygonal) face
	 */	
	public void drawFace(Face<Point_3> f) {
		Halfedge<Point_3> h=f.getEdge();
		Halfedge<Point_3> pEdge=h.getNext();
		
		Point_3 u=h.getOpposite().getVertex().getPoint();
		view.noStroke();
		view.fill(200,200,200,255); // color of the triangle
		
		while(pEdge.getVertex()!=h.getOpposite().getVertex()) {
			Point_3 v=pEdge.getOpposite().getVertex().getPoint();
			Point_3 w=pEdge.getVertex().getPoint();
			this.drawTriangle(u, v, w); // draw a triangle face
			
			pEdge=pEdge.getNext();
		}
	}


	/**
	 * Draw the entire mesh
	 */
	public void draw(int type) {
		//this.drawAxis();
		
		// draw all faces
		view.beginShape(view.TRIANGLES);
		for(Face<Point_3> f: this.polyhedron3D.facets) {
				this.drawFace(f);
		}
		view.endShape();
		
		if(type==1) return; // no rendering of edges
		
		// draw all edges
		view.strokeWeight(2); // line width (for edges)
		view.stroke(20);
		for(Halfedge<Point_3> e: this.polyhedron3D.halfedges) {
			Point_3 p=e.vertex.getPoint();
			Point_3 q=e.opposite.vertex.getPoint();
			
			this.drawSegment(p, q); // draw edge (p,q)
		}
		//view.strokeWeight(1);
		
		if(type==0) return; // no rendering for vertices
		
		view.noStroke();
		view.fill(0f, 0f, 250f);
		for(Vertex<Point_3> v: this.polyhedron3D.vertices) {
			this.drawVertex(v.getPoint());
		}
		view.strokeWeight(1);
	}
	
	/**
	 * Draw the X, Y and Z axis
	 */
	public void drawAxis() {
		double s=1;
		Point_3 p000=new Point_3(0., 0., 0.);
		Point_3 p100=new Point_3(s, 0., 0.);
		Point_3 p010=new Point_3(0.,s, 0.);
		Point_3 p011=new Point_3(0., 0., s);
		
		drawSegment(p000, p100);
		drawSegment(p000, p010);
		drawSegment(p000, p011);
	}


	/**
	 * Return the value after truncation
	 */
	public static double round(double x, int precision) {
		return ((int)(x*precision)/(double)precision);
	}
	
	/**
	 * Compute the scale factor (depending on the max distance of the point set)
	 */
	public double computeScaleFactor() {
		if(this.polyhedron3D==null || this.polyhedron3D.vertices.size()<1)
			return 1;
		double maxDistance=0.;
		Point_3 origin=new Point_3(0., 0., 0.);
		for(Vertex<Point_3> v: this.polyhedron3D.vertices) {
			double distance=Math.sqrt(v.getPoint().squareDistance(origin).doubleValue());
			maxDistance=Math.max(maxDistance, distance);
		}
		return Math.sqrt(3)/maxDistance*150;
	}
	
	/**
	 * Update the scale factor
	 */
	public void updateScaleFactor() {
		this.scaleFactor=this.computeScaleFactor();
	}

	
}
