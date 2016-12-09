import java.util.ArrayList;
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
	
	public void compute() {
		CompRowMatrix M = new CompRowMatrix(polyhedron3D.vertices.size(),1, new int[polyhedron3D.vertices.size()][0]);
		M = M.zero();
		HashMap<Integer,Double> aires = new HashMap<>();
		ArrayList<Vertex<Point_3>> vertices = polyhedron3D.vertices;
		for (Vertex<Point_3> v : vertices) {
			ArrayList<Vertex<Point_3>> neighbours = new ArrayList<>();
			ArrayList<Point_3> circumcenters = new ArrayList<>();
			Halfedge<Point_3> baseHalfedge = v.getHalfedge();
			Halfedge<Point_3> currentHalfedge = v.getHalfedge();
			while (baseHalfedge != currentHalfedge || neighbours.size() == 0) {
				neighbours.add(currentHalfedge.getVertex());
				currentHalfedge = currentHalfedge.getOpposite().getNext();
			}
			for (int i = 0; i < neighbours.size(); i++) {
				//Compute the angle with the prev point
				Vertex<Point_3> current = neighbours.get(i);				
				Vertex<Point_3> prev = (Vertex<Point_3>) prev(i, neighbours);
				Vector_3 vecPrev1 = (Vector_3) prev.getPoint().minus(v.getPoint());
				Vector_3 vecPrev2 = (Vector_3) prev.getPoint().minus(current.getPoint());
				double cosAnglePrev = (double) vecPrev1.innerProduct(vecPrev2) / (Math.sqrt((double) vecPrev1.squaredLength() * (double) vecPrev2.squaredLength()));
				double sinAnglePrev = (double) vecPrev1.crossProduct(vecPrev2).squaredLength() / (Math.sqrt((double) vecPrev1.squaredLength() * (double) vecPrev2.squaredLength()));
				double cotAnglePrev = cosAnglePrev / sinAnglePrev;
				//Compute the angle with the next point
				Vertex<Point_3> next = (Vertex<Point_3>) next(i, neighbours);
				Vector_3 vecNext1 = (Vector_3) next.getPoint().minus(v.getPoint());
				Vector_3 vecNext2 = (Vector_3) next.getPoint().minus(current.getPoint());
				double cosAngleNext = (double) vecNext1.innerProduct(vecNext2) / (Math.sqrt((double) vecNext1.squaredLength() * (double) vecNext2.squaredLength()));
				double sinAngleNext = (double) vecNext1.crossProduct(vecNext2).squaredLength() / (Math.sqrt((double) vecNext1.squaredLength() * (double) vecNext2.squaredLength()));
				double cotAngleNext = cosAngleNext / sinAngleNext;
				M.set(v.index, current.index, (cotAnglePrev + cotAngleNext) / 2 );
				//Compute the circumcenter corresponding to the triangle with the prev point
				Vector_3 a = (Vector_3) v.getPoint().minus(prev.getPoint());
				Vector_3 b = (Vector_3) v.getPoint().minus(current.getPoint());
				Point_3 circumcenter = v.getPoint();
				circumcenter.translateOf(
								a.multiplyByScalar(b.squaredLength())
								.opposite()
								.sum(b.multiplyByScalar(a.squaredLength()))
								.crossProduct((a.crossProduct(b)))
						.multiplyByScalar(1/(2. * (double) a.crossProduct(b).squaredLength())));
				circumcenters.add(circumcenter);
			}
			double aire = 0;
			for (int i = 0; i < circumcenters.size(); i++) {
				double a = Math.sqrt((double) v.getPoint().minus(circumcenters.get(i)).squaredLength());
				double b = Math.sqrt((double) v.getPoint().minus((Point_3) next(i,circumcenters)).squaredLength());
				double c = Math.sqrt((double) circumcenters.get(i).minus((Point_3) next(i,circumcenters)).squaredLength());
				double p = 1/2*(a+b+c);
				aire += 1/4 * Math.sqrt(p*(p-a)*(p-b)*(p-c));
			}
			aires.put(v.index, aire);
		}
		
		CompRowMatrix L = new CompRowMatrix(polyhedron3D.vertices.size(), polyhedron3D.vertices.size(), null);
		for (int i = 0; i < polyhedron3D.vertices.size(); i++) {
			aires.get(i); //si
			M.get(i, j); //mi,j
			//for (int j )
			//L.set();
		}
		
		MTJSparseMatrix Laplacian = new MTJSparseMatrix(L);
		MTJSparseEigenSolver Solver = new MTJSparseEigenSolver(Laplacian);
		Solver.computeEigenvalueDecomposition(20); //le nombre de v propres qu'on veut
		double[] eigenvalues = Solver.getEigenvalues();
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
