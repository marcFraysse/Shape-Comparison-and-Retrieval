import processing.core.*;

import Jcg.geometry.*;
	

public class MeshViewer extends PApplet {

	SurfaceMesh mesh;
	//String filename="OFF/high_genus.off";
	//String filename="OFF/sphere.off";
	//String filename="OFF/cube.off";
	String filename="OFF/MeshsegBenchmark-1.0/data/off/22.off";
	//String filename="OFF/tore.off";
	//String filename="OFF/tri_round_cube.off";
	//String filename="OFF/tri_hedra.off";
	//String filename="OFF/tri_horse.off";
	
	public double comparison(String filename1, String filename2){
		SurfaceMesh mesh1 = new SurfaceMesh(this, filename1);
		SurfaceMesh mesh2 = new SurfaceMesh(this, filename2);
		double[] histo1 = mesh1.compute();
		double[] histo2 = mesh2.compute();
		
		double dist = 0;
		for(int i = 0; i<histo1.length;i++){
			dist += Math.abs(histo1[i]-histo2[i]);
		}
		return dist;
	}
	

	public void setup() {
		  size(800,600,P3D);
		  ArcBall arcball = new ArcBall(this);
		  
		  this.comparison("OFF/MeshsegBenchmark-1.0/data/off/101.off","OFF/MeshsegBenchmark-1.0/data/off/102.off");

	}

		
		/**
		 * For running the PApplet as Java application
		 */
		public static void main(String args[]) {
			//PApplet pa=new MeshViewer();
			//pa.setSize(400, 400);
			PApplet.main(new String[] { "MeshViewer" });
		}
}
