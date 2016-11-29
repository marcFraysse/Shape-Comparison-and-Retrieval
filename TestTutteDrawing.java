import Jcg.geometry.Point_2;
import Jcg.graph.arraybased.ArrayBasedGraph;
import Jcg.graph.arraybased.ArrayBasedGraphLoader;
import Jcg.graph.arraybased.drawing.ArrayBasedGraphDrawing;

/**
 * This class provides input/output methods for testing (planar) drawing algorithms
 */
public class TestTutteDrawing {

	/**
	 * Test graph drawing algorithms
	 */
    public static void main (String[] args) throws InterruptedException {
		System.out.println("Testing Tutte drawing");
		if(args.length<4) {
			System.out.println("Error: wrong input parameters (at least four parameters required)");
			System.out.println("\tfirst parameter: input file");
			System.out.println("\tother parameter: v1 v2 v3 ... (integer values: exterior vertices)");
			System.exit(0);
		}
		
    	String filename=args[0];
    	int k=args.length-1; // number of exterior vertices
    	int[] exteriorVertices=new int[k];
    	System.out.println("Outer face: ");
    	for(int i=0;i<k;i++) {
    		exteriorVertices[i]=Integer.parseInt(args[i+1]);
    		System.out.print(" "+exteriorVertices[i]);
    	}
    	System.out.println();
 
    	ArrayBasedGraph g=null; // the input graph
 		g=ArrayBasedGraphLoader.readGraphFromFile(filename); // load graph from file
		ArrayBasedGraphDrawing<Point_2> d=new ArrayBasedTutteDrawing<Point_2>(g, exteriorVertices);

    	d.computeDrawing();
    	System.out.println("planar representation computed");
    	Thread.sleep(1000);
    	d.draw2D();
    		
    	//ArrayBasedGraphLoader.writeMeshSkeletonToFile("Data/tri_round_cube.off", "round_cube.txt");
    }

}
