package linalg.PColt;

import linalg.EigenSolver;

import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import cern.colt.matrix.tdouble.algo.decomposition.DenseDoubleEigenvalueDecomposition;

/**
 * @author Luca Castelli Aleardi
 * 
 * Eigenvalue decomposition of a Laplacian matrix, based on Parallel Colt library (using sparse matrices)
 */
public class PColtEigenSolver implements EigenSolver {
	PColtMatrix laplacian; // laplacian matrix of the graph (based on Colt library)
	DenseDoubleEigenvalueDecomposition decomposition; // eigenvalue decomposition (based on Colt library)
	
	/**
	 * Create the Laplacian matrix of a graph (with n vertices)
	 * (not efficient implementation: it takes nxn time)
	 * 
	 * @param g
	 * 			the input graph
	 */	
	public PColtEigenSolver(PColtMatrix m) {
		this.laplacian=m;
	}

	public void computeEigenvalueDecomposition(int k) {
		System.out.print("Computing eigenvalue decomposition... ");
		long startTime=System.nanoTime(), endTime; // for evaluating time performances

		this.decomposition=new DenseDoubleEigenvalueDecomposition(this.laplacian.A); // Colt eigenvalue decomposition

		endTime=System.nanoTime();
        double duration=(double)(endTime-startTime)/1000000000.;
        System.out.println("done ("+duration+" seconds)");
	}
	
	public double[] getEigenvalues() {
		return this.decomposition.getRealEigenvalues().toArray();
	}
	
	public double[][] getEigenvectors() {
		DenseDoubleAlgebra algebra=new DenseDoubleAlgebra();
		return algebra.transpose(this.decomposition.getV()).toArray();
	}
	
	public String toString() {
		return this.laplacian.toString();
	}

}
