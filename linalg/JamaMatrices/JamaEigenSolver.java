package linalg.JamaMatrices;

import linalg.EigenSolver;

/**
 * @author Luca Castelli Aleardi
 * 
 * Eigenvalue decomposition of a Laplacian matrix, based on Jama library
 */
public class JamaEigenSolver implements EigenSolver {

	public JamaMatrix laplacian; // laplacian matrix of the graph (based on Jama library)
	Jama.EigenvalueDecomposition decomposition; // eigenvalue decomposition (based on Jama library)
	
	private int k; // not defined at the beginning
	private double[][] firstEigenvectors; // array storing the first k eigenvectors: the first is stored in firstEigenvectors[0]
	private double[] firstEigenvalues; // array storing the first k eigenvectors: the first is stored in firstEigenvectors[0]
	
	/**
	 * Create the Laplacian matrix of a graph
	 * 
	 * @param m  the laplacian matrix (Jama Matrix implementation)
	 */	
	public JamaEigenSolver(JamaMatrix m) {
		this.laplacian=m;
	}

	public int getHeight() {
		return this.laplacian.getHeight();
	}

	public int getWidth() {
		return this.laplacian.getWidth();
	}

	/**
	 * Compute and store the first k smallest eigenvectors/eigenvalues of the matrix
	 * 
	 */	
	public void computeEigenvalueDecomposition(int k) {
		this.k=k;
		System.out.print("Computing eigenvalue decomposition... ");
		long startTime=System.nanoTime(), endTime; // for evaluating time performances

		this.decomposition=this.laplacian.A.eig(); // Jama eigenvalue decomposition
		
		this.firstEigenvectors=new double[k][this.getHeight()]; // first k eigenvectors
		double[][] m=this.decomposition.getV().transpose().getArray(); // array storing all n eigenvectors
		for(int i=0;i<k;i++)
			this.firstEigenvectors[i]=m[i]; // store only first k eigenvectors
		
		this.firstEigenvalues=new double[k]; // first k eigenvalues
		double[] v=this.decomposition.getRealEigenvalues(); // array storing all n eigenvalues
		for(int i=0;i<k;i++)
			this.firstEigenvalues[i]=v[i]; // store only first k eigenvectors
		
		endTime=System.nanoTime();
        double duration=(double)(endTime-startTime)/1000000000.;
		System.out.println("done ("+duration+" seconds)");
	}
	
	/**
	 * Return the already computed (real) eigenvalues
	 */	
	public double[] getEigenvalues() {
		return this.firstEigenvalues;
	}
	
	/**
	 * Return the already computed eigenvectors
	 */	
	public double[][] getEigenvectors() {
		return this.firstEigenvectors;
	}
	
	public String toString() {
		return this.laplacian.toString();
	}
	
}
