package linalg.JamaMatrices;

import linalg.Matrix;

/**
 * @author Luca Castelli Aleardi
 * 
 * Wrapper class for Jama matrices
 */
public class JamaMatrix implements Matrix {

	public Jama.Matrix A; // laplacian matrix of the graph (based on Jama library)
	
	/**
	 * Initialize the matrix
	 * 
	 * @param m
	 * 			the matrix (Jama Matrix)
	 */	
	public JamaMatrix(Jama.Matrix A) {
		this.A=A;
	}

	/**
	 * Create a Jama matrix from a double[][] array
	 * 
	 * @param m  the input array
	 */	
	public JamaMatrix(double[][] m) {
		this.A=new Jama.Matrix(m);
	}
	
	/**
	 * Computes the product of the matrix A and a given vector X
	 * 
	 * @param x
	 *            the input vector
	 * @return the product Ax
	 */
	public double[] times(double[] x) {
		double[][] X=new double[1][]; // row vector
		X[0]=x; // the transpose gives the column vector
		Jama.Matrix result=this.A.times(new Jama.Matrix(X).transpose()); // compute column vector Ax
		return A.times(result).transpose().getArray()[0]; // return the transpose of Ax
	}

	/**
	 * Return the transpose of the matrix
	 * 
	 * @return A'
	 */
	public JamaMatrix transpose() {
		Jama.Matrix result=this.A.transpose(); // compute the transpose of A
		return new JamaMatrix(result);
	}

	/**
	 * Return the product A times B
	 * 
	 * @return A x B
	 */
	public JamaMatrix times(JamaMatrix B) {
		Jama.Matrix result=this.A.times(B.A); // compute A x B
		return new JamaMatrix(result);
	}

	public int getHeight() {
		return this.A.getColumnDimension();
	}

	public int getWidth() {
		return this.A.getRowDimension();
	}

	public String toString() {
		return this.A.toString();
	}
	
}
