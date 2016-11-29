package linalg.PColt;

import java.util.Collection;

import linalg.Matrix;
import Jcg.graph.arraybased.ArrayBasedGraph;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix2D;

/**
 * @author Luca Castelli Aleardi
 * 
 * Wrapper class for sparse matrices (based on Parallel Colt matrix)
 */
public class PColtMatrix implements Matrix {

	public DoubleMatrix2D A; // (sparse) matrix implementation based on Parallel Colt library

	/**
	 * Create the Laplacian matrix of a graph (with n vertices)
	 * 
	 * @param g
	 * 			the input graph
	 */	
	public PColtMatrix(DoubleMatrix2D A) {
		this.A=A;
	}
	
	/**
	 * Create the Laplacian matrix of a graph (with n vertices)
	 * 
	 * @param g
	 * 			the input graph
	 */	
	public PColtMatrix(double[][] M) {
		int n=M.length; // matrix size
		System.out.print("Creating a sparse matrix from a double array"+n+" (using Parallel Colt library)...");
		
		this.A=new SparseDoubleMatrix2D(n, n); // create a sparse matrix of size nxn

		for(int i=0;i<n;i++) {
			for(int j=0;j<n;j++) {
				this.A.set(i, j, M[i][j]);
			}
	   	}
		DenseDoubleAlgebra doubleAlgebra=new DenseDoubleAlgebra();
		doubleAlgebra.inverse(this.A);
		System.out.println("done");
	}

	/**
	 * Computes the product of the matrix A and a given vector X
	 * 
	 * @param x
	 *            the input vector
	 * @return the product Ax
	 */
	public double[] times(double[] x) {
		int size=x.length;
		DoubleMatrix1D X = new DenseDoubleMatrix1D(size);
		DoubleMatrix1D Z = new DenseDoubleMatrix1D(size);
		for (int i=0; i<size; i++) {
			X.set(i, x[i]);
		}
		this.A.zMult(X, Z);
		double[] z=new double[size];
		for (int i=0; i<size; i++) {
			z[i] = Z.get(i);
		}
		
		return z;
	}
	
	/**
	 * Getter returning the height of the matrix
	 * @return the height of the matrix
	 */
	public int getHeight() {
		return this.A.rows();
	}

	/**
	 * Getter returning the height of the matrix
	 * @return the width of the matrix
	 */
	public int getWidth() {
		return this.A.columns();
	}
	
	public String toString() {
		return this.A.toString();
	}

}
