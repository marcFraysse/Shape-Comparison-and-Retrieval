package linalg.JamaMatrices;

import linalg.LinearSolverInterface;
import linalg.Matrix;

/**
 * @author Luca Castelli Aleardi
 * 
 *         Eigenvalue decomposition of a Laplacian matrix, based on Jama library (using dense matrices)
 */
public class JamaLinearSolver implements LinearSolverInterface {

	JamaMatrix m; // matrix (based on Jama library)

	/**
	 * Initialize the linear system solver (direct solver based on Jama library)
	 * 
	 * @param A
	 *            the input matrix
	 */
	public JamaLinearSolver(JamaMatrix A) {
		this.m=A;
	}

	/**
	 * Solve linear system Lx=b (where L is the graph laplacian)
	 * 
	 * @param b
	 *            right hand side vector
	 * 
	 * @return the vector solution x[]
	 */
	public double[] solve(double[] b) {
		double[][] B=new double[1][]; // row vector
		B[0]=b; // the transpose gives the column vector
		Jama.Matrix X=this.m.A.solve(new Jama.Matrix(B).transpose()); // compute column vector Ax
		return X.transpose().getArray()[0]; // return the transpose of Ax
	}

	/**
	 * Return the product Ax
	 * */
	public double[] times(double[] x) {
		return this.m.times(x);
	}
	
	/**
	 * Getter returning the height of the matrix
	 * 
	 * @return the height of the matrix
	 */
	public int getHeight() {
		return this.m.getHeight();
	}

	/**
	 * Return the input matrix (representing the linear system to solve)
	 */
	public Matrix getMatrix() {
		return this.m;
	}

	public String toString() {
		return this.m.toString();
	}

	/**
	 * Check whether Ax=b
	 * 
	 * @param b
	 *            right side vector
	 * @param x
	 *            solution to check
	 * @param prec
	 *            numeric tolerance
	 * 
	 * @return true if dist(Ax, b) < precision
	 **/
	public boolean checkSolution(double[] x, double[] b) {
		throw new Error("To be completed");
	}

}
