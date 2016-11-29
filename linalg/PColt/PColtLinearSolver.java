package linalg.PColt;

import Jcg.graph.arraybased.ArrayBasedGraph;
import cern.colt.matrix.Norm;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.algo.solver.DefaultDoubleIterationMonitor;
import cern.colt.matrix.tdouble.algo.solver.DoubleCG;
import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;
import linalg.LinearIterativeSolver;
import linalg.Matrix;

/**
 * @author Luca Castelli Aleardi
 * 
 * Conjugate Gradient solver, based on Parallel Colt library (using sparse matrices)
 */
public class PColtLinearSolver implements LinearIterativeSolver {

	PColtMatrix laplacian; // laplacian matrix of the graph (based on Parallel Colt library)
	double[] start; // starting vector (initial guess)
	double precision; // numeric tolerance
	DoubleCG itSolver; // iterative solver implemented in Parallel Colt (Conjugate Gradient)
	DefaultDoubleIterationMonitor m;
	
	private int verbosity=NO_STATISTICS; // decide whether to show convergence statistics

	/**
	 * Initialize the Conjugate Gradient solver of Parallel Colt
	 * 
	 * @param A
	 * 			the input matrix
	 */	
	public PColtLinearSolver(PColtMatrix A, double precision) {
		int n=A.getHeight(); // matrix size
		
		this.start=new double[n]; // initial guess for the CG iterator
		this.laplacian=A; // create a sparse matrix of size nxn
		this.precision=precision;
		
		this.m=(DefaultDoubleIterationMonitor)itSolver.getIterationMonitor();
		m.setMaxIterations(10000);
		m.setNormType(Norm.Two); // choose euclidean norm
		this.m.setRelativeTolerance(precision);
	}
	
	/**
	 * Solve linear system Lx=b (where L is the graph laplacian)
	 * 
	 * @param b right hand side vector
	 * 
	 * @return the vector solution x[]
	 */
	public double[] solve(double[] b) {
		int size=b.length;
		double[] x=new double[size]; // solution (and initial guess)
		
		//DoubleDiagonal preconditioner=new DoubleDiagonal(size);
		//preconditioner.setMatrix(this.laplacian);
		//itSolver.setPreconditioner(preconditioner);
		
		DoubleMatrix1D X = new DenseDoubleMatrix1D(size); // solution
		DoubleMatrix1D B = new DenseDoubleMatrix1D(size);
		
		for (int i=0; i<size; i++) {
			B.set(i, b[i]);
			X.set(i, start[i]); // set initial guess
		}
		try {
			itSolver.solve(this.laplacian.A, B, X); // compute solution of the linear system
		} catch (IterativeSolverDoubleNotConvergedException e) {
			e.printStackTrace();
		}
		for (int i=0; i<size; i++) {
			x[i] = X.get(i);
		}
		
		if(this.verbosity>0) {
			int iter=m.iterations();
			double residual=m.residual();
			System.err.println("  Conjugate Gradient solved in "+iter+" turns, dist = "+residual);
		}

		return x;
	}
	
	/** 
	 * Return the product Ax
	 * */
	public double[] times(double[] x) {
		return this.laplacian.times(x);
	}
	
	/** 
	 * Run the iterative solver k times (solve the system Ax=b)
	 * 
	 * @param b
	 * 			right side vector
	 * @param k
	 * 			number of iterations to perform
	 * 
	 * @return	the solution of Ax=b
	 * 
	 * */
	public double[] runKtimes(double[] b,int k) {
		throw new Error("To be implemented");
	}
	
	/**
	 * Set the starting vector (initial guess for the CG method)
	 * 
	 * @param x starting vector
	 */
	public void setStartingVector(double[] start) {
		this.start=start.clone();
	}
	
	public void setPrecision(double precision){
		this.precision = precision;
		this.m.setRelativeTolerance(precision);
	}
	
	public double getPrecision(){
		return this.precision;
	}
	
	/**
	 * Getter returning the height of the matrix
	 * @return the height of the matrix
	 */
	public int getHeight() {
		return this.laplacian.getHeight();
	}
	
	/**
	 * Return the input matrix (representing the linear system to solve)
	 */
	public Matrix getMatrix() {
		return this.laplacian;
	}
	
	public void setVerbosity(int verbosity){
		this.verbosity=verbosity;
	}
	
	public String toString() {
		return this.laplacian.toString();
	}
	
	/** 
	 * Check whether Ax=b
	 * 
	 *  @param b
	 *  		right side vector
	 *  @param x
	 *  		solution to check
	 *  @param prec
	 *  		numeric tolerance
	 *  
	 *  @return true if dist(Ax, b) < precision
	 **/
	public boolean checkSolution(double[] x, double[] b, double prec) {
		throw new Error("To be completed");
	}

}
