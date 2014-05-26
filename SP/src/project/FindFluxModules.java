package project;

import java.io.FileNotFoundException;
import java.util.ArrayList;

import lpsolve.LpSolveException;

import org.sbml.jsbml.*;

import scpsolver.constraints.LinearConstraint;
import scpsolver.constraints.LinearEqualsConstraint;
import scpsolver.lpsolver.LinearProgramSolver;
import scpsolver.lpsolver.SolverFactory;
import scpsolver.problems.LinearProgram;
import scpsolver.util.SparseMatrix;



/**<code>FindFluxModules</code> calculates the flux modules of a metabolic network.
 * <p> 
 * Important!!! The class and all its methods are needing the information of the class SBMLLoad(), represented by {@link #load load}!<p> 
 * 
 * @author guru
 *
 */

public class FindFluxModules{

	/**Initialized class SBMLLoad. Commited by the constructor.
	 */
	SBMLLoad load;
	/**
	 * Saves the maximized value of the biomass reaction (optimized by the method optimize()).
	 */
	double biomassOptValue = 0;
	/**
	 * Is a list of lists of the Id's of the reactions of the modules, at which each inner list is one module.
	 */
	ArrayList<ArrayList<Integer>> adjacency = new ArrayList<ArrayList<Integer>>();
	/**
	 * Number of variable Reactions. Is set by {@link #computeMinModules computeMinModules}
	 */
	int numV = 0;
	/**
	 * If true, the error output is redirected to a file with the choosen path. <p>
	 * Can be set true as a parameter with the start of the program.
	 */
	public boolean ignore;

	//getter & setter
	/**biomassOptValue is initialize by {@link #optimize(LinearProgramSolver solver, SparseMatrix rctMetArr) optimize}. Before that, biomassOptValue is null.
	 * @return {@link #biomassOptValue biomassOptValue}
	 */
	public double getBiomassOptValue() {
		return biomassOptValue;
	}
	/**adjacency is initialize by {@link #dfs(boolean[] isVariable, SparseMatrix dfsArray) dfs}. Before that, adjacency is null.
	 * @return {@link #adjacency adjacency}
	 */
	public ArrayList<ArrayList<Integer>> getAdjacency() {
		return adjacency;
	}
	
	//constructor
	public FindFluxModules(SBMLLoad load, boolean ignore) {
		this.load = load;
		this.ignore = ignore;
	}

	/**<code>findFlux()</code> calls all methods in <code>FindFluxModules.java</code> and passes the arguments between them.
	 * <p> 
	 * Needs correctly initalised {@link #load load}.<p>
	 * @throws LpSolveException
	 * @throws FileNotFoundException
	 */
	public void findFlux() throws LpSolveException, FileNotFoundException {
		
		// build matrix
		System.out.println("Creating Reaction/Metabolites Matrix");
		SparseMatrix rctMetArr = matrixBuild();

		// find optimum using LP
		System.out.println("Loading LP Solver");
		LinearProgramSolver solver = SolverFactory.newDefault();

		System.out.println("Solving LP for biomassreaction");
		double[] vBiomass = optimize(solver, rctMetArr);

		// find min max values
		System.out.println("Finding set V:");
		boolean[] isVariable = minMax(vBiomass, solver, rctMetArr);

		System.out.println("\nFinding minimal modules:");
		SparseMatrix dfsArray = computeMinModules(isVariable, rctMetArr, solver);
		
		System.out.println("Depth-First Search");
		dfs(isVariable,dfsArray);
		System.out.println("printing modules:");
		printModules();
	}
	/**<code>printModules()</code> print the ID's of the modules to the standard output.
	 * <p> 
	 * The modules are separated by lines.<p>
	 *  It needs field {@link #adjacency adjacency} and correctly initalised field {@link #load load}.<p>
	 */
	public void printModules() {
		for (int i = 0; i < adjacency.size(); i++) {
			ArrayList<Integer> t = adjacency.get(i);
			int c = t.size();
			pr("\nModul number " + i + ":\t");
			for (int w = 0; w < c-1; w++) {
				Reaction r = load.getModel().getReaction(t.get(w));
				String id = r.getId();
				System.out.print(id + ", ");
			}
			Reaction r = load.getModel().getReaction(t.get(c-1));
			String id = r.getId();
			System.out.print(id + "\n");
		}
	} 
	
	
	/**<codes>dfs</code> makes an depth first search over a matrix.<p>
	 * <p> The results are saved in the class variable {@link #adjacency adjacency}.
	 * Needs the correctly initalised field {@link #load load}.<p>
	 * @param isVariable Array with a boolean value for each reaction, if true, reaction will be not ignored. Computed by {@link #minMax(double[] firstVector, LinearProgramSolver solver,SparseMatrix rctMetArr) minMax}. Length is equal with the number of reactions in the SBML document. 
	 * @param dfsArray SparseMatrix with the graph of the minimal modules. Computed by {@link #computeMinModules(boolean[] isVariable, SparseMatrix rctMetArr, LinearProgramSolver solver) computeMinModules}.  Length and depth is equal with the number of reactions in the SBML document. 
	 */
	public void dfs(boolean[] isVariable, SparseMatrix dfsArray){
		boolean[] grey = new boolean[load.getNumR()];
		for(int i=0;i<numV;i++){
			if(isVariable[i]){
				if(!grey[i]){					
					ArrayList<Integer> neighbour = new ArrayList<Integer>();
					dfsHelp(i,isVariable,grey,dfsArray,neighbour);
					if(neighbour.size()>1){
						adjacency.add(neighbour);
					}
				}
			}
		}
	}
	
	
	/**<code>dfsHelp()</code> exists only for the recursion of <codes>dfs()</code>.<p>
	 * Needs correctly initalised field {@link #load load}.<p>
	 *  
	 * @param i Which knot is considered at the moment. i is the position in the array with all listed reactions, it is equal with the id's of the reactions in the SBML document. 
	 * @param isVariable Array with a boolean value for each reaction, if true, reaction will be not ignored. Computed by {@link #minMax(double[] firstVector, LinearProgramSolver solver,SparseMatrix rctMetArr) minMax}.  Length is equal with the number of reactions in the SBML document. 
	 * @param grey Marks all visited knots.  Length is equal with the number of reactions in the SBML document. 
	 * @param dfsArray SparseMatrix with the graph of the minimal modules. Computed by {@link #computeMinModules(boolean[] isVariable, SparseMatrix rctMetArr, LinearProgramSolver solver) computeMinModules}
	 * @param neighbour The inner list, that saves the id's of the connected components (modules). In {@link #dfs(boolean[] isVariable, SparseMatrix dfsArray) dfs} they are added to {@link #adjacency adjacency}.
	 */
	private void dfsHelp(int i, boolean[] isVariable, boolean[] grey, SparseMatrix dfsArray, ArrayList<Integer> neighbour){
		grey[i] = true;
		for(int k=0;k<load.getNumR();k++){
			if(isVariable[k]){
				if(!grey[k] && dfsArray.get(i, k) ==1){
					dfsHelp(k, isVariable, grey, dfsArray, neighbour);
				}
			}
		}
		neighbour.add(i);
	}
	
	
	/**<code>computeMinModules</code> calculates with linear programming the dependent 
	 * sets of the stoichiometric vectors of the variable reactions.<p>
	 * Needs correctly initalised field {@link #load load}.<p>
	 * 
	 * @param isVariable Array with a boolean value for each reaction, if true, reaction will be not ignored. Computed by {@link #minMax(double[] firstVector, LinearProgramSolver solver,SparseMatrix rctMetArr) minMax}.  Length is equal with the number of reactions in the SBML document. 
	 * @param rctMetArr SparseMatrix initialized by {@link #matrixBuild() matrixBuild}. Need to has the length of [number of reactions # number of metabolites].
	 * @param solver SCPSolver
	 * @return A SparseMatrix which represents a graph. All reactions are mapped on all reactions. A "1" represents a edge between two reactions, the rest is filled with "0".
	 */
	public SparseMatrix computeMinModules(boolean[] isVariable, SparseMatrix rctMetArr, LinearProgramSolver solver) {

		for (boolean x : isVariable) {
			if (x) {
				numV++;
			}
		}

		boolean[] notIgnore = new boolean[numV];
		if(numV>0){
			notIgnore[0] = true;
		}
		// every round one new lp for all variable reaction
		int numOfReactionOnLHS = 1;
		int reactionNum = 0;
		
		// start at e=1 (countVar=1 and notIgnore[0]=true) because first
		// equation can be skipped
		
		//array for depth first search 
		SparseMatrix dfsArray = new SparseMatrix(load.getNumR(),load.getNumR());
				
		//print
		float n = (float) numV;
		float print = 0;
		float add =  n/100;
		float count = 0;
		if(ignore){
			System.out.println("\t|0%       |10%      |20%      |30%      |40%      |50%      |60%      |70%      |80%      |90%      |100%");
			pr("\t");
		}
		
		for (int e = 1; e < load.getNumR(); e++) {
			if (isVariable[e]) {
				
				//print
				if(ignore){
					if(count == Math.round(print)){
						pr("_");
						print += add;
					}
				}
				else{
					System.out.println("Reaction number: "+count+" of "+numV);
				}
				count++;
				
				// setting all constraints
				ArrayList<LinearConstraint> constraints = new ArrayList<LinearConstraint>();
				int[] rememberRctPos = new int[reactionNum];

				// set ergebnis
				double[] ergebnis = new double[load.getNumS()];
				for (int u = 0; u < load.getNumS(); u++) {
					ergebnis[u] = -rctMetArr.get(e, u);
				}
				
				for (int i = 0; i < load.getNumS(); i++) {

					double[] con = new double[numOfReactionOnLHS];

					// set con
					int skipper = 0;
					for (int k = 0; k < reactionNum; k++) {
						if (notIgnore[k]) {
							con[skipper] = rctMetArr.get(k, i);
							rememberRctPos[skipper] = k;
							skipper++;
						}
					}
					LinearConstraint lc = new LinearEqualsConstraint(con,
							ergebnis[i], null);
					constraints.add(lc);
				}

				double[] objective = new double[numOfReactionOnLHS];
				LinearProgram lp = new LinearProgram(objective);
				lp.addConstraints(constraints);
				double[] result = solver.solve(lp);

				if (lp.isFeasable(result)) {
					
					
					
					ArrayList<Integer> neighbour = new ArrayList<Integer>();
					neighbour.add(e);
					for (int h = 0; h < numOfReactionOnLHS; h++) {
						if (Math.abs(result[h]) > 1e-4) {
							dfsArray.set(e,rememberRctPos[h],1);
							dfsArray.set(rememberRctPos[h],e,1);
						}
						
					}
				} else {
					notIgnore[reactionNum] = true;
					numOfReactionOnLHS++;
				}
				reactionNum++;
			}
		}
		pr("\n\n");
		return dfsArray;
	}

	
	/**<code>minMax()</code> execute for all reactions a flux variability analysis. 
	 * <p>
	 * The method compares every minimisation value and maximisation value of a flux balance analysis of a reaction. If both values are not equal,
	 * the reaction is variable and can be skipped in further computing.<p>
	 * Needs the correctly initalised field {@link #load load}.<p>
	 * <p>
	 * @param firstVector result vector of {@link #optimize(LinearProgramSolver, SparseMatrix) optimize}. Length is equal with the number of reactions in the SBML document. 
	 * @param solver SCPSolver
	 * @param rctMetArr SparseMatrix initialized by {@link #matrixBuild() matrixBuild}. Need to has the length of [number of reactions # number of metabolites].
	 * @return isVariable A boolean array which elements are true for all variable values. It has the length of the number of all reactions.
	 */
	public boolean[] minMax(double[] firstVector, LinearProgramSolver solver,
			SparseMatrix rctMetArr) {
		// set constraints
		ArrayList<LinearConstraint> constraints = new ArrayList<LinearConstraint>();
		for (int i = 0; i < load.getNumS(); i++) {
			double[] constraint = new double[load.getNumR()];
			for (int k = 0; k < load.getNumR(); k++) {
				constraint[k] = rctMetArr.get(k, i);
			}
			Species s = load.getModel().getSpecies(i);
			if (s.getBoundaryCondition()) {
				continue;
			} else {
				LinearConstraint lincon = (new LinearEqualsConstraint(
						constraint, 0, null));
				constraints.add(lincon);

			}
		}

		// biomass as a constraint
		LinearConstraint biomass = new LinearEqualsConstraint(load.getObjectiveFunction(), biomassOptValue-1e-5, "OptBiomass");
		constraints.add(biomass);

		// set objective function and detect differences for each reaction, to
		// find V
		boolean[] isVariable = new boolean[load.getNumR()];
		double[] maxValue = new double[load.getNumR()];
		double[] minValue = new double[load.getNumR()];
		for (int q = 0; q < load.getNumR(); q++) {
			maxValue[q] = firstVector[q];
			minValue[q] = firstVector[q];
		}
		
		//print:
		float print = 0;
		float a = 0;
		float add = 0;
		if(ignore){
			System.out.println("\t|0%       |10%      |20%      |30%      |40%      |50%      |60%      |70%      |80%      |90%      |100%");
			a = (float) load.getNumR();
			pr("\t");
			add = (a/100);
		}
		
		for (int r = 0; r < load.getNumR(); r++) {
			
			if(ignore){
				if(r == Math.ceil(print)){
					pr("_");
					print += add;
				}
			}
			else{
				System.out.println("Reaction number: "+r+" of " + load.getNumR());
			}
			
			if (!isVariable[r]) { // if you do not know if its variable, we test
		
				double[] newObjective = new double[load.getNumR()];
				newObjective[r] = 1;

				LinearProgram lpMax = new LinearProgram(newObjective);
				lpMax.addConstraints(constraints);
				lpMax.setLowerbound(load.getLowerBound());
				lpMax.setUpperbound(load.getUpperBound());
				lpMax.setMinProblem(false);

				LinearProgram lpMin = new LinearProgram(newObjective);
				lpMin.addConstraints(constraints);
				lpMin.setLowerbound(load.getLowerBound());
				lpMin.setUpperbound(load.getUpperBound());
				lpMin.setMinProblem(true);

				double[] solvedMax = solver.solve(lpMax);
				double[] solvedMin = solver.solve(lpMax);
				if (!lpMin.isFeasable(solvedMin)) {
					solvedMin = firstVector;
				}
				if (!lpMax.isFeasable(solvedMax)) {
					solvedMax = firstVector;
				}

				for (int e = 0; e < load.getNumR(); e++) {
					if (Math.abs(firstVector[e] - solvedMax[e]) > 1e-4
							|| Math.abs(firstVector[e] - solvedMin[e]) > 1e-4) {
						isVariable[e] = true;
						maxValue[e] = Math.max(
								Math.max(maxValue[e], solvedMax[e]),
								solvedMin[e]);
						minValue[e] = Math.min(
								Math.min(minValue[e], solvedMin[e]),
								solvedMax[e]);
					}
				}
			}
		}
		pr("\n");
		return isVariable;
	}

	
	/**<code>optimize()</code> execute a flux balance analysis for the biomass reaction in the model.
	 * <p>
	 * The biomass is being maximized.
	 * Needs the correctly initalised field {@link #load load}.<p>
	 * <p> 
	 * @param solver SCPSolver
	 * @param rctMetArr SparseMatrix initialized by {@link #matrixBuild() matrixBuild}. Need to has the length of [number of reactions # number of metabolites].
	 * @return A double array with the result of the linear programm, which has optimized the biomass. 
	 * @throws LpSolveException
	 */
	public double[] optimize(LinearProgramSolver solver, SparseMatrix rctMetArr)
			throws LpSolveException {

		LinearProgram lp = new LinearProgram(load.getObjectiveFunction());

		for (int i = 0; i < load.getNumS(); i++) {
			double[] constraint = new double[load.getNumR()];
			for (int k = 0; k < load.getNumR(); k++) {
				constraint[k] = rctMetArr.get(k, i);
			}
			Species s = load.getModel().getSpecies(i);
			if (s.getBoundaryCondition()) {
				continue;
			} else {
				lp.addConstraint(new LinearEqualsConstraint(constraint, 0, null));
			}
		}

		lp.setLowerbound(load.getLowerBound());
		lp.setUpperbound(load.getUpperBound());

		double[] sol = solver.solve(lp);

		biomassOptValue = sol[load.getBiomassOptValuePos()];
		return sol;
	}

	
	/**<code>matrixBuild()</code> builds a matrix with the stoichiometric values of all reactions.<p>
	 * <p>
	 * Needs the correctly initalised field {@link #load load}.<p>
	 * @return A SparseMatrix which stores the stochiometric values for all metabolite of each reaction. 
	 */
	public SparseMatrix matrixBuild() {
		SparseMatrix rctMetArr = new SparseMatrix(load.getNumR(),
				load.getNumS());

		// calculates the values of rctMetArr

		for (int i = 0; i < load.getNumR(); i++) { // all reactions
			Reaction r = load.getModel().getReaction(i);
			int w = r.getNumReactants();
			int z = r.getNumProducts();

			for (int k = 0; k < load.getNumS(); k++) { // all metabolites
				// compares reactants and products with metabolites (species)
				// and set the stochiometric value in the matrix

				for (int t = 0; t < w; t++) { // all reactants
					SpeciesReference s = r.getReactant(t);
					Species sp = s.getSpeciesInstance();
					String id = sp.getId();
					if (id == load.getMet()[k]) {
						rctMetArr.set(i, k,
								-(r.getReactant(t).getStoichiometry()));
					}
				}

				for (int t = 0; t < z; t++) { // all products
					SpeciesReference s = r.getProduct(t);
					Species sp = s.getSpeciesInstance();
					String id = sp.getId();
					if (id == load.getMet()[k]) {
						rctMetArr.set(i, k, r.getProduct(t).getStoichiometry());
					}
				}
			}
		}

		return rctMetArr;
	}
	/**Shortcut for "System.out.print"
	 * 
	 * @param s
	 */
	public void pr(String s) {
		System.out.print(s);
	}
}
