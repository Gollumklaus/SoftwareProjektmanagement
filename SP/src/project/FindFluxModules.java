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

// boolean = s.getBoundaryCondition();

public class FindFluxModules{


	SBMLLoad load;
	double biomassOptValue = 0;
	ArrayList<ArrayList<Integer>> adjacency = new ArrayList<ArrayList<Integer>>();
	ArrayList<ArrayList<Reaction>> vVectors = new ArrayList<ArrayList<Reaction>>();
	ArrayList<Integer> whichRctIsModule = new ArrayList<Integer>();
	int numV = 0;
	public boolean ignore;

	public double getBiomassOptValue() {
		return biomassOptValue;
	}

	public ArrayList<Integer> getWhichRctIsModule() {
		return whichRctIsModule;
	}

	public ArrayList<ArrayList<Integer>> getAdjacency() {
		return adjacency;
	}

	public FindFluxModules(SBMLLoad load, boolean ignore) {
		this.load = load;
		this.ignore = ignore;
	}

	// findFlux: "main replacement"... call all needed methods
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

	// print modules
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
	public void dfsHelp(int i, boolean[] isVariable, boolean[] grey, SparseMatrix dfsArray, ArrayList<Integer> neighbour){
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
	
	// compute minimal Modules
	public SparseMatrix computeMinModules(boolean[] isVariable, SparseMatrix rctMetArr, LinearProgramSolver solver) {

		for (boolean x : isVariable) {
			if (x) {
				numV++;
			}
		}

		boolean[] notIgnore = new boolean[numV];
		notIgnore[0] = true;

		// every round one new lp for all variable reaction
		int numOfReactionOnLHS = 1;
		int reactionNum = 1;
		// start at e=1 (countVar=1 and notIgnore[0]=true) because first
		// equation can be skipped
		
		//array for deepsearch
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

	
	
	
	
	
	
	
	// calculate min and max values for each reaction
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
					isVariable[r] = true;
				}
				if (!lpMax.isFeasable(solvedMax)) {
					solvedMax = firstVector;
					isVariable[r] = true;
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

	// calculate optimum S*v=0. find all vectors which solve the equation and
	// maximize x.
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

	// build matrix with reaction and metabolites
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

	public void pr(String s) {
		System.out.print(s);
	}
}
