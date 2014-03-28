package project;


import java.io.FileNotFoundException;
import java.util.ArrayList;

import lpsolve.LpSolveException;

import org.sbml.jsbml.*;

import scpsolver.constraints.LinearBiggerThanEqualsConstraint;
import scpsolver.constraints.LinearConstraint;
import scpsolver.constraints.LinearEqualsConstraint;
import scpsolver.constraints.LinearSmallerThanEqualsConstraint;
import scpsolver.lpsolver.LinearProgramSolver;
import scpsolver.lpsolver.SolverFactory;
import scpsolver.problems.LinearProgram;
import scpsolver.util.SparseMatrix;


// boolean = s.getBoundaryCondition();

public class FindFluxModules{
	
	
	
	
	//findFlux: "main replacement"... call all needed methods
	public static void findFlux(SBMLLoad load) throws LpSolveException, FileNotFoundException{
		
		//get values from SBMLLoad
		int numR = load.getNumR();
		int numS = load.getNumS();
		Model m = load.getModel();
		String[] met = load.getMet();
		double[] objectiveFunction = load.getObjectiveFunction();
		int biomassOptValuePos = load.getBiomassOptValuePos();
		double[] low = load.getLowerBound();
		double[] upp = load.getUpperBound();
		
		//build matrix
		SparseMatrix rctMetArr = matrixBuild(numR,numS,met,m);
		
		
		
		//find optimum using LP
		System.out.println("Loading LP Solver");
		LinearProgramSolver solver  = SolverFactory.newDefault(); 
		double biomassOptValue = 0;
		double[] vBiomass = optimize(biomassOptValuePos,biomassOptValue,solver, rctMetArr,numR,numS,m,objectiveFunction,low,upp);
		
		//find min max values
		minMax(biomassOptValue,vBiomass,solver,numR,numS,m,objectiveFunction,low,upp, rctMetArr);
	}	

	//calcualte min and max values for each reaction
	public static void minMax(double biomassOptValue, double[] firstVector, LinearProgramSolver solver, int numR, int numS, Model m, double[] objectiveFunction, double[] low, double[] upp, SparseMatrix rctMetArr){
		
		//set constraints
		ArrayList<LinearConstraint> con = new ArrayList<LinearConstraint>();
			for(int i=0;i<numS;i++){		
				double[] constraint = new double[numR];
				for(int k=0;k<numR;k++){	
					constraint[k] = rctMetArr.get(k,i) ;				
				}
				Species s = m.getSpecies(i);
				if(s.getBoundaryCondition()){
					continue;
				}
				else{
					LinearConstraint lincon = (new LinearEqualsConstraint(constraint, 0,null));
					con.add(lincon);
					
				}
			}
			
			for(int k=0;k<numR;k++){	
				double[] constraint = new double[numR]; 			
				constraint[k]=1;
				if(low[k]<-1003.0 && upp[k]>1003.0){
					low[k]=-1000.0;
					upp[k]=1000;
				}
				else{
					upp[k]=1000;
				}	
				LinearConstraint bigcon = (new LinearBiggerThanEqualsConstraint(constraint, low[k],null));
				LinearConstraint smalcon = (new LinearSmallerThanEqualsConstraint(constraint, upp[k],null));
				con.add(bigcon);
				con.add(smalcon);
			}
			
		//biomass as a constraint
			LinearConstraint biomass = (new LinearEqualsConstraint(objectiveFunction, biomassOptValue,"OptBiomass"));
			con.add(biomass);
			
			
			
			
		//set objective function and detect differences for each reaction, to find V
		double[] isNotConstant = objectiveFunction;
		
		for(int r=0;r<numR;r++){
			if(isNotConstant[r] != 1.0){
			double[] newObjective = new double[numR];
			newObjective[r] = 1;
			
			LinearProgram lp = new LinearProgram(newObjective);
			lp.addConstraints(con);
			lp.setMinProblem(false); 
			double[] solved = solver.solve(lp);

			for(int w=0;w<numR;w++){
				System.out.print(solved[w]+" ");
			}
			System.out.println();
			
int count=0;
for(int zu=0;zu<numR;zu++){
	if(solved[zu] != 0){
		count++;
	}
	
	for(int e=0;e<numR;e++){
		if(firstVector[e] !=  solved[e] ){
			isNotConstant[e] = 1.0;
		}
	}
	
}
System.out.println("#!=0: " + count + " /Number: " + r);	
			} 
		}		
	}
	
	//calculate optimum S*v=0. find all vectors which solve the equation and maximize x.
	public static double[] optimize(int biomassOptValuePos, double biomassOptValue, LinearProgramSolver solver, SparseMatrix rctMetArr, int numR , int numS, Model m, double[] objectiveFunction, double[] low, double[]upp) throws LpSolveException {
		
		LinearProgram lp = new LinearProgram(objectiveFunction); 
	      System.out.println("blau");	      
		for(int i=0;i<numS;i++){		
			double[] constraint = new double[numR];
			for(int k=0;k<numR;k++){	
				constraint[k] = rctMetArr.get(k,i) ;				
			}
			Species s = m.getSpecies(i);
			if(s.getBoundaryCondition()){
				continue;
			}
			else{
				lp.addConstraint(new LinearEqualsConstraint(constraint, 0,null)); 
			}
		}

		for(int k=0;k<numR;k++){	
			double[] constraint = new double[numR]; 			
			constraint[k]=1;
			Reaction r = m.getReaction(k);
	
			if(low[k]<-1003.0 && upp[k]>1003.0){
				low[k]=-1000.0;
				upp[k]=1000;
			}
			else{
				upp[k]=1000;
			}			
			lp.addConstraint(new LinearBiggerThanEqualsConstraint(constraint, low[k],null)); 
			lp.addConstraint(new LinearSmallerThanEqualsConstraint(constraint, upp[k],null));
		}

		lp.setMinProblem(false); 
		
		double[] sol = solver.solve(lp);
		
		for(int i=0;i<numR;i++){
			System.out.print(sol[i]+" ");
		}
		biomassOptValue = sol[biomassOptValuePos];
	
System.out.println("\nbiomassOptValue: " + biomassOptValue);
int count=0;
for(int zu=0;zu<numR;zu++){
	if(sol[zu] != 0){
		count++;
	}
}
System.out.println("#!=0: " + count);

		
		return sol;
	}

	//build matrix with reaction and metabolites	
	public static SparseMatrix matrixBuild(int numR,int numS, String[] met, Model m){
		SparseMatrix rctMetArr= new SparseMatrix(numR,numS);
		//array with reactions and metabolites 
		//double[][] rctMetArr = new double[numR][numS];
		int count =0;
		int correct = 0;
		

		
		
		System.out.println("Creating Reaction/Metabolites Matrix...");
		//calculates the values of rctMetArr
		for(int i=0;i<numR;i++){					//all reactions
			for(int k=0;k<numS;k++){				//all metabolites
			//compares reactants and products with metabolites (species) and set the stochiometric value in the matrix
				Reaction r = m.getReaction(i);
				int w = r.getNumReactants();		
					for(int t=0;t<w;t++){			//all reactants
						SpeciesReference s = r.getReactant(t);
						Species sp = s.getSpeciesInstance();
						String id = sp.getId();
						if(id == met[k]){
							rctMetArr.set(i,k,-(r.getReactant(t).getStoichiometry())); correct++;
						}
					}
				int z = r.getNumProducts();
					for(int t=0;t<z;t++){			//all products
						SpeciesReference s = r.getProduct(t);
						Species sp = s.getSpeciesInstance();
						String id = sp.getId();
						if(id == met[k]){
							rctMetArr.set(i,k,r.getProduct(t).getStoichiometry()); correct++;
						}
					}
			}
		}
		
/*
		//prints matrix
		for(int i=0;i<numR;i++){
			for(int k=0;k<numS;k++){
					System.out.print(rctMetArr[i][k]);
			}
			System.out.println();
		}
*/	
		
		//tests
		for(int i=0;i<numR;i++){					//all reactions
			for(int k=0;k<numS;k++){
				if(rctMetArr.get(i,k)!=0){
					count++;
				}
			}
		}
		double wus=0;
		int error=0;
		
		for(int i=0;i<numR;i++){
			Reaction r = m.getReaction(i);
			try{
				wus += r.getNumReactants();
			}
			catch(NullPointerException tt){
				error++;
			}
			
			try{
				wus += r.getNumProducts();
			}
			catch(NullPointerException tt){
				error++;
			}
		}
		
//System.out.println(count+":"+correct+":"+wus+":"+error);

		return rctMetArr;
	}

	
	
	
}


