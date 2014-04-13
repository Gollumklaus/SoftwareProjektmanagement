package project;


import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.Set;

import lpsolve.LpSolveException;

import org.sbml.jsbml.*;

import scpsolver.constraints.LinearConstraint;
import scpsolver.constraints.LinearEqualsConstraint;
import scpsolver.lpsolver.LinearProgramSolver;
import scpsolver.lpsolver.SolverFactory;
import scpsolver.problems.LinearProgram;
import scpsolver.util.SparseMatrix;

// boolean = s.getBoundaryCondition();

public class FindFluxModules implements Modules{
	
	private class ModuleImpl implements Module {

		@Override
		public Set<Reaction> getReactions() {
			// TODO Auto-generated method stub
			return null;
		}

		@Override
		public Map<Species, Double> getInterface() {
			// TODO Auto-generated method stub
			return null;
		}
		
	}
	
SBMLLoad load;
double biomassOptValue = 0;
ArrayList<ArrayList<Reaction>> adjacency = new ArrayList<ArrayList<Reaction>>();
ArrayList<ArrayList<Reaction>> vVectors = new ArrayList<ArrayList<Reaction>>();
ArrayList<Integer> whichRctIsModule = new ArrayList<Integer>();
int numV = 0;


	public double getBiomassOptValue(){
		return biomassOptValue;
	}
	public ArrayList<Integer> getWhichRctIsModule(){
		return whichRctIsModule;
	}
	public ArrayList<ArrayList<Reaction>> getAdjacency(){
		return adjacency;
	}

	public FindFluxModules(SBMLLoad load){
		this.load = load;
	}
	
	
	//findFlux: "main replacement"... call all needed methods
	public void findFlux() throws LpSolveException, FileNotFoundException{
		
		//build matrix
		System.out.println("Creating Reaction/Metabolites Matrix");
		SparseMatrix rctMetArr = matrixBuild();
				
		//find optimum using LP
		System.out.println("Loading LP Solver");
		LinearProgramSolver solver  = SolverFactory.newDefault(); 
		
		System.out.println("Solving LP for biomassreaction");
		double[] vBiomass = optimize(solver, rctMetArr);
		
		//find min max values
		System.out.println("Finding set V");
		boolean[] isConstant = minMax(vBiomass,solver,rctMetArr);
		
		
		System.out.println("Finding minimal modules");
		computeMinModules(isConstant, rctMetArr,solver);
		
		System.out.println("printing modules");
		printModules();
	}	
	
	
	//print modules
	public void printModules(){
		for(int i=0;i<numV;i++){
			ArrayList<Reaction> t = adjacency.get(i);
			int c = t.size();
			for(int w=0;w<c;w++){
				Reaction r = t.get(w);
				String p = r.getId();
				System.out.println(p);
			}
			//pr("\n");
		}
	}
	//compute minimal Modules
	public void computeMinModules(boolean[] isConstant, SparseMatrix rctMetArr, LinearProgramSolver solver){	
	
		for(boolean x: isConstant){
			if(!x){numV++;}
		}
									//OUTPUTPRINT: print numV & isConstant (C = constant, . = variable)
										pr("final isConstant:  ");for(int i=0;i<load.getNumR();i++){String w=null;if(isConstant[i]){w="C";}else{w=".";}System.out.print(w);}System.out.println("\nrct not const"+numV+"\t alle rct"+load.getNumR());
		
		SparseMatrix x = new SparseMatrix(numV, load.getNumS());
		double[] objectiv = new double[numV];
		LinearProgram lp = new LinearProgram(objectiv);
		lp.setMinProblem(false);
		
		ArrayList<LinearConstraint> constraints = new ArrayList<LinearConstraint>();
		int posInX = 0;
		int posInModules = 0;
		
		for(int i=0;i<load.getNumR();i++){
			if(!isConstant[i]){
				//set constraints
				double[] sR = new double[load.getNumS()];
				for(int k=0;k<load.getNumS();k++){
					sR[k] = - rctMetArr.get(i,k);
				}
				for(int t=0;t<load.getNumS();t++){
					double[] columnOfX = new double[load.getNumS()];
					for(int s=0;s<numV;s++){
						columnOfX[t] = x.get(s,t);
					}
					
					LinearConstraint lc = new LinearEqualsConstraint(columnOfX, sR[t], null);
					constraints.add(lc);
					 
				}
			
				double[] solved = solver.solve(lp);
				
				
				if(lp.isFeasable(solved)){
					ArrayList<Reaction> neighbour = new ArrayList<Reaction>();
					for(int h=0;h<numV;h++){
						if(Math.abs(solved[h])>1e-4){
							Reaction r = load.getModel().getReaction(h);
							neighbour.add(r);
						}
					}
					adjacency.add(posInModules,neighbour);
					posInModules++;
					whichRctIsModule.add(i);
				}
				else{
					for(int u=0;u<load.getNumS();u++){
						x.set(posInX,u,sR[u]);
					}
					posInX++;
				}
			}
		}	
	}
	
	

	//calculate min and max values for each reaction
	public boolean[] minMax(double[] firstVector, LinearProgramSolver solver, SparseMatrix rctMetArr){
		
		//set constraints
		ArrayList<LinearConstraint> constraints = new ArrayList<LinearConstraint>();
			for(int i=0;i<load.getNumS();i++){	
				double[] constraint = new double[load.getNumR()];
				for(int k=0;k<load.getNumR();k++){	
					constraint[k] = rctMetArr.get(k,i) ;				
				}
				Species s = load.getModel().getSpecies(i);
				if(s.getBoundaryCondition()){
					continue;
				}
				else{
					LinearConstraint lincon = (new LinearEqualsConstraint(constraint, 0,null));
					constraints.add(lincon);
					
				}
			}
			

			
		//biomass as a constraint
			LinearConstraint biomass = new LinearEqualsConstraint(load.getObjectiveFunction(), biomassOptValue,"OptBiomass");
			constraints.add(biomass);
			
		//set objective function and detect differences for each reaction, to find V
		boolean[] isConstant = new boolean[load.getNumR()];
		isConstant[load.getBiomassOptValuePos()] = true;
	
	
	
		for(int r=0;r<load.getNumR();r++){
			if(!isConstant[r]){
								//OUTPUTPRINT: print which rct 
									pr("reaction Number: ");System.out.println(r);
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
	
			
				for(int e=0;e<load.getNumR();e++){
					if(Math.abs(firstVector[e]-solvedMax[e])>1e-4 || Math.abs(firstVector[e] - solvedMin[e])>1e-4){
						isConstant[e] = true;
					}
				}
										//OUTPUTPRINT: print solvedMax/solvedMin & #numbers in both arrays unequal zero
											pr("\tsolvedMax: ");for(int pp=0;pp<load.getNumR();pp++){System.out.print(solvedMax[pp]+" ");	}pr("\n");
											pr("\tsolvedMin: ");for(int pp=0;pp<load.getNumR();pp++){System.out.print(solvedMax[pp]+" ");	}pr("\n");
											int c=0;for(int i=0;i<load.getNumR();i++){if(Math.abs(solvedMax[i])<1e-4){c++;}}System.out.println("\tMax ungleich null: "+c);
											int t=0;for(int i=0;i<load.getNumR();i++){if(Math.abs(solvedMin[i])<1e-4){t++;}}System.out.println("\tMin ungleich null: "+t);
											pr("\tisConstant: ");for(int i=0;i<load.getNumR();i++){String w=null;if(isConstant[i]){w="C";}else{w=".";}System.out.print(w);}
											pr("\n");} 
		}	

		return isConstant;
	}
	
	//calculate optimum S*v=0. find all vectors which solve the equation and maximize x.
	public double[] optimize(LinearProgramSolver solver, SparseMatrix rctMetArr) throws LpSolveException {
		
		LinearProgram lp = new LinearProgram(load.getObjectiveFunction()); 
     
		for(int i=0;i<load.getNumS();i++){		
			double[] constraint = new double[load.getNumR()];
			for(int k=0;k<load.getNumR();k++){	
				constraint[k] = rctMetArr.get(k,i) ;				
			}
			Species s = load.getModel().getSpecies(i);
			if(s.getBoundaryCondition()){
				continue;
			}
			else{
				lp.addConstraint(new LinearEqualsConstraint(constraint, 0,null)); 
			}
		}

		lp.setLowerbound(load.getLowerBound());
		lp.setUpperbound(load.getUpperBound());
		
		

	
		double[] sol = solver.solve(lp);

		biomassOptValue = sol[load.getBiomassOptValuePos()];
		for(int pp=0;pp<load.getNumR();pp++){System.out.print(sol[pp]+" ");	}pr("\n");
		return sol;
	}
	
	//build matrix with reaction and metabolites	
	public SparseMatrix matrixBuild(){
		SparseMatrix rctMetArr= new SparseMatrix(load.getNumR(),load.getNumS());
	

		//calculates the values of rctMetArr
	
		for(int i=0;i<load.getNumR();i++){					//all reactions
			Reaction r = load.getModel().getReaction(i);
			int w = r.getNumReactants();	
			int z = r.getNumProducts();
			
			for(int k=0;k<load.getNumS();k++){				//all metabolites
			//compares reactants and products with metabolites (species) and set the stochiometric value in the matrix
					
					for(int t=0;t<w;t++){			//all reactants
						SpeciesReference s = r.getReactant(t);
						Species sp = s.getSpeciesInstance();
						String id = sp.getId();
						if(id == load.getMet()[k]){
							rctMetArr.set(i,k,-(r.getReactant(t).getStoichiometry())); 
						}
					}
				
					for(int t=0;t<z;t++){			//all products
						SpeciesReference s = r.getProduct(t);
						Species sp = s.getSpeciesInstance();
						String id = sp.getId();
						if(id == load.getMet()[k]){
							rctMetArr.set(i,k,r.getProduct(t).getStoichiometry()); 
						}
					}
			}
		}
		
		return rctMetArr;
	}

	public Modules Modules(){
		throw new RuntimeException("Noch nicht implementiert.");
	}
//	@Override
//	public Set<Reaction> getReactions() {
//		// TODO Auto-generated method stub
//		return null;
//	}
//	@Override
//	public Map<Species, Double> getInterface() {
//		// TODO Auto-generated method stub
//		return null;
//	}
	@Override
	public Collection<Module> getModules() {
		// TODO Auto-generated method stub
		
		ArrayList<Module> modules = new ArrayList<Module>();
		
		modules.add(new ModuleImpl());
		
		return modules;
	}
	public void pr(String s){
		System.out.print(s);
	}
}


