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
import scpsolver.constraints.Constraint;
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
ArrayList<ArrayList<Integer>> adjacency = new ArrayList<ArrayList<Integer>>();
ArrayList<ArrayList<Reaction>> vVectors = new ArrayList<ArrayList<Reaction>>();
ArrayList<Integer> whichRctIsModule = new ArrayList<Integer>();
int numV = 0;


	public double getBiomassOptValue(){
		return biomassOptValue;
	}
	public ArrayList<Integer> getWhichRctIsModule(){
		return whichRctIsModule;
	}
	public ArrayList<ArrayList<Integer>> getAdjacency(){
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
		boolean[] isVariable = minMax(vBiomass,solver,rctMetArr);
		
		
		System.out.println("Finding minimal modules");
		computeMinModules(isVariable, rctMetArr,solver);
		
		System.out.println("printing modules");
		printModules();
	}	
	
	
	//print modules
	public void printModules(){
		for(int i=0;i<adjacency.size();i++){
			ArrayList<Integer> t = adjacency.get(i);
			int c = t.size();
			pr("\nModul number "+ i+":\t");
			for(int w=0;w<c;w++){
				System.out.print(t.get(w)+" ");
			}
			//pr("\n");
		}
	}
	
	
	//compute minimal Modules
	public void computeMinModules(boolean[] isVariable, SparseMatrix test, LinearProgramSolver solver){	
	
		for(boolean x: isVariable){
			if(x){numV++;}
		} 
							
		boolean[] notIgnore = new boolean[numV];
		notIgnore[0] = true; 
		//every round one new lp for all variable reaction
		int numOfReactionOnLHS = 1;
		int reactionNum = 1;
					//start at e=1 (countVar=1 and notIgnore[0]=true) because first equation can be skipped
		for(int e=1;e<load.getNumR();e++){
			if(isVariable[e]){
				//print position
				//pr("LP: "+ reactionNum + "\n");
				
				//setting all constraints
				ArrayList<LinearConstraint> constraints = new ArrayList<LinearConstraint>();
				int[] rememberRctPos = new int[reactionNum];
				
				//set ergebnis
				double[] ergebnis =  new double[load.getNumS()];
				for(int u=0;u<load.getNumS();u++){
					ergebnis[u] = -test.get(e,u);
				} 
				
				for(int i=0;i<load.getNumS();i++){
					
					double[] con = new double[numOfReactionOnLHS];
					
					//set con
					int skipper = 0;
					for(int k=0;k<reactionNum;k++){
						if(notIgnore[k]){
							con[skipper] = test.get(k, i);
							rememberRctPos[skipper] = k;
							skipper++;
						}
					}
					LinearConstraint lc = new LinearEqualsConstraint(con, ergebnis[i], null);
					constraints.add(lc);
				}
			
				double[] objective = new double[numOfReactionOnLHS];
				LinearProgram lp = new LinearProgram(objective);
				lp.addConstraints(constraints);
				double[] result = solver.solve(lp);
				
					/*									//PRINT
														ArrayList<Constraint> litt = lp.getConstraints();
														for(Constraint p : litt){
															double[] weigth = ((LinearConstraint) p).getC();
															for(double kl : weigth){System.out.print(kl + " ");};															
															System.out.println( " = " + p.getRHS());														
														}
														pr("\nResult: ");
														for(double j:result){System.out.print(j+"    ");}pr("\n");
				*/
														
				if(lp.isFeasable(result)){
					ArrayList<Integer> neighbour = new ArrayList<Integer>();
					neighbour.add(e);
					for(int h=0;h<numOfReactionOnLHS;h++){
						if(Math.abs(result[h])>1e-4){
							neighbour.add(rememberRctPos[h]);
						}
					}
					adjacency.add(neighbour);
					whichRctIsModule.add(e);
				}
				else{
					notIgnore[reactionNum] = true;
					numOfReactionOnLHS++;					 
				}
				reactionNum++;	
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
		boolean[] isVariable = new boolean[load.getNumR()];
		double[] maxValue = new double[load.getNumR()];
		double[] minValue = new double[load.getNumR()];
		for(int q=0;q<load.getNumR();q++){
			maxValue[q] = firstVector[q];
			minValue[q] = firstVector[q];
		}
		
		
		for(int r=0;r<load.getNumR();r++){
			if(!isVariable[r]){	//if u dont know if its variable, we test
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
														System.out.println("solved Max  " +solvedMax.length);
														System.out.println("solved Min  " +solvedMax.length);
														System.out.println("numR  "+load.getNumR());
														System.out.println("numS  "+load.getNumS());
				if(!lpMin.isFeasable(solvedMin)){
					solvedMin = firstVector;
				}
				if(!lpMax.isFeasable(solvedMax)){
					solvedMax = firstVector;
				}
				
				for(int e=0;e<load.getNumR();e++){
					if(Math.abs(firstVector[e]-solvedMax[e])>1e-4 || Math.abs(firstVector[e] - solvedMin[e])>1e-4){
						isVariable[e] = true;
						maxValue[e] = Math.max(Math.max(maxValue[e], solvedMax[e]), solvedMin[e]);
						minValue[e] = Math.min(Math.min(minValue[e], solvedMin[e]), solvedMax[e]);
					}
				}
										//OUTPUTPRINT: print solvedMax/solvedMin & #numbers in both arrays unequal zero
											pr("\tsolvedMax: ");for(int pp=0;pp<load.getNumR();pp++){System.out.print(solvedMax[pp]+" ");	}pr("\n");
											pr("\tsolvedMin: ");for(int pp=0;pp<load.getNumR();pp++){System.out.print(solvedMin[pp]+" ");	}pr("\n");
											int c=0;for(int i=0;i<load.getNumR();i++){if(Math.abs(solvedMax[i])<1e-4){c++;}}System.out.println("\tMax gleich null: "+c);
											int t=0;for(int i=0;i<load.getNumR();i++){if(Math.abs(solvedMin[i])<1e-4){t++;}}System.out.println("\tMin gleich null: "+t);
											pr("\tisConstant: ");for(int i=0;i<load.getNumR();i++){String w=null;if(!isVariable[i]){w="C";}else{w=".";}System.out.print(w);}
											pr("\n");} 
		}	
		for(int uh=0;uh<load.getNumR();uh++){
			if(isVariable[uh]){
			System.out.println(maxValue[uh]+" : "+minValue[uh]);
		
			}
		}	
		pr("{");
		for(int z=0;z<load.getNumR();z++){
			System.out.print(isVariable[z] + ",");
		}
		pr("}");
		return isVariable;
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


