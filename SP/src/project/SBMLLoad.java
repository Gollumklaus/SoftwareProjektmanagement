package project;

import java.io.File;
import java.io.IOException;

import javax.xml.stream.XMLStreamException;

import org.sbml.jsbml.KineticLaw;
import org.sbml.jsbml.LocalParameter;
import org.sbml.jsbml.Model;
import org.sbml.jsbml.Reaction;
import org.sbml.jsbml.SBMLDocument;
import org.sbml.jsbml.SBMLReader;
import org.sbml.jsbml.Species;


/**<code>SBMLLoad</code> loads the SBML document and provide all loaded variable due to getter and setter methods to the public.
 * <p><p>
 * methods: <p>
 * <code>public void loadSBML(String path)</code><p>
 * @author guru
 */

public class SBMLLoad {

	// these variables get set by loadSBML
	/**Number of all reactions in the SBML document.
	 */
	private int numR;
	/**Number of all metabolites in the SBML document.
	 */
	private int numS;
	/**The model loaded from the SBML Document.
	 */
	private Model model;
	/**Id's of all reactions.
	 */
	private String[] rct;
	/**Id's of all metabolites.
	 */
	private String[] met;
	/**The SBML document as a file.
	 */
	private File file;
	/**The objectiveFunction shows, which reaction should be optimized.<p>
	 * the double[] consist of zeros for reaction with no optimization and ones for reactions with optimization<p>
	 *<p>
	 *Should have the length of numR!
	 */
	private double[] objectiveFunction;
	// upper bound and lower bound
	/**Consist of the lowerBounds of the SBMLDocument.<p>
	 * <p>
	 *Should have the length of numR!
	 */
	private double[] lowerBound;
	/**Consist of the upperBounds of the SBMLDocument.<p>
	 * <p>
	 *Should have the length of numR!
	 */
	private double[] upperBound;
	/**
	 * Saves the position of the biomass reaction in the array objectiveFunction.
	 */
	private int biomassOptValuePos;

	// getter
	/**numR is initialize by {@link #loadSBML(String path) loadSBML}.
	 * Before that, the variable is null.
	 * @return {@link #numR numR}
	 */
	public int getNumR() {
		return numR;
	}
	/**numS is initialize by {@link #loadSBML(String path) loadSBML}.
	 * Before that, the variable is null.
	 * @return {@link #numS numS}
	 */
	public int getNumS() {
		return numS;
	}
	/**model is initialize by {@link #loadSBML(String path) loadSBML}.
	 * Before that, the variable is null.
	 * @return {@link #model model}
	 */
	public Model getModel() {
		return model;
	}
	/**rct is initialize by {@link #loadSBML(String path) loadSBML}.
	 * Before that, the variable is null.
	 * @return {@link #rct rct}
	 */
	public String[] getRct() {
		return rct;
	}
	/**met is initialize by {@link #loadSBML(String path) loadSBML}.
	 * Before that, the variable is null.
	 * @return {@link #met met}
	 */
	public String[] getMet() {
		return met;
	}
	/**file is initialize by {@link #loadSBML(String path) loadSBML}.
	 * Before that, the variable is null.
	 * @return {@link #file file}
	 */
	public File getFile() {
		return file;
	}
	/**objectiveFunction is initialize by {@link #loadSBML(String path) loadSBML}.
	 * Before that, the variable is null.
	 * @return {@link #objectiveFunction objectiveFunction}
	 */
	public double[] getObjectiveFunction() {
		return objectiveFunction;
	}
	/**lowerBound is initialize by {@link #loadSBML(String path) loadSBML}.
	 * Before that, the variable is null.
	 * @return {@link #lowerBound lowerBound}
	 */
	public double[] getLowerBound() {
		return lowerBound;
	}
	/**upperBound is initialize by {@link #loadSBML(String path) loadSBML}.
	 * Before that, the variable is null.
	 * @return {@link #upperBound upperBound}
	 */
	public double[] getUpperBound() {
		return upperBound;
	}
	/**biomassOptValuePos is initialize by {@link #loadSBML(String path) loadSBML}.
	 * Before that, the variable is null.
	 * @return {@link #biomassOptValuePos biomassOptValuePos}
	 */
	public int getBiomassOptValuePos() {
		return biomassOptValuePos;
	}

	// load model
	/**<code>loadSBML()</code> load the Model, the number of Reactions/Species, the lower/upper bounds and the objective coefficients from the provided SBML document.
	 * <p>The method sets also the objective function.
	 * <p> 
	 * @param path the hole path of the SBML document
	 * @throws XMLStreamException
	 * @throws IOException
	 */
	public void loadSBML(String path) throws XMLStreamException, IOException {
		file = new File(path);
		SBMLDocument d = SBMLReader.read(file);
		model = d.getModel();

		numR = model.getNumReactions();
		numS = model.getNumSpecies();

		rct = new String[model.getNumReactions()];
		met = new String[model.getNumSpecies()];

		for (int i = 0; i < model.getNumReactions(); i++) {
			Reaction r = model.getReaction(i);
			String id = r.getId();
			rct[i] = id;
		}

		for (int i = 0; i < model.getNumSpecies(); i++) {
			Species s = model.getSpecies(i);
			String id = s.getId();
			met[i] = id;
		}

		// set objective function, set lower and upper bound
		objectiveFunction = new double[numR];
		lowerBound = new double[numR];
		upperBound = new double[numR];

		for (int i = 0; i < numR; i++) {
			Reaction r = model.getReaction(i);
			KineticLaw k = r.getKineticLaw();
			LocalParameter p1 = k.getLocalParameter("OBJECTIVE_COEFFICIENT");
			LocalParameter p2 = k.getLocalParameter("LOWER_BOUND");
			LocalParameter p3 = k.getLocalParameter("UPPER_BOUND");

			double v = p1.getValue();

			lowerBound[i] = p2.getValue();
			upperBound[i] = p3.getValue();

			if (v == 0.0) {
				objectiveFunction[i] = 0;
			} else {
				objectiveFunction[i] = 1;
				biomassOptValuePos = i;
			}

		}

		// test
		if (model != null && numR != 0 && numS != 0 && (rct.length > 0)
				&& (met.length > 0)) {
		} else {
			System.out.println("Error load data in loadSBML.");
			System.exit(1);
		}

	}

}
