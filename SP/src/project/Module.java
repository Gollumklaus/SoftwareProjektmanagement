package project;

import org.sbml.jsbml.Reaction;

import java.util.*;

/**
 * Represents a single flux module.
 * A flux module consists of a set of reactions and also has a fixed interface flux.
 * 
 * @author arne
 *
 */
public interface Module {

	/**
	 * Return the reactions in the module
	 * @return
	 */
	public Set<Reaction> getReactions();
	
	
	/**
	 * Return the interface flux vector.
	 * @return
	 */
	public Map<Reaction, Double> getInterface();
}
