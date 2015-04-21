package be.iminds.cloudspeller.postprocessing_Single;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.PriorityQueue;


import be.iminds.cloudspeller.output.MotifBLSRestrictionsWithCutoff;

public class BestMotifContainer {

	private int k;
	private MotifBLSRestrictionsWithCutoff restrictions;
	private PriorityQueue<GraphFmax> best = new PriorityQueue<GraphFmax>();
	
	
	public BestMotifContainer(int k, MotifBLSRestrictionsWithCutoff restrictions) {
		this.restrictions=restrictions;
		this.k=k;
	}

	public List<String> getAll() {
		
		List<String> l = new ArrayList<String>();
		Iterator<GraphFmax> iterator = best.iterator();
		while (iterator.hasNext()){
			l.add(iterator.next().getGraph());
		}
		return l;
	}

	public List<String> extractAll() {

		List<String> l = new ArrayList<String>();

		while (best.size()>0){
			l.add(best.poll().getGraph());
		}
		return l;
	}



	public void add(String g) {
		
		int F = restrictions.getMaxFamilies(g);
		
		if (F>0){
			best.offer(new GraphFmax(g,F));
		
			if (best.size()>k){
				best.poll();
			}
		}
	}

}

class GraphFmax implements Comparable<GraphFmax> {

	private int F;
	private String g;
	
	public GraphFmax(String g, int F){
		this.F=F;
		this.g=g;
	}
	
	public String getGraph() {
		return g;
	}

	@Override
	public int compareTo(GraphFmax o) {
		
		return this.F-o.F;
	}
	
}


