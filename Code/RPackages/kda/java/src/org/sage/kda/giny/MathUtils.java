/*
  Copyright (c) 2010, Sage Bionetworks

  This library is free software; you can redistribute it and/or modify it
  under the terms of the GNU Lesser General Public License as published
  by the Free Software Foundation; either version 2.1 of the License, or
  any later version.

  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY, WITHOUT EVEN THE IMPLIED WARRANTY OF
  MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  The software and
  documentation provided hereunder is on an "as is" basis, and the
  Institute for Systems Biology and the Whitehead Institute
  have no obligations to provide maintenance, support,
  updates, enhancements or modifications.  In no event shall the
  Institute for Systems Biology and the Whitehead Institute
  be liable to any party for direct, indirect, special,
  incidental or consequential damages, including lost profits, arising
  out of the use of this software and its documentation, even if the
  Institute for Systems Biology and the Whitehead Institute
  have been advised of the possibility of such damage.  See
  the GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with this library; if not, write to the Free Software Foundation,
  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
*/
package org.sage.kda.giny;

import java.util.Arrays;

import org.apache.commons.math.special.Gamma;


public class MathUtils {
	
	public static double[] qvals(double[] obs, double[] fdr){
		
		int N = obs.length;
		assert(fdr.length == N);
		double[] qvals = new double[N];
		
		int[] sortIdxs = indexSort(obs);
		
		qvals[sortIdxs[N-1]] = fdr[sortIdxs[N-1]];
		for(int i = N-2; i > -1; --i){
			qvals[sortIdxs[i]] = Math.min(fdr[sortIdxs[i]],
											qvals[sortIdxs[i+1]]);
		}
		return qvals;
	}
	
	public static double[] fdr(double[] obs, double[] background){
		int[] obsIdxs = indexSort(obs);
		int[] permIdxs = indexSort(background);
		double[] fdr = new double[obs.length];
		int backgroundIdx = 0;
		for(int i = 0; i < fdr.length; ++i){
			int idx = findIdxOfNull(obs[obsIdxs[i]], backgroundIdx, background, permIdxs);
			double top = (((double)idx) + 1) / (background.length + 1);
			double bottom = ((double)i+1) / obs.length;
			fdr[obsIdxs[i]] = Math.min(top / bottom, 1.0); 
		}
		return fdr;
	}
	
	// utility function for fdr
	private static int findIdxOfNull(double val, int startIdx, double[] background, int[] backgroundIdxs){
		int idx = startIdx;
		do{
			if(val < background[backgroundIdxs[idx]])
				return idx;
		}while((idx++ < background.length));
		return backgroundIdxs.length;
	}
	
	private static int[] indexSort(double[] values){
		SortPair[] sortPair = new SortPair[values.length];
		for(int i = 0; i < sortPair.length; ++i){
			sortPair[i] = new SortPair(values[i], i);
		}
		Arrays.sort(sortPair);
		int[] sortedIdxs = new int[sortPair.length];
		for(int i = 0; i < sortedIdxs.length; ++i){
			sortedIdxs[i] = sortPair[i].originalIndex;
		}
		return sortedIdxs;
	}
	
	public static double mean(int[] values){
		double acc = 0;
		for(int i = 0; i < values.length; ++i){
			acc += values[i];
		}
		return acc / values.length;
	}
	
	/**
	 * Standard deviation - unbiased
	 * @param values
	 * @return
	 */
	public static double sd(int[] values){
		double mean = mean(values);
		return sd(values, mean);
	}
	
	/**
	 * Standard deviation - unbiased
	 * @param values
	 * @return
	 */
	public static double sd(int[] values, double mean){
		double acc = 0.0;
		double diff;
		for(int i = 0; i < values.length; ++i){
			diff = values[i] - mean;
			acc += (diff * diff);
		}
		return Math.sqrt(acc / (values.length-1));
	}
	
	 /**
     * constructor with as arguments strings containing numbers.
     *
     * @param x    number of red balls selected in samples
     * @param bigX number of samples
     * @param n    number of red balls
     * @param bigN number of all balls
     */
	public static double hyperUpperCumulative(int x, int bigX, int n, int bigN) {
		if(bigN < 2)
			return 1.0;
		
		double sum = 0;
		//mode of distribution, integer division (returns integer <= double result)!
		int mode = (bigX+1)*(n+1)/(bigN+2) ;
		if(x >= mode){
            int i = x ;
            while ((bigN - n >= bigX - i) && (i <= Math.min(bigX, n))) {	
                double pdfi = Math.exp(Gamma.logGamma(n+1)-Gamma.logGamma(i+1)-Gamma.logGamma(n-i+1) + Gamma.logGamma(bigN-n+1)-Gamma.logGamma(bigX-i+1)-Gamma.logGamma(bigN-n-bigX+i+1)- Gamma.logGamma(bigN+1)+Gamma.logGamma(bigX+1)+Gamma.logGamma(bigN-bigX+1)) ;
                sum = sum+pdfi;
                i++;
            }	
		}else{
	       int i = x - 1;
	       while ((bigN - n >= bigX - i) && (i >= 0)) {
	    	   double pdfi = Math.exp(Gamma.logGamma(n+1)-Gamma.logGamma(i+1)-Gamma.logGamma(n-i+1) + Gamma.logGamma(bigN-n+1)-Gamma.logGamma(bigX-i+1)-Gamma.logGamma(bigN-n-bigX+i+1)- Gamma.logGamma(bigN+1)+Gamma.logGamma(bigX+1)+Gamma.logGamma(bigN-bigX+1)) ;
                    sum = sum+pdfi;
                    i--;
            }	
                sum = 1-sum;
         }
         return sum;
	}
        
	
	static class SortPair implements Comparable<SortPair>
	{
	  private int originalIndex;
	  private double value;

	  public SortPair(double value, int originalIndex){
	    this.value = value;
	    this.originalIndex = originalIndex;
	  }

	  public int compareTo(SortPair o){
	    return Double.compare(value, o.getValue());
	  }

	  public int getOriginalIndex(){
	    return originalIndex;
	  }

	  public double getValue(){
	    return value;
	  }
	}

}
