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

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class SetUtils {
	
	public static boolean[] OR(boolean[] a1, boolean[] a2){
		if(a1.length != a2.length)
			throw new IllegalStateException("Logical arrays not of same size");
		boolean[] r = new boolean[a1.length];
		for(int i = 0; i < a1.length; ++i)
			r[i] = a1[i] || a2[i];
		return r;	                 
	}
	
	public static <T> boolean isSetInSets(Set<T> set, Collection<Set<T>> collection){
		for(Set<T> check : collection){
			if(check.containsAll(set))
				return true;
		}
		return false;
	}
	
	public static <T> int[] match(List<? extends T> list, Set<? extends T> col){
		List<Integer> matches = new ArrayList<Integer>();
		for(int i = 0; i < list.size(); ++i){
			if(col.contains(list.get(i)))
				matches.add(i);
		}
		int[] idxs = new int[matches.size()];
		for(int i = 0; i < idxs.length; ++i){
			idxs[i] = matches.get(i);
		}
		return idxs;
	}
	
	public static <T> Set<T> union(Collection<? extends T> set1, Collection<? extends T> set2){
		Set<T> union = new HashSet<T>(set1);
		union.addAll(set2);
		return union;
	}
	
	public static <T> Set<T> diff(Collection<? extends T> set1, Collection<? extends T> set2){
		Set<T> diff = union(set1, set2);
		diff.removeAll(intersect(set1, set2));
		return diff;
	}
	
	public static <T> Set<T> intersect(Collection<? extends T> set1, Collection<? extends T> set2){
		Set<T> intersection = new HashSet<T>(set1);
		intersection.retainAll(set2);
		return intersection;
	} 
}
