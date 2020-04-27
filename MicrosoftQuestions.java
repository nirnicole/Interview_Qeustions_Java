import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Stack;

//writen for JAVA8
public class MicrosoftQuestions {

	//how many 'a' chars we can insert without 3 a consecutive
	public static int maxAInserts(String S) {
		
		if(S==null)
		return -1;

		if(S.equals(""))
			return 2;

		int acounter = 0;
		int strSize = S.length();
		int result = 0;

		S.toLowerCase();

		for(int i=0;i<strSize;i++) {
			char c = S.charAt(i);
			if(c == 'a')
				{
				acounter++;
				if(acounter == 3) return -1;
				}
			else {
				result += 2;
				result -= acounter;
				acounter = 0;
			}
		}

		if(S.charAt(strSize - 1) != 'a')
				result+=2;
			else
			{
				result += 2;
				result -= acounter;
			}
			return result;
		}

	  //maximum value when inserting '5' digit
	public static int maxPossibleValue(int num) {

		final int K=5;
		char[] chars = ("" + num).toCharArray();
		    int index=0;

		    if(num>=0)
			    while(index<chars.length && (K+48)<(int)chars[index] ) index++;
		    else 
			{
			index=1;	//skip the minus
			    while(index<chars.length && (K+48)>(int)chars[index] ) index++;
			}  

		    char[] result = new char[chars.length+1];

		    for(int i=0,charsIndex=0 ;i<result.length;i++,charsIndex++)
		    {
			if(i==index)
				{
				result[i] = K+48;
				charsIndex--;
				}
			else result[i] = chars[charsIndex];
		    }

		    return Integer.parseInt(String.valueOf(result));
	}

	//K days after given one.
	public  static String dayOfTheWeek(String S, int K) {
		// write your code in Java SE 8

		 if(K<0 || K>500)
			 return "Error: invalid K value";

		 String[] week = new String[] {"Mon", "Tue", "Wed", "Thu", "Fri","Sat", "Sun"};
		 final int DAYS_WEEK=7;
		 boolean error=true;

		 for(String dayCheck : week)
			 if(dayCheck.equalsIgnoreCase(S))
				 error=false;

		 if(error)
			 return "Error: invalid Day value";

		 Map<String, Integer> map = new HashMap<>();

		 for( int day=0; day<week.length ; day++)
			 map.put(week[day], day);

		 int resultDay = (map.get(S) + K) % DAYS_WEEK;

		 return week[resultDay];
		}
	

    //Given an integer n, return any array containing n unique integers such that they add up to 0.
    public int[] sumZero(int n) {
        int[] res = new int[n];
        
        if(n==1)
            return res; //{0}
        
        int start= -1 * (n/2);  //int, so we get floor

        
        if(n%2==1)
        {            
            for(int i=0; i<res.length ; i++)
                res[i] = start++;
        }
        else{
            for(int i=0; i<res.length ; i++)
                if(start!=0)
                    res[i] = start++;
                else {
                    start++;
                    i--;
                }
             }
        
        return res;
    }
    
	static Node prevNode=null;

	static void connectLeafs(Node root)  
	{  
		
	    // If node is null, return  
	    if (root == null)  
	        return;  
	    
	    // If node is leaf node, print its data  
	    if (root.left == null && root.right == null)  
	    {  
	        if(prevNode!=null)
	        	prevNode.right = root;
	        prevNode=root;
	        return;  
	    }  
	
	    // If left child exists, check for leaf  
	    // recursively  
	    if (root.left != null)  
	    	connectLeafs(root.left);  
	    
	    // If right child exists, check for leaf  
	    // recursively  
	    if (root.right != null)  
	    	connectLeafs(root.right);  
	  
} 
	//  -5,-2,1,5,8,10,15
	//M=3 // 
	
	//Find set of m-elements with difference of any two elements is divisible by k
	//Given an array of n-positive integers and a positive integer k, find a set of exactly m-elements such that difference of any two element is equal to k.
	//Largest M-aligned Subset
	public static int[] largestMalignedSubset(int arr[], int m, int n) {
		
		/*algorithem:
		 *we will insert the values of the array to a 2 dim list that will contain 
		 */
		
        // syntax check
        if (null == arr || arr.length < Math.max(2, n)) 
            return null;
        
        Map<Integer, List<Integer>> map = new HashMap<>();
        for (int i : arr) 
        {
            int rem = i % m;
            
            if (rem < 0) 
                rem = m + rem;
            
            if (!map.containsKey(rem)) 
                map.put(rem, new ArrayList<Integer>());
            
            map.get(rem).add(i);
        }
        
        int[] result = {}; //we will take only the max, the empty group is the minimum case.
        
        for (Map.Entry<Integer, List<Integer>> entry : map.entrySet()) 
        	{
        	if (entry.getValue().size() <= n && entry.getValue().size() > result.length ) 
        		{
        	    result = new int[entry.getValue().size()];
        	    for(int i=0; i< entry.getValue().size() ; i++)
        	    	result[i] = entry.getValue().get(i);
        	    }
            }
        
        return result;
	}
	
	
	//Find set of m-elements with difference of any two elements is divisible by k
	//Given an array of n-positive integers and a positive integer k, find a set of exactly m-elements such that difference of any two element is equal to k.
	//**brute approache**
	public static int[] findLargestSetDiffdivisibleByM(int arr[], int m) {
		
		int[] res;
		int finalMax=0, finalIndex=-1;
		int currentMax=1, chainIndex=0;
		
		//algorithem - find all pairs, chain pairs, return elements appeared in chain.
		//algorithem 2 - sort, start with first index, update max index and max chain created
		
		/* Sort the elements */
        Arrays.sort(arr); //O(nlogn)

		for(int currentIndex=0; currentIndex<arr.length; currentIndex++) 
		{
			chainIndex=currentIndex;
			for(int indexToCompare= currentIndex+1; indexToCompare<arr.length ; indexToCompare++ )
			{
				if( (Math.abs(arr[indexToCompare]-arr[chainIndex])) % m == 0 )
						{
						chainIndex = indexToCompare;
						currentMax++;
						}

			}
			if(currentMax > finalMax) {
				finalMax = currentMax;
				finalIndex = currentIndex;
			}
			currentMax=1;
		}
		
		res = new int[finalMax];
		int i=0;
		res[i++]= arr[finalIndex];
		chainIndex = finalIndex;
		
		for(int indexToCompare=chainIndex+1; indexToCompare < arr.length; indexToCompare++)
		{
			if( (Math.abs(arr[indexToCompare]-arr[chainIndex])) % m == 0 )
				{
				chainIndex = indexToCompare;
				res[i++] = arr[chainIndex];
				}
		}
		
		return res;
	}
	
	//Minimum Moves to Obtain String Without 3 Identical Consecutive Letters
	public static int minMovesWithoutIdenticalLetters(String s) {

		if(s.length()<0 || s.length()>200000)
	    	return -1;
	    if(s.length()<3)
	    	return 0;
		
		int moves = 0;
	    for (int i = 0 ; i < s.length(); i++) {
	        int runLength = 1;
	        for (; i + 1 < s.length() && s.charAt(i) == s.charAt(i + 1); i++) {
	            runLength++;
	        }
	        moves += runLength / 3;
	    }
	    return moves;
	}
	
	
	//Longest Substring Without 3 Contiguous Occurrences of Letter
	public static String longestSubstringWithout3ContiguousOccurrencesofLetter(String s) {

		if(s.length()<3)
	    	return s;

	    for (int i = 0 ; i < s.length(); i++) {
	        int runLength = 1;
	        for (; i + 1 < s.length() && s.charAt(i) == s.charAt(i + 1) && runLength<3; i++) {
	            runLength++;
	        }
	      if(runLength==3)
	    	  return s.substring(0, i);
	    }
	    return s;
	}
	
	// Longest Semi-Alternating Substring
	public static int  longestSemiAlternatingSubstring(String s) {

		if(s.length()<1 || s.length()>200000)
	    	return -1;
		for(int i=0;i<s.length(); i++)
			if(s.charAt(i) != 'a' && s.charAt(i) != 'b')
				return -1;

		int max=0, lastIndex=0;
		
	    for (int i = 0 ; i < s.length(); i++) {	        
	    	
	        int runLength = 1;
	        for (; i + 1 < s.length() && s.charAt(i) == s.charAt(i + 1) && runLength<3; i++) 
	        {
	            runLength++;
	        }
	        if(runLength==3)
	        {
		    	if(max<i-lastIndex)
	    		{
	    		max = i-lastIndex;
	    		lastIndex = i-1; //(we reached 3 consec.  letters "aaa" so we need to go back one index "a->aa" and try again.)
	    		}
	        }
	    }
	    //if max never been reached the whole string can be semi alter
	    return max==0 ? s.length() : max;
	}
	
	//String Without 3 Identical Consecutive Letters
	public static String stringClearedOfAll3IdenticalConsecutiveLetters(String s) {

		if(s.length()<1 || s.length()>200000)
	    	{
			System.out.println("error: invalid string length.");
			return null;
	    	}
		if(s.length()<3)
	    	return s;
	
		StringBuilder bufferString = new StringBuilder();
		
	    for (int i = 0 ; i < s.length(); i++) {
	        int runLength = 1;
	        bufferString.append(s.charAt(i));
	        for (; i + 1 < s.length() && s.charAt(i) == s.charAt(i + 1); i++) {
	            runLength++;
	            if(runLength<3)
	            	bufferString.append(s.charAt(i));
	        }
	    }
	    
	    return bufferString.toString();
	}
	
	//Lexicographically smallest string formed by removing at most one character.
	public static String  lexicographicallySmallestString(String s) {
		
		if(s.length()<1)
			return s;
		
		int max=0,maxIndex=0;
		String sLower = s.toLowerCase();
		
		for(int i=0; i<sLower.length() ; i++)
			if((int)sLower.charAt(i)>max)
			{
				max = (int)sLower.charAt(i);
				maxIndex = i;
			}
		
		return sLower.substring(0, maxIndex) + sLower.substring(maxIndex+1);
		
	}
	
	//Numbers With Equal Digit Sum
	public	static int numbersWithEqualDigSum(int[] A) {
	
		int ERROR_VALUE=-1;
		int MAX_SUM = 81;	// 81 is the maximum sum of 10 digit number, and we only need to save one pair.
	    int MAX_NUM_VALUE = 1000000000;
		int MAX_N = 200000;

		if(A.length<2 || A.length> MAX_N)
			return ERROR_VALUE;
		
		//space complexity O(82 * 2) = O(1)   
		int[][] hashTable = new int[MAX_SUM+1][2];
	    
		//time complexity at most O(10*n) = O(n)   
		for(int num : A)
		{
			if(num<1 || num>MAX_NUM_VALUE)
				return ERROR_VALUE;
			
			int sum=0;
			int tempNum = num;
			//time complexity can reach O(10) dor a 10 digit number,  in other cases it will bee O(logk) (log base 10 of k digits)
			while(tempNum>0)
			{
				sum+=tempNum%10;
				tempNum /=10;
			}
			if(hashTable[sum][0]<num)
				{
				hashTable[sum][1]=hashTable[sum][0];
				hashTable[sum][0]=num;
				}
			else if(hashTable[sum][1]<num)
				hashTable[sum][1]=num;
			// else -> do nothing with num, it is irelevent as smaller then our max current pair		
		}
		
		int max=ERROR_VALUE;
		for(int i=0 ; i<hashTable.length ; i++)
		{
			if(hashTable[i][0]==0 || hashTable[i][1]==0) continue;
			if((hashTable[i][0] + hashTable[i][1]) > max)
				max = hashTable[i][0] + hashTable[i][1];
		}

		if(max==0)
			max=ERROR_VALUE;

	return max;
	}
  
  
	public static int[] doubleAndSort(int[] arr) {
		
		int[] res = new int[arr.length];
		int[] temp = new int[arr.length];
		
		//O(n)
		for(int i=0; i<arr.length ; i++)
		{
			if(i<arr.length-1 && arr[i]==arr[i+1])
				{
				temp[i]=arr[i]*2;
				arr[i+1]=0;
				}
			else temp[i]=arr[i];
			
			//ignoring arr[i]=0;+
		}
		
		int resIndex=0;
		for(int number : temp)
			if(number!=0)
				res[resIndex++] = number;
		
		return res;
	}



	public static int maxIndexAscendingValuesGap(int[] arr) {
		
		int i=0, j=arr.length-1;
				
		while( i<arr.length && i<j )
			{
			if(arr[i]<arr[j])
				return j-i;
			
			if(arr[i+1]<arr[i])
				i++;
			
			if(arr[i]<arr[j])
				return j-i;
			
			if(arr[j-1]>arr[j])
				j--;
			
			if(arr[i]<arr[j])
				return j-i;	
			}
		
		return -1;
	}

	
	public static int maxNetworkRank(int[] A, int[] B , int N) {
		
		if(A.length!=B.length || A.length<1 || N<1)
			{
			System.out.print("\nError: illegal arrays length");
			return -1;
			}
		
		int[] hashTable = new int[N+1];

		for(int i=0; i< A.length ; i++)
		{
			if(A[i]==B[i] || A[i]<1 || B[i]<1 || A[i]>N || B[i]>N )
			{
				System.out.print("\nError: illegal edge on index " + i);
				return -1;
			}
			hashTable[A[i]]++;
			hashTable[B[i]]++;
		}
		
		int max=0, viceMax=0;
		for(int i=0; i<hashTable.length; i++)
		{
			if(hashTable[i]>hashTable[max])
			{
				viceMax=max;
				max=i;
			}
			else if(hashTable[i]>hashTable[viceMax])
				viceMax = i;

		}

		for(int i=0; i< A.length ; i++)
			if( (A[i]==max && B[i]==viceMax) || (B[i]==max && A[i]==viceMax) )
				return hashTable[max] + hashTable[viceMax] -1;
		
		return hashTable[max] + hashTable[viceMax];
		
	}
	
	
	//Min Deletions to Make Frequency of Each Letter Unique
	public static int minDeletionstoMakeFrequencyofEachLetterUnique(String s) {
		
		String stringTemp = s.toLowerCase();
		//'a' = 97 in ascii  	, hash table will consist 0-25 places for every abc letter
		int[] abcHashingTable = new int[26];
		
		for(int index=0 ; index<stringTemp.length() ; index++)
			abcHashingTable[ (int)stringTemp.charAt(index) - 97]++;
		
		int[] accurenceTable = new int[s.length()+1];
		
		for(int index=0 ; index<abcHashingTable.length ; index++)
			if(abcHashingTable[index]>0)
				accurenceTable[abcHashingTable[index]]++;
		
		
		int counter=0;
		
		for(int i=accurenceTable.length-1 ; i>0 ; i--)
		{
			if(accurenceTable[i]>1)
			{
				accurenceTable[i-1]+= accurenceTable[i]-1;
				counter+= accurenceTable[i]-1;
				accurenceTable[i]=1;
			}
		}

		return counter;
	}
	
	/*
	 * aabbffddeaee

aaa
eee
bb
ff
dd

1	2	3
0	3	2
	
0	4	1	->1	
3	1	1	->3
1	1	1	->2
Sum=6
	 * */
	 
	
	//Min Adj Swaps to Make Palindrome
	public static int minAdjSwapstoMakePalindrome(String s) {
		
		String stringTemp = s.toLowerCase();
		//'a' = 97 in ascii  	, hash table will consist 0-25 places for every abc letter
		int[] hashingTable = new int[26];
		
		for(int index=0 ; index<stringTemp.length() ; index++)
			hashingTable[ (int)stringTemp.charAt(index) - 97]++;
		
		int counterOnes = 0;
		
		for(int i=0; i<hashingTable.length ; i++)
			if(hashingTable[i]%2==1)
				counterOnes++;
		
		//string cant be handeld as a polyndrom
		if(counterOnes>1)
			return -1;

		//sorting the polyndrom
		int l=0, r=stringTemp.length()-1, steps=0;
	
		char[] charedString = stringTemp.toCharArray();
		
		while(l<r)
		{
			
			if(charedString[l] == charedString[r]) {
			//do nothing
			}
			else if(hashingTable[ (int)charedString[r] - 97]==1)
			{
				int tempr = r;
				
				while(charedString[l] != charedString[tempr])
					tempr--;
				
				while(tempr<r) {
					char tempChar = charedString[tempr+1];
					charedString[tempr+1] = charedString[tempr];
					charedString[tempr] = tempChar;
					tempr++;steps++;
				}
			}
			//if l index is the uniq and in case their just different we go from the left side, it doesnt matter from where we go if its not the uniq.
			else {
				int templ = l;
				
				while(charedString[templ] != charedString[r])
					templ++;
				
				while(templ>l) {
					char tempChar = charedString[templ-1];
					charedString[templ-1] = charedString[templ];
					charedString[templ] = tempChar;
					templ--;steps++;
				}
			}
			//advance to next indexes after all been handeld
			l++;r--;

		}
				return steps;
	
	}//end of method
	
	
	//Max Possible Value
	public static int maxPossibleValue(int num) {
		
		final int K=5;
		char[] chars = ("" + num).toCharArray();
	    int index=0;
	    
	    if(num>=0)
		    while(index<chars.length && (K+48)<(int)chars[index] ) index++;
	    else 
	    	{
	    	index=1;
		    while(index<chars.length && (K+48)>(int)chars[index] ) index++;
	    	}  
	    
	    char[] result = new char[chars.length+1];
	    
	    for(int i=0,charsIndex=0 ;i<result.length;i++,charsIndex++)
	    {
	    	if(i==index)
	    		{
	    		result[i] = K+48;
	    		charsIndex--;
	    		}
	    	else result[i] = chars[charsIndex];
	    }
	    
	    return Integer.parseInt(String.valueOf(result));
	}
	
	
	//Count Visible Nodes in Binary Tree
	// Time: O(n)
	// Space: O(n)
	public int countVisibleNodes(TreeNode root) {
	    return countVisibleNodes(root, Integer.MIN_VALUE);
	}

	private int countVisibleNodes(TreeNode node, int maxSoFar) {
	    if (node == null) return 0;

	    int count = 0;

	    if (node.val >= maxSoFar) {
	        count = 1;
	        maxSoFar = node.val;
	    }

	    return count + countVisibleNodes(node.left, maxSoFar) + countVisibleNodes(node.right, maxSoFar);
	}
	
	public static class TreeNode {
        int val;
        TreeNode left;
        TreeNode right;

        TreeNode(int val) {
            this.val = val;
        }
	}
	
	//Largest number X which occurs X times
	public static int findTheMax(int[] nums) {
        HashMap<Integer, Integer> map = new HashMap<>(); 
        for (int num : nums) {
            map.put(num, map.getOrDefault(num, 0) + 1);
        }
        int max = 0;
        for (int key : map.keySet()) {
            if (key == map.get(key)) {
                max = Math.max(max, key);
            }
        }
        return max;
    }
	
	//Fair Indexes
	public static int getNumOfFairIndexes(int[] A, int[] B) {
		int res = 0, sumA = 0, sumB = 0;
		for(int i=0;i<A.length;i++) {
			sumA += A[i];
			sumB += B[i];
		}

	     if (sumA != sumB || sumA % 2 != 0 || sumB % 2 != 0)    // if total sum of arrays are not equal or not even, then can't divide
	        return 0;

		int tmpA = 0, tmpB = 0;
		for(int i=0;i<A.length-1;i++) {
			tmpA += A[i];
			tmpB += B[i];
			if(sumA == 2 * tmpA && tmpA == tmpB)     // Only need to check either sumA or sumB is twice of tmpA or tmpB
				res++;
		}
		return res;
	}
	
	// Microsoft | OA 2020 | Count Visible Nodes in Binary Tree
	// https://leetcode.com/discuss/interview-question/546703/
    
	    // or, you can simply count the number of visibles without this O(n) space
	    static ArrayList<Integer> visible = new ArrayList<>();
	    
	    public static void dfs(TreeNode node, int max) {
	        if (node == null) {
	            return;
	        }
	        
	        if (node.val >= max) {
	            visible.add(node.val);
	            max = Math.max(node.val, max);
	        }
	        
	        dfs(node.left, max);
	        dfs(node.right, max);
	    }

	    
	    // class that stores the value of count 
	    public static class count { 
	        int c = 0; 
	        int finalValue=0;
	    }
	    
	    public static int KthNodeInBST(TreeNode root, int k) {
	        count c = new count(); // object of class count 
	    	KthNodeInBST(root,k, c);
	    	return c.finalValue;
	    }
	    
	    private static void KthNodeInBST(TreeNode root,int k,  count index) {
	    	
	    	if(root==null|| index.c>=k)
	    		return ;
	 
	    	KthNodeInBST(root.right, k,index);

	    	index.c++;

	    	if(k==index.c)
	    		{index.finalValue = root.val; return;}
	    	
	    	KthNodeInBST(root.left,k, index);
  
	 	    }
	    
	    
	    public static int KthSmallestNodeInBST(TreeNode root, int k) {
	        count c = new count(); // object of class count 
	        KthSmallestNodeInBST(root,k, c);
	    	return c.finalValue;
	    }
	    
	    private static void KthSmallestNodeInBST(TreeNode root,int k,  count index) {
	    	
	    	if(root==null|| index.c>=k)
	    		return ;
	 
	    	KthSmallestNodeInBST(root.left, k, index);

	    	index.c++;

	    	if(k==index.c)
	    		{index.finalValue = root.val; return;}
	    	
	    	KthSmallestNodeInBST(root.right, k, index);
  
	 	    }
	
	    
        static Stack<TreeNode> stack = new Stack<TreeNode>();
    
        public static int kthSmallest(TreeNode root,int k) {
            
            myKthSmallestNodeInBST(root,k);
            TreeNode res = stack.pop();
            stack.clear();
            return res.val;
        }

	    public static void myKthSmallestNodeInBST(TreeNode root,int k) {
	    	
		    if(root==null || stack.size()>k)
	    		return ;
	 	      	       
	        myKthSmallestNodeInBST(root.left, k);
	        
	        if(stack.size()>=k)
	        	return;
	        else 	 stack.push(root);
	    	
	    	myKthSmallestNodeInBST(root.right, k);
  
	 	    }
	    
    
	    
	    
	    public static int[][] merge(int[][] intervals) {
	        /*
	         Input: [[1,3],[2,6],[8,10],[15,18]]
	         Output: [[1,6],[8,10],[15,18]]
	         Explanation: Since intervals [1,3] and [2,6] overlaps, merge them into [1,6].
	        */
	    	
	    	ArrayList<int[]> list = new ArrayList<int[]>(0);
	    	
			Arrays.sort(intervals, (i1, i2) -> Integer.compare(i1[0], i2[0]));

	    	
	    	for(int[] item : intervals)
	    		list.add(item);   	
	    	
	   
	    	/*
	    	Collections.sort(list, new Comparator<int[]>() {
	    	    @Override
	    	    public int compare(int[] o1, int[] o2) {
	    	        if(o1[0]>o2[0] || (o1[0]==o2[0] && o1[1]>o2[1]) || (o1[0]==o2[0] && o1[1]==o2[1]))
	    	        	return 1;
	    	        else return -1;
	    	    }
	    	});	   
	    	
	   
	    	System.out.println();
	    	for(int[] item : intervals)
	    		System.out.println(item[0]);
	    	System.out.println();
	    	
	    	for(int[] item : list)
	    		System.out.println(item[0]);
    	*/
	    	
	        for(int interval=0; interval<list.size()-1 ; interval++)
	        {
	        	if(list.get(interval)[1] >= list.get(interval+1)[0])
	        		{
	        		if(list.get(interval)[1] > list.get(interval+1)[1])
	        			list.get(interval)[1] = list.get(interval+1)[1];
	        		list.remove(interval+1);
	        		interval--;
	        		}
	        }
	     
	    	int[][] res = new int[list.size()][2];
	    	
	    	for(int interval=0; interval<list.size() ; interval++)
	    	{
	    		res[interval][0] = list.get(interval)[0];
	    		res[interval][1] = list.get(interval)[1];
	    	}
	    	
	    return res;
	    }
	    
	    /*
	     * "123"
	     * "22"
	     * 	0	2*1	2*2	2*3
	     * 	2*1	2*2	2*3	0
	     * 	2	4+2	4+6	6
	     * 2	6+1	0	6
	     * 2	7	0	6
	     * "2706"
	     * 
	     */

	    
        public static String multiply(String num1, String num2) {
            int[] arr = new int[num1.length()+num2.length()];
            for(int i = num1.length()-1; i >= 0; i--) {
                int carry = 0;
                for(int j = num2.length()-1; j >= 0; j--) {
                    arr[i+j+1] += carry + (num1.charAt(i)-'0') * (num2.charAt(j)-'0');
                    carry = arr[i+j+1] / 10;
                    arr[i+j+1] %= 10;
                }
                arr[i] = carry;
            }
            StringBuilder builder = new StringBuilder();
            int i = 0;
            while(i < arr.length && arr[i] == 0) i++;
            if(i >= arr.length) return "0";
            for(int j = i; j < arr.length; j++) {
                builder.append(arr[j]);
            }
            return builder.toString();
        }
        
        public static String multiplyLargeRec(int num1,int num2) {
        	if(num2==1)
        		return Integer.toString(num1);
        	if(num2==0)
        		return "1";
        	if(num2<0)
        		throw new RuntimeException("\nError:\nnot supporting negetive powering.\n");
        	
        	return multiply(Integer.toString(num1), multiplyLargeRec(num1, num2-1)); 
        }
	    
	    public static class Singleton {
	    	  private static Singleton INSTANCE=null;
	    	 
	    	  // Private constructor suppresses generation of a (public) default constructor
	    	  private Singleton() {}
	    	  	//	    	    	can throw new RuntimeException("\nError:\ncant create another instance, \nonly one instance of this singletone is allowed.\n");
	    	 
	    	  public static Singleton getInstance() {
	    	    if(INSTANCE==null)
	    	    	INSTANCE = new Singleton();

	    	    return INSTANCE;
	    	  }
	    	}
	    
	    public static int[] sortBinarArraySwappong(int[] arr) {
	    	
	    	int start=0, end=arr.length-1 ;
	    	while( start<end )
	    	{
	    		while(arr[start]==1)
	    			start++;
	    		while(arr[end]==0)
	    			end--;
	    		if(start<end)
	    			{
	    			arr[end]=0;
	    			arr[start]=1;
	    			}
	    	}
	    	
	    	return arr;
	    }
	    
		public static void insertNodeInSortedList(Node root, Node nodeToReplace, int place) {
		    if (nodeToReplace == null) throw new RuntimeException("Error: null node asked.");
		    Node preiorNode=root;
		    
		    while(place>1 && preiorNode!=null)
		    {
		    preiorNode = preiorNode.right;
		    place--;
		    }
		    
		    if(place>1 || preiorNode==null) throw new RuntimeException("Error: null nodes reached.");
		    
		    Node nextNode = preiorNode.right;
		    		    
		    //preiorNode.left = the same;
		    preiorNode.right = nodeToReplace;
		    nodeToReplace.left = preiorNode;
		    nodeToReplace.right = nextNode;
		    nextNode.left = nodeToReplace;
		    //nextNode.right = the same;
		
		}
	    
	    public static int[] findMaxConSubarray2(int[] arr)
	    {
	        
	        int max=0;
	        int sum=0;
	        int startIndex=0 , endIndex=0;
	        int[] res = new int[2];     // {0,0}
	        
	        while(endIndex<arr.length)
	         {
	            if(sum>=0)
	            {
	                sum+=arr[endIndex];

	                if(sum>max)
	                    {
	                	max=sum; 
	                    res[0] = startIndex;
	                    res[1] = endIndex;
	                    }
		            endIndex++;
	            }
	            else{
	                sum=0;
	                startIndex = endIndex;
	            }
	         }
	        
	        return res;
	    }
	    
}//end of class
