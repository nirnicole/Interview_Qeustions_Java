import java.util.Arrays;

public class DynamicProgramming {
	
	//dynamic programming - memoized
	public static int coinChange(int[] coins, int amount) {
	    if (amount < 1) return 0;
	    return coinChange(coins, amount, new int[amount]);
	  }

	  private static int coinChange(int[] coins, int rem, int[] count) {
	    if (rem < 0) return -1;
	    if (rem == 0) return 0;
	    if (count[rem - 1] != 0) return count[rem - 1];
	    int min = Integer.MAX_VALUE;
	    for (int coin : coins) {
	      int res = coinChange(coins, rem - coin, count);
	      if (res >= 0 && res < min)
	        min = 1 + res;
	    }
	    if(min == Integer.MAX_VALUE) 
	    		min=-1;
	    
	    count[rem - 1] = min;
	  
	    return count[rem - 1];
	        }
	
	  
	  //fibonacci with DP -- memoized
	  public static int fibonacci(int n) {
	  
		  int[] history = new int[n+1];
		  
		  return fibonacci(n , history);
	  }
	
	  private static int fibonacci(int n, int[] history) {
	  	
	  	int res=0;
	  	
	  	if(history[n]!=0)
	  		return history[n];
	  	
	  	if(n==1 || n==2)
	  		res = 1;
	  	else
	  		res = fibonacci(n-1,history) + fibonacci(n-2,history);
	  	
	  	history[n] = res;
	  	
	  	return res;
	  }
	
	  
	  //fibonacci with DP -- Bottom-UP (***NOT RECURSION***)
	      public static int fibonacciBU(int n) {
	      
	    	int[] history = new int[n+1];
	    	  
	      	if(n==1 || n==2)
	      		return 1;
	      		      	
	      	history[1] = 1;
	      	history[2] = 1;
	      	
	      	for(int i =3 ; i<=n ; i++)
	      		history[i] = history[i-1] + history[i-2];
	      		
	      	return history[n];
	      }
	      
	      
	      //
	      public static int solution2(int[] A) {
	    	  
	          if(A.length==1)
	              return Math.abs(A[0]);
	              
	          int[] absArray = new int[A.length];
	          
	          for(int i=0 ; i<A.length ; i++)
	        	  absArray[i] = Math.abs(A[i]);
	          
	          Arrays.sort(absArray);
	          
	          for(int i=0 ; i<A.length ; i++)
        	      	  System.out.println(absArray[i]);
	          
	          int sum=0;
	          int r=0, l = absArray.length-1;
	          
	          while(r<l)
	          {
	        	  System.out.println("\n r= " + absArray[r] + " l = " + absArray[l]  + "\t\t sum is: " + sum);
	        	  if(sum<absArray[l])
	        		  sum+=absArray[r++];
	        	  else 
	        		  sum-=absArray[l--];
	        	  
	        	  if(r==l)
	        		  sum-=absArray[l--];
	          }
	   
		      return sum;
	      }
	      	
	}
