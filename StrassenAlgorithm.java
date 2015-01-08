package strassen;

public class StrassenAlgorithm {
	public static void main(String[] args) {
		//2X2 matrix
		int[][] matrixA2 = {{1,2}, {3,4}};
		int[][] matrixB2 = {{5,6}, {4,7}};
		//4X4 matrix
		int[][] matrixA4 = {{1,2,3,4},{5,6,7,8}, {1,2,3,4}, {4,3,2,1}};
		int[][] matrixB4 = {{2,4,6,8}, {1,3,5,7}, {2,6,4,3}, {9,8,7,6}};
		//8X8 matrix
		int[][] matrixA8 = {{1,2,3,4,5,6,7,8}, {2,3,4,5,6,7,8,9}, {4,3,5,6,4,5,3,2}, {5,6,4,5,3,1,2,3},
				{4,5,6,3,4,2,3,4}, {1,7,3,6,5,8,9,5}, {1,2,3,4,5,6,7,2}, {4,2,5,6,4,4,4,4}};
		int[][] matrixB8 = {{4,5,6,3,4,2,3,4}, {1,7,3,6,5,8,9,5}, {1,2,3,4,5,6,7,2}, {4,2,5,6,4,4,4,4},
		{1,2,3,4,5,6,7,8}, {2,3,4,5,6,7,8,9}, {4,3,5,6,4,5,3,2}, {5,6,4,5,3,1,2,3}};
		
		//print out 2X2 matrix A and B
		System.out.println("2X2 Matrix A");
		System.out.println();
		printArray(matrixA2);
		System.out.println();
		System.out.println("2X2 Matrix B");
		System.out.println();
		printArray(matrixB2);
		System.out.println();
		//print out multiplications of A and B
		System.out.println("Matrix Multiplication of A and B");
		System.out.println();
		printArray(matrixMultiplication(matrixA2, matrixB2));
		System.out.println();
		System.out.println("Strassen Multiplication of A and B");
		System.out.println();
		printArray(Strassen(matrixA2, matrixB2));
		
		//print out 4x4 matrix A and B
		System.out.println();
		System.out.println("4X4 Matrix A");
		System.out.println();
		printArray(matrixA4);
		System.out.println();
		System.out.println("4X4 Matrix B");
		System.out.println();
		printArray(matrixB4);
		System.out.println();
		//print out multiplications of A and B
		System.out.println("Matrix Multiplication of A and B");
		System.out.println();
		printArray(matrixMultiplication(matrixA4, matrixB4));
		System.out.println();
		System.out.println("Strassen Multiplication of A and B");
		System.out.println();
		printArray(Strassen(matrixA4, matrixB4));
		
		//print out 8x8 matrix A and B
		System.out.println();
		System.out.println("8X8 Matrix A");
		System.out.println();
		printArray(matrixA8);
		System.out.println();
		System.out.println("8X8 Matrix B");
		System.out.println();
		printArray(matrixB8);
		//print out multiplications of A and B
		System.out.println();
		System.out.println("Matrix Multiplication of A and B");
		System.out.println();
		printArray(matrixMultiplication(matrixA8, matrixB8));
		System.out.println();
		System.out.println("Strassen Multiplication of A and B");
		System.out.println();
		printArray(Strassen(matrixA8, matrixB8));
		System.out.println();
	}
	/**Strassen matrix multiplication**/
	public static int[][] Strassen(int[][] A, int[][] B){
		int size = A.length; 
	    int[][] C = new int[size][size];//make the final array size
	    
	    if (size==1){ //base case
	    	C[0][0] = A[0][0]*B[0][0];
        }
	    else{
	    	//make new arrays for the sections of matrices A and B
	      	int[][] A11 = getQuadrant(A, 0, 0);
            int[][] A12 = getQuadrant(A, 0, size/2);
            int[][] A21 = getQuadrant(A, size/2, 0);
	        int[][] A22 = getQuadrant(A, size/2, size/2);
	        int[][] B11 = getQuadrant(B, 0, 0);
	        int[][] B12 = getQuadrant(B, 0, size/2);
	        int[][] B21 = getQuadrant(B, size/2, 0);
	        int[][] B22 = getQuadrant(B, size/2, size/2);
	            
	       //M1 = (A11 + A22)(B11 + B22)  
	       int [][] M1 = Strassen(add(A11, A22), add(B11, B22));
	       //M2 = (A21 + A22) B11
	       int [][] M2 = Strassen(add(A21, A22), B11);
	       //M3 = A11 (B12 - B22)
	       int [][] M3 = Strassen(A11, subtract(B12, B22));
	       // M4 = A22 (B21 - B11)
	       int [][] M4 = Strassen(A22, subtract(B21, B11));
	       //M5 = (A11 + A12) B22
	       int [][] M5 = Strassen(add(A11, A12), B22);
	       //M6 = (A21 - A11) (B11 + B12)
	       int [][] M6 = Strassen(subtract(A21, A11), add(B11, B12));
	       //M7 = (A12 - A22) (B21 + B22)
	       int [][] M7 = Strassen(subtract(A12, A22), add(B21, B22));
	    
	       //C11 = M1 + M4 - M5 + M7
	       int[][] C11add= add(M1, M4);
	       int [][] C11subtract = subtract(C11add, M5);
	       int[][] C11 = add(C11subtract, M7);
	       //C12 = M3 + M5
	       int [][] C12 = add(M3, M5);
	       //C21 = M2 + M4
	       int [][] C21 = add(M2, M4);
	       //C22 = M1 - M2 + M3 + M6
	       int[][] C22add= add(M1, M3);
	       int [][] C22subtract = subtract(C22add, M2);
	       int[][] C22 = add(C22subtract, M6);
	
	       //Combine the matrix
	       combine(C11, C, 0 , 0);
	       combine(C12, C, 0 , size/2);
	       combine(C21, C, size/2, 0);
	       combine(C22, C, size/2, size/2);
	    }
       return C;
	}
	/**Gets the quadrants of the matrix**/
	public static int[][] getQuadrant(int[][] matrix, int num1, int num2){
		int size = matrix.length;
		int quadrant[][] = new int[size/2][size/2];
		for(int i = 0, j = num1; i < quadrant.length; i++, j++){
            for(int k = 0, l = num2; k < quadrant.length; k++, l++){
            	quadrant[i][k] = matrix[j][l];
            }
		}return quadrant;
	}	
	/**Subtracts two matrices**/
	public static int[][] subtract(int[][] A, int[][] B){
		//make a new matrix for the subtracted values
		int size = A.length;
        int[][] sub = new int[size][size];
        //for the values in the matrices subtract them
        for (int i = 0; i < size; i++){
            for (int j = 0; j < size; j++){
                sub[i][j] = A[i][j] - B[i][j];
            }
        }
        //return the subtracted values matrix
        return sub;
	}
	/**Adds two matrices**/
	public static int[][] add(int[][] A, int[][] B){
			int size = A.length;
			//make a new matrix for the added values
	        int[][] add = new int[size][size];
	      //for the values in the matrices add them
	        for (int i = 0; i < size; i++){
	            for (int j = 0; j < size; j++){
	                add[i][j] = A[i][j] + B[i][j];
	            }
	        }
	        return add;
	}
	/**Combines two matrices**/
	public static void combine(int[][] small, int[][] big, int num1, int num2){
		//for values in the small matrix that do not exceed the length of the small matrix
		//put the values from the small matrix in the big matrix
		for(int i = 0, j = num1; i < small.length; i++, j++){ 
            for(int k = 0, l = num2; k < small.length; k++, l++){
            	big[j][l] = small[i][k];
            }
		}
	}
	/**regular matrix multiplication**/
	public static int[][] matrixMultiplication(int[][] A, int[][] B) {
	       int rowsA = A.length;
	       int columnsA = A[0].length;
	       int columnsB = B[0].length;
	       int[][] C = new int[rowsA][columnsB];
	       for (int i = 0; i < rowsA; i++) {
	           for (int j = 0; j < columnsB; j++) {
	               for (int k = 0; k < columnsA; k++) {
	                   C[i][j] = C[i][j] + A[i][k] * B[k][j];
	               }
	           }
	       }
	       return C;
	   }
	/**print array**/
	public static void printArray(int[][] array){
		//print all the values in the array for a row
		for(int i = 0; i<array.length; i++){
			for(int j = 0; j<array.length; j++){
				System.out.print(array[i][j] + " ");
			}
			//start a new line for each row
			System.out.println();
		}
	}
}
