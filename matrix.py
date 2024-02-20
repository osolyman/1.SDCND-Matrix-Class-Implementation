import math
from math import sqrt
import numbers

def zeroes(height, width):
        """
        Creates a matrix of zeroes.
        """
        g = [[0.0 for _ in range(width)] for __ in range(height)]
        return Matrix(g)

def identity(n):
        """
        Creates a n x n identity matrix.
        """
        I = zeroes(n, n)
        for i in range(n):
            I.g[i][i] = 1.0
        return I

class Matrix(object):

    # Constructor
    def __init__(self, grid):
        self.g = grid
        self.h = len(grid)               # rows
        self.w = len(grid[0])            # columns

    #
    # Primary matrix math methods
    #############################
 
    def determinant(self):
        """
        Calculates the determinant of a 1x1 or 2x2 matrix.
        """
        if not self.is_square():
            raise(ValueError, "Cannot calculate determinant of non-square matrix.")
        if self.h > 2:
            raise(NotImplementedError, "Calculating determinant not implemented for matrices largerer than 2x2.")
        
        if self.h == 1:               # in 1x1 matrix, then we return that 1 element
            return self.g[0][0]
        
        elif self.h == 2:             # in 2x2 matrix, then we make the calculation ad - bc
            a = self.g[0][0]
            b = self.g[0][1]
            c = self.g[1][0]
            d = self.g[1][1]
            
            return a*d - b*c

    def trace(self):
        """
        Calculates the trace of a matrix (sum of diagonal entries).
        """
        if not self.is_square():
            raise(ValueError, "Cannot calculate the trace of a non-square matrix.")

        d_sum = 0    

        for r in range(self.h):          # we iterate just through the rows, because our matrix should be square
            d_sum += self.g[r][r]
        
        return d_sum

    def inverse(self):
        """
        Calculates the inverse of a 1x1 or 2x2 Matrix.
        """
        if not self.is_square():
            raise(ValueError, "Non-square Matrix does not have an inverse.")
        if self.h > 2:
            raise(NotImplementedError, "inversion not implemented for matrices larger than 2x2.")
   
        result = []
    
        # if we have just 1 element, we divide 1/element, if it's not 0
        if self.h == 1:
            if self.g[0][0] == 0:
                raise ValueError("Inverse does not exist for 0 value.")  
            result.append([1 / self.g[0][0]])
        
        # if our matrix is 2x2
        elif self.h == 2:
            
            # first we call the determinant method
            det = self.determinant()
            if det == 0:
                 raise ValueError('The matrix is not invertible.')
                    
            else:
                # if the determinant is not 0, we reshape the matrix like [[d, -b], [-c, a]]
                result = [[self.g[1][1], -self.g[0][1]], [-self.g[1][0], self.g[0][0]]]
            
            # now we iterate through the grid and multiply the reversed determinant with each element in the grid
            for i in range(len(result)):
                for j in range(len(result[0])):
                    result[i][j] = (1/det) * result[i][j]
         
        return Matrix(result)        

    def T(self):
        """
        Returns a transposed copy of this Matrix.
        """
        result = []
        
        # to transpose, we do the outer iteration in columns
        for c in range(self.w):
            new_row = []
            
            # then the inner iteration is in rows
            for r in range(self.h):
                new_row.append(self.g[r][c])
            
            result.append(new_row)
        
        return Matrix(result)
    
    def dot_product(self, vec1, vec2):
        """
        Calculates the dot product of 2 vectors.
        
        for example: A [1, 2, 3] * B [1, 2, 3]
        
                    result = 1*1 + 2*2 + 3*3 = 14
        """
        # first checking if both vectors have the same length
        if len(vec1) != len(vec2):
            raise ValueError("Vectors must have the same length for dot product. ")

        result = 0
        
        # now iterating through any of them
        for i in range(len(vec1)):
            
            # multiply the same index element in each vector together and updating the result
            result += vec1[i] * vec2[i]
        
        return result

    def is_square(self):
        return self.h == self.w

    #
    # Begin Operator Overloading
    ############################
    def __getitem__(self,idx):
        """
        Defines the behavior of using square brackets [] on instances
        of this class.

        Example:

        > my_matrix = Matrix([ [1, 2], [3, 4] ])
        > my_matrix[0]
          [1, 2]

        > my_matrix[0][0]
          1
        """
        return self.g[idx]

    def __repr__(self):
        """
        Defines the behavior of calling print on an instance of this class.
        """
        s = ""
        for row in self.g:
            s += " ".join(["{} ".format(x) for x in row])
            s += "\n"
        return s

    def __add__(self,other):
        """
        Defines the behavior of the + operator
        """
        if self.h != other.h or self.w != other.w:
            raise(ValueError, "Matrices can only be added if the dimensions are the same") 
        
        result = []
        
        # iterate through each number in the any of both matrices and adding the elements with same index together
        for r in range(self.h):
            new_row = []
            for c in range(self.w):
                new_row.append(self.g[r][c] + other.g[r][c])
            result.append(new_row)
        
        return Matrix(result)

    def __neg__(self):
        """
        Defines the behavior of - operator (NOT subtraction)

        Example:

        > my_matrix = Matrix([ [1, 2], [3, 4] ])
        > negative  = -my_matrix
        > print(negative)
          -1.0  -2.0
          -3.0  -4.0
        """
        
        result = []
        
        # the same idea like add but here we just assign our elements with -, we don't do subtraction
        for r in range(self.h):
            new_row = []
            for c in range(self.w):
                new_row.append(-self.g[r][c])
            result.append(new_row)
            
        return Matrix(result)

    def __sub__(self, other):
        """
        Defines the behavior of - operator (as subtraction)
        """
        
        if self.h != other.h or self.w != other.w:
            raise(ValueError, "Matrices can only be subtracted if the dimensions are the same") 
            
        result = []
        
        # the same idea like add but we subtract instead of adding
        for r in range(self.h):
            new_row = []
            for c in range(self.w):
                new_row.append(self.g[r][c] - other.g[r][c])
            result.append(new_row)
        
        return Matrix(result)

    def __mul__(self, other):
        """
        Defines the behavior of * operator (matrix multiplication)
        """
        if self.w != other.h:
            raise ValueError("Number of columns in matrix 1 must match the number of rows in matrix 2 for multiplication. ")
        
        result = []
        
        # using transpose to multiply two matrices
        transpose_other = other.T()

        for r in range(self.h):
            new_row = []
            for c in range(other.w):
                dp = self.dot_product(self.g[r], transpose_other[c])
                new_row.append(dp)
            result.append(new_row)
            
        return Matrix(result)

    def __rmul__(self, other):
        """
        Called when the thing on the left of the * is not a matrix.

        Example:

        > identity = Matrix([ [1,0], [0,1] ])
        > doubled  = 2 * identity
        > print(doubled)
          2.0  0.0
          0.0  2.0
        """
        if isinstance(other, numbers.Number):
            result = []

            for r in range(self.h):
                new_row = []
                for c in range(self.w):
                    new_row.append(other * self.g[r][c])
                result.append(new_row)
            return Matrix(result) 
            