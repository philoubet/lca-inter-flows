This is a sparse representation of the technology (A), intervention (B) and characterization (C) matrices.  This means that only the non-zero elements are stored.  Each row is a non-zero coefficient in the matrix.  The first two fields are the row and column position of the element, and the third field is the value of the element.  The row and column positions are zero-based.  

Examples:
row		column		coefficient
0		3039		-0.197490092
The coefficient at position (0, 3039) is equal to -0.197490092.  

For the meaning of the row and column indexes, consult the following files:
A_public, rows: ie_index.csv
A_public, columns: ie_index.csv
B_public, rows: ee_index.csv
B_public, columns: ie_index.csv
C, rows: LCIA_index.csv
C, columns: ee_index.csv

There is no uncertainty information in the C matrix because none is not provided by LCIA method developers.  
In the A matrix, off-diagonal coefficients have the opposite sign compared to what is found in the spold files.  The uncertainty information for triangular, uniform and undefined distribution has been adjusted accordingly.  
The lognormal distribution does not accept negative values.  However, the positive average can be entered, with the same variance as specified in the exchange, and sampling multiplied by -1 after its creation.  