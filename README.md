# ADBMS_Project2
# Template supplied by: John Paparrizos
# Code Implementaions supplied by: Jonathon Ross, Benji Ofori, Gabe Richner

# Description:
This project focuses on testing the efficiency of various optimized binary search algorithms and band join operations. These algorithms include techniques like arithmetic optimization, SIMD (Single Instruction Multiple Data), bit masking, and loop unrolling. The project demonstrates how these optimizations can significantly enhance performance, which is critical for database operations.


# Instructions:
- make the executable with the command <make>
- run the executable with the command ./dbms5242 N X Y Z R 
    - N: Number of elements in the array being searched
        This variable will specify the size of the data array that will be used for binary search operations
    - X: Size of the out table for band join
        This variable defines how many elements are in the outer array, these elements will be used in the band join operation.
    - Y: size of the output result. 
        This variable will have a varying number of successful matches between inner and outer tables. The result of the band join will be (0 < number <= Y).
        This output of Y is determined by the size of (N) and (X) and the bound of (Z)
    - Z: bound that defines how close the inner and outer values have to be for a match
        for a match to occur (Z > outer[i] - data[j]) Otherwise, its a miss.
        A higher bound increaes the change the two elements are withing the bound and causing a match.
    - R: Number of repeates for the experiment
         This helps toward computing the the average run time over the multiple results to help create more reliable metrics.
- Example:
    - ./db5242 10000 50 200 5 20
        N = 10000 random data array elements
        X = 50 elements in outer array
        Y = 200 maximum number of results to output
        Z = the bound between outer and data array has to be within 5.
        R = The program will run 20 times and calculate the average time.