# ADBMS_Project2
# Template supplied by: John Paparrizos
# Code Implementaions supplied by: Jonathon Ross, Benji Ofori, Gabe Richner

# Description:
This project focuses on testing the efficiency of various optimized binary search algorithms and band join operations. These algorithms include techniques like arithmetic optimization, SIMD (Single Instruction Multiple Data), bit masking, and loop unrolling. The project demonstrates how these optimizations can significantly enhance performance, which is critical for database operations.


# Instructions:
- make the executable with the command <make>
- run the executable with the command ./dbms N X Y Z R 
    - N: Number of elements in the array being searched
    - X: outer table band join 
    - Y: size of the output result. Depends on the number of successful matches between the inner and outer tables. Y is influenced by variable N and X.
    - Z: bound that defines how close outer value is to the inner value to be considered a match
    - R: Number of repeates for the experiment, Help toward computing the the average run time over the multiple results
- Example:
    - ./db5242 10000 50 200 5 20