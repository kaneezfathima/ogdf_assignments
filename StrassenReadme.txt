Environment used:
Ubuntu 12.04 (precise) 64-bit
RAM: 15.4 GiB
Processor: Intel® Core™ i5-3320M CPU @ 2.60GHz × 4 
Compiler: gcc version 4.6.3 (Ubuntu/Linaro 4.6.3-1ubuntu5) 

Instructions to Compile and Run
g++ strassenMultiplication.cpp -L OGDF-snapshot/_release -lOGDF -I OGDF-snapshot/include -l pthread -std=c++0x
./a.out

Summary:

Runs a simple 4 X 4 matrix multiplication and prints result to output.

Function multiplyArray2D demonstrates multiplication of 2 Array2D objects.

Also runs a set of testcases doing multiplication of 2 random matrices 
of size 2^n X 2^n where n=1,2..10. Results are compared against normal 
multiplication, performance in terms of time taken is output too.
It can be noticed that strassen multiplication outperforms normal O(N^3) matrix
multiplication when N increases beyond 2^9. For smaller matrices there is no 
noticeable difference due to Strassen Algorithm.

Notes:
Matrix is represented as a Tree with 4 child nodes for internal node, 
(each representing submatrices). Recursive multiplication of submatrices is done,
with leaf node storing matrices of size <= Threshold provided. 
This way each element of matrix is stored only once and memory used, time taken is
minimised.

Templates are used throughout to allow for different types like int, float, long.

Using c++11 smart pointers like unique_ptr, shared_ptr was explored..
But performance detoriated by more than 5X due to additional book keeping, data used
by smart pointers. So only raw pointers have been used.

Attempted to follow the naming conventions, programming style of OGDF.

