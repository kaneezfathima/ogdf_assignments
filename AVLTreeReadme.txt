Environment used:
Ubuntu 12.04 (precise) 64-bit
RAM: 15.4 GiB
Processor: Intel® Core™ i5-3320M CPU @ 2.60GHz × 4 
Compiler: gcc version 4.6.3 (Ubuntu/Linaro 4.6.3-1ubuntu5) 

Instructions to Compile and Run:
g++ avl_template.cpp -std=c++0x
./a.out

Summary:
A demo test case runs with creation of tree with ~30 random integers into an AVL tree
and displays the tree on console. Resultant tree after each of 2 insertions and 2 removals
is shown to demonstrate the balancing of tree.
After that Performance tests run with about 10^5, 10^6, 10^7 insertions,removals, lookups of random
sequence of integers into avltree, this is compared to insertions, removals, lookups in Red Black Tree
(stl map is used). It can be seen that insert, deletes are faster in red black trees as number of operations
increase, whereas lookups are faster in avl tree because avl tree is more rigidly balanced.
i.e AVL tree height is ~1.44logn whereas Red Black tree is ~2logn
To prevent running out of memory, elements inserted and deleted range between 1 to 10000, so a tree,
map would contain atmost 10000 elements at any instant.
Raw pointers are used again, since shared_ptr, unique_ptr detoriate performance.
