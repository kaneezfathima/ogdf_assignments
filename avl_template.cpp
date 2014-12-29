/*
 * $ Author:Kaneez Fathima $
 * $ Date: 2014-03-15 $
 *
 * This file contains a templatised implementation of AVL tree.
 * Each node in Tree contains a Key, Value and link to left subtree, right subtree
 * This way the tree can be used for storing <key, value> pairs. A compare function
 * can be passed for custom data types. Default compare uses '<','==' operators to compare
 * 2 keys. User can choose to either overload '<', '==' operators or pass in a compare functor
 * during AvlTree instance creation. Test cases are added comparing avl tree insertion,
 * removals, lookups with stl's map(which uses Red Black Trees).
*/
#include<iostream>
#include<cmath>
#include<cstdlib>
#include<random>
#include<ctime>
#include<exception>
#include<iomanip>
#include<deque>
#include<sstream>
#include<map>
using namespace std;
/*
 *Default compare functor
*/
template < typename T1> struct defaultCompare{
	int operator () (const T1 &a, const T1 &b){
		if (a == b)
		{
			return 0;
		}
		else
		{
			return (a < b ? -1 : 1);
		}
	}
};

/*
* Custom exception to be thrown when we run out of memory,
* similar to OGDF Out of memory macro assertion.
*/
class OutOfMemoryException: public exception
{
	virtual const char* what() const throw()
	{
		return "Out of memory";
	}
} OOM;

/*
 * Structure of AVLNode<T1,T2,T3>
*/
template <typename T1,typename T2> 
struct AVLNode{
	T1 m_key;
	T2 m_val;
	int m_ht;
	AVLNode < T1,T2> *rlink;
	AVLNode < T1,T2> *llink;
};

/*
 *AVL Tree Class
*/
template < typename T1,typename T2,typename T3=defaultCompare<T1> >
class Avl{
	public:
		Avl ();
		~Avl ();
		//! Inserts a key and value into AVL tree.
		AVLNode < T1,T2 > *insert (const T1 & m_key,const T2 & value);
		//! Deletes a node with given key.
		void remove (const T1 & m_key);
		//! Search for a node with given key m_key.
		AVLNode < T1,T2 > *lookup (T1 & m_key);
		//! Prints Tree in order.
		void print();

	private:
		//! Inserts key, value when a root is given and returns root.
		AVLNode < T1,T2 > *insertNode (AVLNode < T1,T2 > *position, const T1 & m_key,
			const T2 & value, AVLNode < T1,T2 >* &result);
		//! Prints subtree rooted with p in order.
		void print (AVLNode < T1,T2 > *p, int indent);
		//! Left Rotates subtree rooted with p.
		AVLNode < T1,T2 > *leftRotate (AVLNode < T1,T2 > *p);
		//! Searches for m_key in a tree with a given root and returns node if found else nullptr ;
		AVLNode < T1,T2 > *search (AVLNode < T1,T2 > *root, const T1 & m_key);
		//! Deletes a node from a tree with given key m_key if found. 
		AVLNode < T1,T2 > *deleteNode (AVLNode < T1,T2 > *root, const T1 & m_key);
		//! Deletes the entire tree rooted at root. Used by destructor.
		void deleteTree (AVLNode < T1,T2 > *root);
		//! Get Balance factor of a given node.
		int balance (AVLNode < T1,T2> *p);
		//! To find maximum of 2 integers.
		int maximum (const int &a,const int &b);
		//! Get height of tree.
		int getHeight (AVLNode < T1,T2 > *p);
		//! Finds the min value in tree, i.e leftmost node;
		AVLNode < T1,T2 > *minValueNode (AVLNode <T1,T2 > *p);
		//! Right rotates subtree rooted with p
		AVLNode < T1,T2 > *rightRotate (AVLNode < T1,T2 > *p);
		T3 cmp;
		AVLNode < T1,T2> *root;
};

template < typename T1,typename T2,typename T3 > 
Avl < T1,T2,T3 >::Avl (){
	root = nullptr;
}

template < typename T1,typename T2,typename T3 > 
Avl < T1,T2,T3 >::~Avl (){
	deleteTree(root);
}

template < typename T1,typename T2,typename T3 > 
void Avl < T1,T2,T3 >::deleteTree(AVLNode < T1,T2 > *root){
	if (root==nullptr) 
	{
		return;
	}
	deleteTree(root->llink);
	deleteTree(root->rlink);
	delete root;
}

template < typename T1,typename T2,typename T3 > 
void Avl < T1,T2,T3 >::print(){
	print(root);
}

/*
 * Searches for a AVLNode<T1> in the tree with the given m_key;
 * if no such AVLNode<T1> exists, nullptr shall be returned
*/
template < typename T1,typename T2,typename T3 > 
AVLNode < T1,T2> *Avl < T1,T2,T3 >::lookup (T1 & m_key){
	return search (root, m_key);
}


template < typename T1,typename T2,typename T3 >
AVLNode < T1,T2 > *Avl < T1,T2,T3 >::search (AVLNode < T1,T2 > *current,
	const T1 & m_key){
	if (current == nullptr)
	{
		return nullptr;
	}
	if (Avl<T1,T2,T3>:: cmp( m_key, current->m_key)==0)	
		//!m_key == current->m_key
	{
		return current;
	}
	else if (Avl<T1,T2,T3>:: cmp(m_key,current->m_key)==(-1))
		//!m_key < current->m_key
	{
		search (current->llink, m_key);
	}
	else if (Avl<T1,T2,T3>:: cmp(m_key,current->m_key)==(1))
		//!m_key > current->m_key
	{
		search (current->rlink, m_key);
	}
}

template < typename T1,typename T2,typename T3 > 
AVLNode < T1,T2 > *Avl < T1,T2,T3 >::insert (const T1 & m_key,const T2 & value){
	AVLNode < T1,T2 > *result;
	root = insertNode (root, m_key, value, result);
	return result;
}

template < typename T1,typename T2,typename T3 > 
int Avl < T1,T2,T3 >::getHeight (AVLNode < T1,T2> *p)
{
	if (p == nullptr)
	{
		return 0;
	}
	return p->m_ht;
}

template < typename T1,typename T2,typename T3 > 
int Avl < T1,T2,T3 >::maximum (const int &a, const int &b)
{
	return (a > b) ? a : b;
}

template < typename T1,typename T2,typename T3 > 
int Avl < T1,T2,T3 >::balance (AVLNode < T1,T2 > *p)
{
	if (p == nullptr)
	{
		return 0;
	}
	return (getHeight (p->llink) - getHeight (p->rlink));
}

template < typename T1,typename T2,typename T3 > 
AVLNode < T1,T2> *Avl < T1,T2,T3 >::leftRotate (AVLNode < T1,T2 > *p)
{
	AVLNode < T1,T2> *q = p->rlink;
	AVLNode < T1,T2> *z = q->llink;
	q->llink = p;
	p->rlink = z;
	//update heights
	p->m_ht = maximum (getHeight (p->llink), getHeight (p->rlink)) + 1;
	q->m_ht = maximum (getHeight (q->llink), getHeight (q->rlink)) + 1;
	return q;
}

template < typename T1,typename T2,typename T3 >
AVLNode < T1,T2> *Avl < T1,T2,T3 >::rightRotate (AVLNode < T1,T2 > *p)
{
	AVLNode < T1,T2> *q = p->llink;
	AVLNode < T1,T2> *z = q->rlink;
	q->rlink = p;
	p->llink = z;
	//update heights
	p->m_ht = maximum (getHeight (p->llink), getHeight (p->rlink)) + 1;
	q->m_ht = maximum (getHeight (q->llink), getHeight (q->rlink)) + 1;
	return q;
}

template < typename T1,typename T2,typename T3 >
AVLNode < T1,T2> *Avl < T1,T2,T3 >::insertNode (AVLNode < T1,T2> *current,
	const T1 &m_key, const T2 &value, AVLNode < T1,T2>* &result) {
	if (current == nullptr)
	{
		//create a new AVLNode<T1>
		AVLNode < T1,T2 > *newNode = new AVLNode < T1,T2 > ();
		if ( newNode == 0) 
		{
			throw OOM;
		}
		newNode->m_key = m_key;
		newNode->m_val=value;
		newNode->llink = nullptr;
		newNode->rlink = nullptr;
		newNode->m_ht = 1;
		result = newNode;
		return newNode;
	}
	if (Avl<T1,T2,T3>:: cmp( m_key, current->m_key)==0)
		//! cmp(current->m_key,m_key ) ==0
	{
		current->m_val=value;
		result = current;
		return current;
	}

	if (Avl<T1,T2,T3>:: cmp(m_key,current->m_key)==(-1))
		//! m_key < current->m_key
	{
		current->llink = insertNode (current->llink, m_key, value, result);
	}
	else
	{
		current->rlink = insertNode (current->rlink, m_key, value, result);
	}
	current->m_ht =
		maximum (getHeight (current->llink), getHeight (current->rlink)) + 1;
	int difference = balance (current);
	//left left
	if( (difference > 1) && (Avl<T1,T2,T3>:: cmp(m_key,current->llink->m_key)==(-1)))
		//!m_key < current->llink->m_key
		return rightRotate (current);
	//right right
	if ((difference < -1) && (Avl<T1,T2,T3>:: cmp(m_key,current->rlink->m_key)==(1)))
		//! m_key > current->rlink->m_key
		return leftRotate (current);
	//left right
	if ((difference > 1) && (Avl<T1,T2,T3>:: cmp(m_key,current->llink->m_key)==(1)))
		//! m_key > current->llink->m_key
	{
		current->llink = leftRotate (current->llink);
		return rightRotate (current);
	}
	//right left
	if ((difference < -1) && (Avl<T1,T2,T3>:: cmp(m_key,current->rlink->m_key)==(-1)))
		//! m_key < current->rlink->m_key
	{
		current->rlink = rightRotate (current->rlink);
		return leftRotate (current);
	}
	return current;
}

template < typename T1,typename T2,typename T3 > 
void Avl < T1,T2,T3 >::remove (const T1 & m_key){
	if (root == nullptr)
	{
		return;
	}
	root = deleteNode (root, m_key);
	return;
}

template < typename T1,typename T2,typename T3 >
AVLNode < T1,T2> *Avl < T1,T2,T3 >::minValueNode (AVLNode < T1,T2> *p){
	AVLNode < T1,T2> *current = p;
	while (current->llink != nullptr)
	{
		current = current->llink;
	}
	return current;
}

template < typename T1,typename T2,typename T3 >
AVLNode < T1,T2> *Avl < T1,T2,T3 >::deleteNode (AVLNode < T1,T2> *root, const T1 & m_key){
	if (root == nullptr)
		return root;
	//! If the m_key to be deleted is smaller than the root's m_key,
	// then it lies in left subtree
	if (Avl<T1,T2,T3>:: cmp(m_key,root->m_key)==(-1))
		//! m_key < root->m_key
		root->llink = deleteNode (root->llink, m_key);
	//! If the m_key to be deleted is greater than the root's m_key,
	//! then it lies in right subtree
	else if (Avl<T1,T2,T3>:: cmp(m_key,root->m_key)==(1))
		//! m_key > root->m_key
		root->rlink = deleteNode (root->rlink, m_key);
	//! If m_key is same as root's m_key, then This is the AVLNode<T1>
	//! to be deleted
	else
	{
		// AVLNode<T1> with only one child or no child
		if ((root->llink == nullptr) || (root->rlink == nullptr))
		{
			AVLNode < T1,T2> *temp = root->llink ? root->llink : root->rlink;
			// No child case
			if (temp == nullptr)
			{
				temp = root;
				root = nullptr;
			}
			else{
				// One child case
				// Copy the contents of the non-empty child
				root->m_key=temp->m_key;
				root->m_val=temp->m_val;
				root->llink=temp->llink;
				root->rlink=temp->rlink;
			}			
			delete temp;
		}
		else{
			AVLNode < T1,T2> *temp = minValueNode (root->rlink);
			//! Copy the inorder successor's data to this AVLNode<T1>*/
			root->m_key = temp->m_key;
			root->m_val=temp->m_val;
			// Delete the inorder successor
			root->rlink = deleteNode (root->rlink, temp->m_key);
		}
	}
	//! If the tree had only one AVLNode<T1> then return
	if (root == nullptr)
		return root;
	//! STEP 2: UPDATE HEIGHT OF THE CURRENT NODE
	root->m_ht = maximum (getHeight (root->llink), getHeight (root->rlink)) + 1;
	//! STEP 3: GET THE BALANCE FACTOR OF THIS NODE
	//! (to check whether this AVLNode<T1> became unbalanced)
	int difference = balance (root);
	//! If this AVLNode<T1> becomes unbalanced, then there are 4 cases.
	// Left Left Case
	if (difference > 1 && balance (root->llink) >= 0){
		return rightRotate (root);
	}
	// Left Right Case
	if (difference > 1 && balance (root->llink) < 0){
		root->llink = leftRotate (root->llink);
		return rightRotate (root);
	}
	// Right Right Case
	if ((difference < (-1)) && balance (root->rlink) <= 0)
		return leftRotate (root);
	// Right Left Case
	if ((difference < (-1)) && balance (root->rlink) > 0)
	{
		root->rlink = rightRotate (root->rlink);
		return leftRotate (root);
	}
	return root;
}

template < typename T1,typename T2,typename T3 >
void Avl < T1,T2,T3 >::print (AVLNode < T1,T2> *p, int indent = 0){
	if (p != nullptr)
	{
		if (p->llink)
		{
			print (p->llink, indent + 4);
		}
		cout << p->m_key << ":";
		cout<< p->m_val<<"\t"<<endl;
		if (p->rlink)
		{
			print (p->rlink, indent + 4);
		}
	}
}

//! Test functor to compare doubles.
struct doubleCompare{
	int operator () (const double &a, const double &b){
		if (abs(a - b)<0.0001)
		{
			return 0;
		}
		else
		{
			return (a < b ? -1 : 1);
		}
	}
};

// Test utility function. Converts an integer value to string.
string intToString(int val) {
	ostringstream ss;
	ss << val;
	return ss.str();
}

// Test function used to print avl tree. Print the arm branches (eg, /    \ ) on a line.
void printBranches(int branchLen, int nodeSpaceLen, int startLen,
		int nodesInThisLevel, const deque< AVLNode<int,int>* >& nodesQueue, ostream& out) {
	deque<AVLNode<int,int>*>::const_iterator iter = nodesQueue.begin();
	for (int i = 0; i < nodesInThisLevel / 2; i++) {  
		out << ((i == 0) ? setw(startLen-1) : setw(nodeSpaceLen-2)) << "" << ((*iter++) ? "/" : " ");
		out << setw(2*branchLen+2) << "" << ((*iter++) ? "\\" : " ");
	}
	out << endl;
}

// Test function to Print the branches and node (eg, ___10___ )
void printNodes(int branchLen, int nodeSpaceLen, int startLen,
	int nodesInThisLevel, const deque<AVLNode<int,int>*>& nodesQueue, ostream& out) {
	deque<AVLNode<int,int>*>::const_iterator iter = nodesQueue.begin();
	for (int i = 0; i < nodesInThisLevel; i++, iter++) {
		out << ((i == 0) ? setw(startLen) : setw(nodeSpaceLen)) << ""
		<< ((*iter && (*iter)->llink) ? setfill('_') : setfill(' '));
		out << setw(branchLen+2) << ((*iter) ? intToString((*iter)->m_key) : "");
		out << ((*iter && (*iter)->rlink) ? setfill('_') : setfill(' '))
		<< setw(branchLen) << "" << setfill(' ');
	}
	out << endl;
}

// Test function to Print the leaves only (just for the bottom row)
void printLeaves (int indentSpace, int level, int nodesInThisLevel,
	const deque<AVLNode<int,int>*>& nodesQueue, ostream& out) {
	deque<AVLNode<int,int>*>::const_iterator iter = nodesQueue.begin();
	for (int i = 0; i < nodesInThisLevel; i++, iter++) {
		out << ((i == 0) ? setw(indentSpace+2) : setw(2*level+2))
		<< ((*iter) ? intToString((*iter)->m_key) : "");
	}
	out << endl;
}

// Formatting of a binary tree to the output stream
// @param level  Controls how wide you want the tree to sparse 
// (eg, level 1 has the minimum space between nodes, while level 2 has a larger space between nodes)
// @param indentSpace  Change this to add some indent space to the left (eg, indentSpace of 0
// means the lowest level of the left node will stick to the left margin)
void printInt(AVLNode<int,int> *root, int level, int indentSpace, ostream& out) {
	if(root==nullptr)
		return;
	int h = root->m_ht;
	int nodesInThisLevel = 1;

	//! eq of the length of branch for each node of each level
	int branchLen = 2*((int)pow(2.0,h)-1) - (3-level)*(int)pow(2.0,h-1);
	//! distance between left neighbor node's right arm and right
	//! neighbor node's left arm
	int nodeSpaceLen = 2 + (level+1)*(int)pow(2.0,h);
	//! starting space to the first node to print of each level
	//! (for the left most node of each level only)
	int startLen = branchLen + (3-level) + indentSpace;

	deque<AVLNode<int,int>*> nodesQueue;
	nodesQueue.push_back(root);
	for (int r = 1; r < h; r++) {
		printBranches(branchLen, nodeSpaceLen, startLen, nodesInThisLevel, nodesQueue, out);
		branchLen = branchLen/2 - 1;
		nodeSpaceLen = nodeSpaceLen/2 + 1;
		startLen = branchLen + (3-level) + indentSpace;
		printNodes(branchLen, nodeSpaceLen, startLen, nodesInThisLevel, nodesQueue, out);

		for (int i = 0; i < nodesInThisLevel; i++) {
			AVLNode<int,int> *currNode = nodesQueue.front();
			nodesQueue.pop_front();
			if (currNode) {
				nodesQueue.push_back(currNode->llink);
				nodesQueue.push_back(currNode->rlink);
			} else {
				nodesQueue.push_back(NULL);
				nodesQueue.push_back(NULL);
			}
		}
		nodesInThisLevel *= 2;
	}
	printBranches(branchLen, nodeSpaceLen, startLen, nodesInThisLevel, nodesQueue, out);
	printLeaves(indentSpace, level, nodesInThisLevel, nodesQueue, out);
}

template <> 
void Avl < int,int>::print(){
	printInt(root, 1, 0, cout);
}

//! Test function that displays a tree of ~30 random elements 
//! and does two inserts and deletes.
void demoTestCase(int num_inserts) {
	int seed = time(NULL);
	srand(seed);
    Avl< int, int > obj_map;
	AVLNode < int, int > *root_1;
	for(int i=0;i<num_inserts;i++)
	{
		int key=rand()%30;
	    root_1=obj_map.insert(key,key);
	}	
	cout<<"=========Tree=======";
	obj_map.print();
	int newkey = rand()%30+20;
	obj_map.insert(newkey,newkey);
	cout<<"Tree after inserting "<<newkey<<endl;
	obj_map.print();
	newkey = rand()%30+40;
	obj_map.insert(newkey,newkey);
	cout<<"Tree after inserting "<<newkey<<endl;
	obj_map.print();
	newkey = rand()%30;
	obj_map.remove(newkey);
	cout<<"Tree after deleting "<<newkey<<endl;
	obj_map.print();
	newkey = rand()%30;
	obj_map.remove(newkey);
	cout<<"Tree after deleting "<<newkey<<endl;
	obj_map.print();
}

// Test function that creates a tree of ~10000 random elements
// and does num_operations of look ups.
void runLookupTestCases(int num_operations){
	Avl< int, int > obj_map;
	AVLNode < int, int > *inserted_node;
	map<int, int> *intMap = new map<int,int>();
	if ( intMap == 0) 
	{
		throw OOM;
	}
	for(int i=0;i<num_operations;i++)
	{
		int key=rand()%10000;
		inserted_node=obj_map.insert(key,key);
		(*intMap)[key]=key;
	}
	//! We use this seed for both avl and red black trees
	//! so that same sequence of random numbers is generated.
	int seed = time(NULL);
	clock_t begin, end;
	srand(seed);
	begin = clock();
	for(int i=0;i<num_operations;i++)
	{
		int key=rand()%10000;
		inserted_node=obj_map.lookup(key);
	}
	end = clock();
	cout<<num_operations<<" random lookups in AVL tree takes ";
	cout<<double(end - begin)/CLOCKS_PER_SEC<<" seconds."<<endl ;
	srand(seed);
	map<int,int>::iterator it;
	begin = clock();
	for(int i=0;i<num_operations;i++)
	{
		int key=rand()%10000;
        it=intMap->find(key);
	}
	end = clock();
	cout<<num_operations<<" random lookups in Red Black tree takes ";
	cout<<double(end - begin)/CLOCKS_PER_SEC<<" seconds."<<endl ;
}

// Test function that creates a tree of ~10000 random elements
// and does num_operations of insertions with ~50% removals.
void runInsertionRemovalTestCases(int num_operations){
	//! We use this seed for both avl and red black trees
	//! so that same sequence of random numbers is generated.
	int seed = time(NULL);
	clock_t begin, end;
	Avl< int, int > obj_map;
	AVLNode < int, int > *inserted_node;
	srand(seed);
	begin = clock();
	for(int i=0;i<num_operations;i++)
	{
		int key=rand()%10000;
		if (!(key%2))
		{
			obj_map.remove(key);
		}
		else
			inserted_node=obj_map.insert(key,key);
	}
	end = clock();
	cout<<num_operations<<" random insertions with ~50% removals in AVL tree takes ";
	cout<<double(end - begin)/CLOCKS_PER_SEC<<" seconds."<<endl ;
	map<int, int> *intMap=new map<int,int>();
	if ( intMap == 0) 
	{
		throw OOM;
	}
	map<int,int>::iterator it;
	srand(seed);
	begin = clock();
	for(int i=0;i<num_operations;i++)
	{
		int key=rand()%10000;
		if ((key%2)) 
		{
			(*intMap)[key]=key;
		}
		else
		{
        	it=intMap->find(key);
			if(it!=intMap->end())
				intMap->erase(it);
		}
	}
	end = clock();
	cout<<num_operations<<" random insertions with ~50% removals in Red Black tree takes ";
	cout<<double(end - begin)/CLOCKS_PER_SEC<<" seconds."<<endl ;	
}

int main (){
	int n, number,v;
	char a;
	AVLNode < int,int > *root;
	AVLNode < int, int > *root_1;
	Avl < int,int> obj;
	Avl < int, int > obj_map;
	int seed = time(NULL);
	srand(seed);
	int cases=100000;
	demoTestCase(100);
	for (int i=0;i<3;i++) 
	{
		runInsertionRemovalTestCases(cases);
		runLookupTestCases(cases);
		cases*=10;
	}
	return 0;
}
