/*
 * $ Author:Kaneez Fathima $
 * $ Date: 2014-03-15 $
 *
 * This file contains a templatised implementation of Strassen 
 * multiplication of two square matrices of size 2^n X 2^n.
 * Given two matrices a and b, where a11, a12, a21, a22 represent
 * 4 sub matrices of a, b11, b12, b21, b22 represent 4 submatrices of
 * b and c11, c12, c21, c22 represent 4 submatrices of resultant Matrix
 * c. The Strassen Matrix Multiplication algorithm is: 
 *
 *      q7 = (a12-a22)(b21+b22)
 *      q6 = (-a11+a21)(b11+b12)
 *      q5 = (a11+a12)b22
 *      q4= a22(-b11+b21)
 *      q3 = a11(b12-b22)
 *      q2 = (a21+a22)b11
 *      q1 = (a11+a22)(b11+b22)
 *      c11 = q1+q4-q5+q7
 *      c12 = q3+q5
 *      c21 = q2+q4
 *      c22 = q1+q3-q2+q6
*/
#include<iostream>
#include<cstdlib>
#include<ctime>
#include<stdexcept>
#include<iomanip>
#include<ogdf/basic/Array2D.h>

using namespace ogdf;
using namespace std;

#define OGDF_NUM_SUB_MATRICES 4

/**
* Matrix is represented as quad tree containing 4 sub matrices
* The actual data is in **d at the leaf nodes of the quad tree. Each internal
* node has 4 children representing sub matrices of size n/2 X n/2.
* m_basesize are set for the leave nodes indicating size of the matrix
*/
template<typename T> struct Matrix{
	T **d;
	int m_basesize;
	Matrix<T> **m; 
	//! Adds 2 matrices represented in Quadtree form. Used by Strassen Multiplication algorithm.
	static Matrix<T>* add(Matrix<T> *left, Matrix<T> *right, int nsize, int threshold);
	//! Subtracts 2 matrices represented in Quadtree form. Used by Strassen Multiplication algorithm.
	static Matrix<T>* sub(Matrix<T> *left, Matrix<T> *right, int nsize, int threshold);
	//! Multiplies 2 matrices represented in Quadtree form. Actual Strassen Multiplication core
	//! algorithm is in here.
	static Matrix<T>* multiply(Matrix<T> *left, Matrix<T> *right, int nsize, int threshold);
	static void allocateMatrix(Matrix<T> *res, int nsize);
	~Matrix();
};


template<typename T>
class StrassMatrix{
	//! Points to the matrix in Quadtree form.
	Matrix<T> *m_matrix;
	//! Size of the matrix.
	int m_size;
	//! Threshold Size below which matrix is not divided as submatrices
	//! and is represented as 
	int m_threshold;
	//! Converts a given matrix in Quad Tree form into a 2 dimensional array.
	void flattenMatrix(T **r, Matrix<T> *m_matrix, int startrow, int startcol, int nsize);
	//! Converts a given matrix in Quad Tree form into an Array2D.
	void flattenMatrix(Array2D<T> *r, Matrix<T> *m_matrix, int startrow, int startcol, int nsize);
	//! Creates a matrix Quadtree(Strass matrix) from a 2 dimensional matrix.
	Matrix<T>* createStrassMatrix(T **p, int startrow, int startcol, int nsize); 
	//! Creates a matrix Quadtree from Array2D instance.
	Matrix<T>* createStrassMatrix (Array2D<T> &p, int startrow, int startcol, int nsize);
	StrassMatrix(Matrix<T> *m, int nsize, int threshold);
	//! Made the constructor private to prevent instance creation without matrix.
	StrassMatrix(){}
	public:
	//! Creates a Strass Matrix from a 2 dimensional matrix, which gets represented
	//! as Quadtree internally.
	/**
	 * @param p is the pointer to 2D matrix.
	 * @param nsize is the number of rows, columns of square matrix.
	 * @param threshold is the maximum allowed size of
	 *  2D matrices stored at leaf nodes. default is 10.
	 * if threshold is >= current matrix size it just is same as
	 * using normal matrix operation.
	 */
	StrassMatrix(T **p, int nsize, int threshold=10);
	//! Creates a Strass Matrix from a Array2D instance, which gets represented
	//! as Quadtree internally.
	/**
	 * @param p is the reference to Array2d Object.
	 * @param threshold is the maximum allowed size of
	 *  2D matrices stored at leaf nodes. default is 10.
	 */
	StrassMatrix(Array2D<T> &p, int threshold=10);
	~StrassMatrix();
	//! Converts the data in Strass Matrix form into a 2 dimensional array.
	// @param r is pointer to 2D array.
	void strassMatrixtoMatrix(T** &r);
	//! Converts the data in Strass Matrix form into an Array2D object.
	// @param r is pointer to Array2D object.
	void strassMatrixtoArray2D(Array2D<T>* &r);
	int getThreshold() { return m_threshold; }
	//! Adds 2 Strass Matrices and returns reference to resultant Strass Matrix.
	template<typename Y>
		friend StrassMatrix<Y>& operator+(StrassMatrix<Y> &left, StrassMatrix<Y> &right);
	//! Subtracts 2 Strass Matrices and returns reference to resultant Strass Matrix.
	template<typename Y>
		friend StrassMatrix<Y>& operator-(StrassMatrix<Y> &left, StrassMatrix<Y> &right);
	//! Multiplies 2 Strass Matrices and returns reference to resultant Strass Matrix.
	template<typename Y>
		friend StrassMatrix<Y>& operator*(StrassMatrix<Y> &left, StrassMatrix<Y> &right);

};

template<typename T>
StrassMatrix<T>::StrassMatrix(Matrix<T> *m, int nsize, int threshold){
	m_size=nsize;
	m_matrix=m;
	m_threshold=threshold;
}

template<typename T>
StrassMatrix<T>::StrassMatrix(T **p, int nsize, int threshold){
	m_size=nsize;
	m_threshold=threshold;
	m_matrix=createStrassMatrix(p,0,0,m_size);
} 

template<typename T>
StrassMatrix<T>::StrassMatrix(Array2D<T> &p, int threshold){
	if ((p.high1()-p.low1())!=(p.high2()-p.low2())) {
		throw invalid_argument( "Array2D must be a square matrix." );
	}
	m_size=p.high1()-p.low1()+1;
	m_threshold=threshold;
	m_matrix=createStrassMatrix(p,p.low1(),p.low2(),m_size);
} 

template<typename T> StrassMatrix<T>::~StrassMatrix(){
	if (m_matrix!=NULL) {
		delete m_matrix;
	}
}

template<typename T>
void StrassMatrix<T>::strassMatrixtoMatrix(T** &r) {
	r = new T*[m_size];
	if (r==0)
		OGDF_THROW(InsufficientMemoryException);
	for (int i=0; i<m_size; i++) { 
		r[i]= new T[m_size];
		if (r[i]==0)
			OGDF_THROW(InsufficientMemoryException);
	}
	flattenMatrix(r,m_matrix,0,0,m_size);
}

template<typename T>
void StrassMatrix<T>::flattenMatrix(Array2D<T> *r, Matrix<T> *m_matrix, int startrow, int startcol, int nsize ) {
	if (nsize>m_threshold) {
		flattenMatrix(r,m_matrix->m[0],startrow,startcol,nsize/2);
		flattenMatrix(r,m_matrix->m[1],startrow,startcol+nsize/2,nsize/2);
		flattenMatrix(r,m_matrix->m[2],startrow+nsize/2,startcol,nsize/2);
		flattenMatrix(r,m_matrix->m[3],startrow+nsize/2,startcol+nsize/2,nsize/2);    
	} else {
		for(int i=0;i<nsize;i++) {
			for(int j=0;j<nsize;j++) {
				(*r)(startrow+i,startcol+j)=m_matrix->d[i][j];
			} 
		}
	}
}
		
template<typename T>
void StrassMatrix<T>::strassMatrixtoArray2D(Array2D<T>* &r) {
	r = new Array2D<T>(0,m_size-1,0,m_size-1);
	if (r==0)
		OGDF_THROW(InsufficientMemoryException);
	flattenMatrix(r,m_matrix,0,0,m_size);
}

template<typename T>
Matrix<T>* StrassMatrix<T>::createStrassMatrix (T **p, int startrow, int startcol, int nsize) {
	Matrix<T> *m_matrix=new Matrix<T>;
	if (m_matrix==0)
		OGDF_THROW(InsufficientMemoryException);
	if (nsize>m_threshold) {
		m_matrix->m=new Matrix<T>*[OGDF_NUM_SUB_MATRICES];
		if (m_matrix->m==0)
			OGDF_THROW(InsufficientMemoryException);
		m_matrix->m[0]=createStrassMatrix(p,startrow,startcol,nsize/2);
		m_matrix->m[1]=createStrassMatrix(p,startrow,startcol+nsize/2,nsize/2);
		m_matrix->m[2]=createStrassMatrix(p,startrow+nsize/2,startcol,nsize/2);
		m_matrix->m[3]=createStrassMatrix(p,startrow+nsize/2,startcol+nsize/2,nsize/2);
	} else {
		Matrix<T>::allocateMatrix(m_matrix, nsize);
		for(int i=0;i<nsize;i++) {
			for(int j=0;j<nsize;j++) {
				m_matrix->d[i][j]=p[startrow+i][startcol+j];
			} 
		}
	}
	return m_matrix;
}

template<typename T>
Matrix<T>* StrassMatrix<T>::createStrassMatrix (Array2D<T> &p, int startrow, int startcol, int nsize) {
	Matrix<T> *m_matrix=new Matrix<T>;
	if (m_matrix==0)
		OGDF_THROW(InsufficientMemoryException);
	if (nsize>m_threshold) {
		m_matrix->m=new Matrix<T>*[OGDF_NUM_SUB_MATRICES];
		if (m_matrix->m==0)
			OGDF_THROW(InsufficientMemoryException);
		m_matrix->m[0]=createStrassMatrix(p,startrow,startcol,nsize/2);
		m_matrix->m[1]=createStrassMatrix(p,startrow,startcol+nsize/2,nsize/2);
		m_matrix->m[2]=createStrassMatrix(p,startrow+nsize/2,startcol,nsize/2);
		m_matrix->m[3]=createStrassMatrix(p,startrow+nsize/2,startcol+nsize/2,nsize/2);
	} else {
		Matrix<T>::allocateMatrix(m_matrix, nsize);
		for(int i=0;i<nsize;i++) {
			for(int j=0;j<nsize;j++) {
				m_matrix->d[i][j]=p(startrow+i,startcol+j);
			} 
		}
	}
	return m_matrix;
}

template<typename T>
void StrassMatrix<T>::flattenMatrix(T **r, Matrix<T> *m_matrix, int startrow, int startcol, int nsize ) {
	if (nsize>m_threshold) {
		flattenMatrix(r,m_matrix->m[0],startrow,startcol,nsize/2);
		flattenMatrix(r,m_matrix->m[1],startrow,startcol+nsize/2,nsize/2);
		flattenMatrix(r,m_matrix->m[2],startrow+nsize/2,startcol,nsize/2);
		flattenMatrix(r,m_matrix->m[3],startrow+nsize/2,startcol+nsize/2,nsize/2);    
	} else {
		for(int i=0;i<nsize;i++) {
			for(int j=0;j<nsize;j++) {
				r[startrow+i][startcol+j]=m_matrix->d[i][j];
			} 
		}
	}
}

template<typename T>
StrassMatrix<T>& operator+(StrassMatrix<T> &left, StrassMatrix<T> &right) {
	if ((left.m_threshold!=right.m_threshold)||(left.m_size!=right.m_size)) {
		throw invalid_argument( "Cannot Add Unequally sized of matrices." );
	}
	Matrix<T> *m3=Matrix<T>::add(left.m_matrix,right.m_matrix,left.m_size,left.m_threshold);
	StrassMatrix<T>* result = new StrassMatrix<T>(m3,left.m_size,left.m_threshold);
	if (result==0)
		OGDF_THROW(InsufficientMemoryException);
	return *result;
}

template<typename T>
StrassMatrix<T>& operator-(StrassMatrix<T> &left,StrassMatrix<T> &right) {
	if ((left.m_threshold!=right.m_threshold)||(left.m_size!=right.m_size)) {
		throw invalid_argument( "Cannot Subtract Unequally sized of matrices." );
	}
	Matrix<T> *m3=Matrix<T>::sub(left.m_matrix,right.m_matrix,left.m_size,left.m_threshold);
	StrassMatrix<T>* result = new StrassMatrix<T>(m3,left.m_size,left.m_threshold);
	if (result==0)
		OGDF_THROW(InsufficientMemoryException);
	return *result;
}

template<typename T>
StrassMatrix<T>& operator*(StrassMatrix<T> &left,StrassMatrix<T> &right) {
	if ((left.m_threshold!=right.m_threshold)||(left.m_size!=right.m_size)) {
		throw invalid_argument( "Cannot Multiply Unequally sized of matrices." );
	}
	Matrix<T> *m3=Matrix<T>::multiply(left.m_matrix,right.m_matrix,left.m_size,left.m_threshold);
	StrassMatrix<T>* result = new StrassMatrix<T>(m3,left.m_size,left.m_threshold);
	if (result==0)
		OGDF_THROW(InsufficientMemoryException);
	return *result;
}

template<typename T>
Matrix<T>* Matrix<T>::add(Matrix<T> *left, Matrix<T> *right, int nsize, int threshold) {
	Matrix<T> *res=new Matrix<T>;
	if (res==0)
		OGDF_THROW(InsufficientMemoryException);
	if (nsize > threshold) {
		res->m=new Matrix<T>*[OGDF_NUM_SUB_MATRICES];
		if (res->m==0)
			OGDF_THROW(InsufficientMemoryException);
		res->m[0]=add(left->m[0], right->m[0], nsize/2, threshold);
		res->m[1]=add(left->m[1], right->m[1], nsize/2, threshold);
		res->m[2]=add(left->m[2], right->m[2], nsize/2, threshold);
		res->m[3]=add(left->m[3], right->m[3], nsize/2, threshold);
	} else {
		allocateMatrix(res, nsize);
		for (int i = 0; i < nsize; i++) {
			for (int j = 0; j < nsize; j++) {
				res->d[i][j] = left->d[i][j] + right->d[i][j];
			}
		}
	}
	return res;
}

template<typename T>
Matrix<T>* Matrix<T>::sub(Matrix<T> *left, Matrix<T> *right, int nsize, int threshold) {
	Matrix<T> *res=new Matrix<T>;
	if (res==0)
		OGDF_THROW(InsufficientMemoryException);
	if (nsize > threshold) {
		res->m=new Matrix<T>*[OGDF_NUM_SUB_MATRICES];
		if (res->m==0)
			OGDF_THROW(InsufficientMemoryException);
		res->m[0]=sub(left->m[0], right->m[0], nsize/2, threshold);
		res->m[1]=sub(left->m[1], right->m[1], nsize/2, threshold);
		res->m[2]=sub(left->m[2], right->m[2], nsize/2, threshold);
		res->m[3]=sub(left->m[3], right->m[3], nsize/2, threshold);
	} else {
		allocateMatrix(res, nsize);
		for (int i = 0; i < nsize; i++) {
			for (int j = 0; j < nsize; j++) {
				res->d[i][j] = left->d[i][j] - right->d[i][j];
			}
		}
	}
	return res;
}

template<typename T>
Matrix<T>* Matrix<T>::multiply(Matrix<T> *a, Matrix<T> *b, int nsize, int threshold) {
	Matrix<T> *tmp=new Matrix<T>;
	if (tmp==0)
		OGDF_THROW(InsufficientMemoryException);
	Matrix<T> *c=new Matrix<T>;
	if (c==0)
		OGDF_THROW(InsufficientMemoryException);
	if (nsize > threshold) {
		tmp->m=new Matrix<T>*[OGDF_NUM_SUB_MATRICES];
		if (tmp->m==0)
			OGDF_THROW(InsufficientMemoryException);
		c->m=new Matrix<T>*[OGDF_NUM_SUB_MATRICES];
		if (c->m==0)
			OGDF_THROW(InsufficientMemoryException);
		tmp->m[0]=Matrix<T>::sub(a->m[1],a->m[3],nsize/2, threshold);
		tmp->m[1]=Matrix<T>::add(b->m[2],b->m[3],nsize/2, threshold);
		c->m[0]=Matrix<T>::multiply(tmp->m[0],tmp->m[1],nsize/2, threshold);
		tmp->m[0]=Matrix<T>::sub(a->m[2],a->m[0],nsize/2, threshold);
		tmp->m[1]=Matrix<T>::add(b->m[0],b->m[1],nsize/2, threshold);
		c->m[3]=Matrix<T>::multiply(tmp->m[0],tmp->m[1],nsize/2, threshold);
		tmp->m[0]=Matrix<T>::add(a->m[0],a->m[1],nsize/2, threshold);
		c->m[1]=Matrix<T>::multiply(tmp->m[0],b->m[3],nsize/2, threshold);
		c->m[0]=Matrix<T>::sub(c->m[0],c->m[1],nsize/2, threshold);
		tmp->m[0]=Matrix<T>::sub(b->m[2],b->m[0],nsize/2, threshold);
		c->m[2]=Matrix<T>::multiply(a->m[3],tmp->m[0],nsize/2, threshold);
		c->m[0]=Matrix<T>::add(c->m[2],c->m[0],nsize/2, threshold);
		tmp->m[0]=Matrix<T>::sub(b->m[1],b->m[3],nsize/2, threshold);
		tmp->m[1]=Matrix<T>::multiply(a->m[0],tmp->m[0],nsize/2, threshold);
		c->m[1]=Matrix<T>::add(tmp->m[1],c->m[1],nsize/2, threshold);
		c->m[3]=Matrix<T>::add(tmp->m[1],c->m[3],nsize/2, threshold);
		tmp->m[0]=Matrix<T>::add(a->m[2],a->m[3],nsize/2, threshold);
		tmp->m[1]=Matrix<T>::multiply(tmp->m[0],b->m[0],nsize/2, threshold);
		c->m[2]=Matrix<T>::add(tmp->m[1],c->m[2],nsize/2, threshold);
		c->m[3]=Matrix<T>::sub(c->m[3],tmp->m[1],nsize/2, threshold);
		tmp->m[0]=Matrix<T>::add(a->m[0],a->m[3],nsize/2, threshold);
		tmp->m[1]=Matrix<T>::add(b->m[0],b->m[3],nsize/2, threshold);
		tmp->m[2]=Matrix<T>::multiply(tmp->m[0],tmp->m[1],nsize/2, threshold);
		c->m[0]=Matrix<T>::add(tmp->m[2],c->m[0],nsize/2, threshold);
		c->m[3]=Matrix<T>::add(tmp->m[2],c->m[3],nsize/2, threshold);
		delete tmp->m[0];
		delete tmp->m[1];
		delete tmp->m[2];
	} else {
		allocateMatrix(c, nsize);
		T sum;
		int i, j, k;
		for (i = 0; i < nsize; i++) {
			for (j = 0; j < nsize; j++) {
				for (sum = 0, k = 0; k < nsize; k++)
					sum += a->d[i][k] * b->d[k][j];
					c->d[i][j] = sum;
			}
		}
	}
	return c;
}

template<typename T>
void Matrix<T>::allocateMatrix(Matrix<T> *res, int nsize){
	//res->d = unique_ptr<unique_ptr<T[]>[]>(new unique_ptr<T[]>[nsize]);
	res->d = new T*[nsize];
	if (res->d==0)
		OGDF_THROW(InsufficientMemoryException);
	for(int i=0; i<nsize; i++){ 
		//res->d[i] = unique_ptr<T[]>(new T[ncol]);
		res->d[i] = new T[nsize];
		if (res->d[i]==0)
			OGDF_THROW(InsufficientMemoryException);
	}
	res->m_basesize=nsize;
	res->m = NULL;
}

template<typename T> Matrix<T>::~Matrix() {
	if (m==NULL) {
		if (d!=NULL) {
			for (int i=0;i<m_basesize;i++)
				if(d[i]!=NULL)
					delete[] d[i];
			delete[] d;
		}
	} else {
		for (int i=0;i<OGDF_NUM_SUB_MATRICES;i++)
			if(m[i]!=NULL)
				delete m[i]; 
		delete[] m;
	}
}

//! Test function to create a matrix of random values.
template<typename T>
void createMatrix (T** &p, int nsize) {
	p = new T*[nsize];
	if (p==0)
		OGDF_THROW(InsufficientMemoryException);
	for(int i=0; i<nsize; i++){ 
		p[i]= new T[nsize];
		if (p[i]==0)
			OGDF_THROW(InsufficientMemoryException);
	}
	for(int i=0;i<nsize;i++) {
		for(int j=0;j<nsize;j++) {
			p[i][j]=rand()%100+(T)((rand()%100)/100);
		}
	}
}

//! Test function to print a matrix.
template<typename T> void print(T **p, int nsize) {
	for(int i=0;i<nsize;i++) {
		for(int j=0;j<nsize;j++) {
			cout<<p[i][j]<<"\t";
		}
		cout<<endl;
	}
}

//! Test function to print an Array2D instance.
template<typename T> void print(Array2D<T> &array2d) {
	int nrow = array2d.high1()-array2d.low1()+1;
	int ncol = array2d.high2()-array2d.low2()+1;
	for(int i=0;i<nrow;i++) {
		for(int j=0;j<ncol;j++) {
			cout<< std::fixed << array2d(i,j)<<"\t ";
		}
		cout<<endl;
	}
}

//! Test function to fill an Array2D with random values.
void fillArray2D (Array2D<double> &array2d) {
	static std::mt19937 rng(std::time(nullptr));
	std::normal_distribution<double> normal(0,100) ;
	int nrow = array2d.high1()-array2d.low1()+1;
	int ncol = array2d.high2()-array2d.low2()+1;
	for(int i=0;i<nrow;i++) {
		for(int j=0;j<ncol;j++) {
			array2d(i+array2d.low1(),j+array2d.low2())=normal(rng);
		}
	}
}

//! Test function to create a matrix from Array2D.
template<typename T>
void fillArray2DinMatrix (Array2D<T> &array2d, T **&r) {
	int nrow = array2d.high1()-array2d.low1()+1;
	int ncol = array2d.high2()-array2d.low2()+1;
	if (r==0)
		OGDF_THROW(InsufficientMemoryException);
	r = new T*[nrow];
	for(int i=0; i<nrow; i++) { 
		r[i]= new T[ncol];
		if (r[i]==0)
			OGDF_THROW(InsufficientMemoryException);
	}
	for(int i=0;i<nrow;i++) {
		for(int j=0;j<ncol;j++) {
			r[i][j]=array2d(i+array2d.low1(),j+array2d.low2());
		}
	}
}

//! Test function to check if 2 matrices are equal.
template<typename T>
bool matricesAreEqual(T **matrix1, T **matrix2, int nrow, int ncol) {
	for(int i=0;i<nrow;i++) {
		for(int j=0;j<ncol;j++) {
			if((matrix1[i][j]-matrix2[i][j])>0.0001) {
				cout<<matrix1[i][j]<<"\t"<<matrix2[i][j]<<endl;
				return false;
			}
		}
	}
	return true;
}
// Enum used by tests.
enum MatrixOperations{
	kAddition,
	kSubtraction,
	kMultiplication,
};

//! Generates 2 random matrices of given size and Tests Matrix operations on Strass
//! Matrices by comparing results with normal matrix operations.
/**
* @param matrixSize is size of random matrices to generated.
* @param threshold is max size at leaf node of Strass Matrix Quadtree.
* @param operation is the type of Matrix operation to test, kAddition, kSubtraction or
*  kMultiplication. 
*/
template<typename T>
bool TestMatrixOperation(int matrixSize, int threshold, MatrixOperations operation) {
	int nrow = matrixSize;
	int ncol = matrixSize; 
	Array2D<T> firstArray2d(0,nrow-1,0,ncol-1);
	fillArray2D(firstArray2d);
	Array2D<T> secondArray2d(0,nrow-1,0,ncol-1);
	fillArray2D(secondArray2d);
	StrassMatrix<T> *smArray= new StrassMatrix<T>(firstArray2d,threshold);
	StrassMatrix<T> *smArray2= new StrassMatrix<T>(secondArray2d,threshold);
	T **arrayValues;
	T **arrayValues2;
	// With threshold param set to matrix Size, this is effectively just creating
	// normal matrix and operations are usual addition, subtraction and O(N^3)
	// multiplication.
	StrassMatrix<T> *smArray4 = new StrassMatrix<T>(firstArray2d,matrixSize);
	StrassMatrix<T> *smArray5 = new StrassMatrix<T>(secondArray2d,matrixSize);
	clock_t begin, end, begin2, end2;
	if (operation == kAddition) {
		StrassMatrix<T> smArray3=*smArray + *smArray2;
		StrassMatrix<T> smArray6=*smArray4 + *smArray5;
		smArray3.strassMatrixtoMatrix(arrayValues);
		smArray6.strassMatrixtoMatrix(arrayValues2);
	}
	if (operation == kSubtraction) {
		StrassMatrix<T>	smArray3=(*smArray - *smArray2);
		StrassMatrix<T>	smArray6=(*smArray4 - *smArray5);
		smArray3.strassMatrixtoMatrix(arrayValues);
		smArray6.strassMatrixtoMatrix(arrayValues2);
	}
	if (operation == kMultiplication) {
		begin = clock();
		StrassMatrix<T>	smArray3=(*smArray * *smArray2);
		end = clock();
		begin2 = clock();
		StrassMatrix<T>	smArray6=(*smArray4 * *smArray5);
		end2 = clock();
		Array2D<T> *array2d;
		smArray3.strassMatrixtoArray2D(array2d);
		fillArray2DinMatrix(*array2d,arrayValues);
		smArray6.strassMatrixtoMatrix(arrayValues2);
		// Just for comparision the time taken by each type of multiplication 
		// is printed out.
		cout<<"Strass Multiplication for matrix size "<<matrixSize<<" takes ";
		cout<<double(end - begin)/CLOCKS_PER_SEC<<" seconds."<<endl ;
		cout<<"Normal Multiplication for matrix size "<<matrixSize<<" takes ";
		cout<<double(end2 - begin2)/CLOCKS_PER_SEC<<" seconds."<<endl ;
	}
	// arrayValues has strassen matrix results and arrayValues2 has normal matrix
	// operation results.
	if (!matricesAreEqual(arrayValues, arrayValues2, nrow, ncol)) {
		return false;
	}
	for(int j=0;j<matrixSize;j++) {
		delete[] arrayValues[j];
		delete[] arrayValues2[j];
	}
	delete[] arrayValues;
	delete[] arrayValues2;
	delete smArray;
	delete smArray2;
	delete smArray4;
	delete smArray5;
	return true;
}

// Runs testcases by generating random matrices of sizes 2^0, 2^1..2^10.
// Threshold size at leaf node are varied from 1 to 150.
template<typename T>
void RunTestCases() {
	int threshold[] = {1,1,2,2,4,4,5,10,20,40,70,150};
	for(int i=1,j=0; i<=1024 && j<=10; i<<=1,j++) {
		if (!TestMatrixOperation<T>(i, threshold[j], kAddition)) {
			cout<<"Strass Matrix Addition failed."<<i<<endl;
			return;
		}
		if (!TestMatrixOperation<T>(i, threshold[j], kSubtraction)) {
			cout<<"Strass Matrix Subtraction failed."<<i<<endl;
			return;
		}
		if (!TestMatrixOperation<T>(i, threshold[j], kMultiplication)) {
			cout<<"Strass Matrix Multiplication failed."<<i<<endl;
			return;
		}
	}
	cout<<"Strass Matrix Addition testcases passed."<<endl;
	cout<<"Strass Matrix Subtraction testcases passed."<<endl;
	cout<<"Strass Matrix Multiplication testcases passed."<<endl;
}

//! Multiplies 2 Array2D objects and returns a resultant Array2D
//! object. This can be used to overload operator* for Array2D.
//! Square matrices of size 2^n X 2^n are expected as input.
template<typename T>
Array2D<T>& multiplyArray2D(Array2D<T> &a1, Array2D<T> &a2) {
	StrassMatrix<T> *smArray= new StrassMatrix<T>(a1,1);
	StrassMatrix<T> *smArray2= new StrassMatrix<T>(a2,1);
	StrassMatrix<T> smArray3= *smArray * *smArray2;
	Array2D<T> *a3;
	smArray3.strassMatrixtoArray2D(a3);
	return *a3;
}

int main(int argc, char *argv[]) {
	int nsize,ncol;
	nsize=4;
	ncol=4;
	cout<<"Multiplying two 4X4 matrices"<<endl;
	Array2D<double> firstArray2d(0,nsize-1,0,ncol-1);
	fillArray2D(firstArray2d);
	cout<<"first matrix"<<endl<<endl;
	print<double>(firstArray2d);
	cout<<endl<<endl;
	Array2D<double> secondArray2d(0,nsize-1,0,ncol-1);
	fillArray2D(secondArray2d);
	cout<<"second matrix"<<endl<<endl;
	print<double>(secondArray2d);
	cout<<endl<<endl;
	// Demo of multiplying two 4 X 4 matrices
	Array2D<double> a3 = multiplyArray2D(firstArray2d, secondArray2d);
	cout<<"third matrix"<<endl<<endl;
	print<double>(a3);
	RunTestCases<double>();
}
