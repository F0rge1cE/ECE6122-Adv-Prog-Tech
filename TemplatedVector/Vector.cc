// Implementation of the templated Vector class
// ECE4893/8893 lab 3
// Xueyang Xu

#include <iostream> // debugging
#include "Vector.h"

// Your implementation here
// Fill in all the necessary functions below
using namespace std;

// Default constructor
template <typename T>
Vector<T>::Vector()
{
	elements = NULL;
	count = 0;
	reserved = 0;
}

// Copy constructor
template <typename T>
Vector<T>::Vector(const Vector& rhs)
{
	count = rhs.count;
	reserved = rhs.reserved;
	elements = ( T* ) malloc(sizeof(T) * rhs.count);
	for (size_t i = 0; i < count; ++i)
	{
		new (&elements[i]) T(rhs[i]); // In-place new
	}
}

// Assignment operator
template <typename T>
Vector<T>& Vector<T>::operator=(const Vector& rhs)
{
	// Clear the left hand vector
	for (size_t i = 0; i < count; ++i){
		elements[i].~T();
	}
	free(elements);

	// Copy RHS vector to the left
	count = rhs.count;
	reserved = rhs.reserved;
	elements = ( T* )malloc(sizeof(T) * rhs.count);
	for (size_t i = 0; i < count ; ++i){
		new (&elements[i]) T(rhs[i]); // In-place new
	}
	return *this;
}

#ifdef GRAD_STUDENT
// Other constructors
template <typename T>
Vector<T>::Vector(size_t nReserved)
{ // Initialize with reserved memory
	count = 0;
	reserved = nReserved;
	elements = ( T* )malloc(sizeof(T) * nReserved);
}

template <typename T>
Vector<T>::Vector(size_t n, const T& t)
{ // Initialize with "n" copies of "t"
	count = n;
	reserved = 0;
	elements = ( T* )malloc(sizeof(T) * n);
	for (size_t i = 0; i < n; ++i){
		new (&elements[i])  T(t);
	}
}

template <typename T>
void Vector<T>::Reserve(size_t n)
{ // Reserve extra memory
	T* new_elements = ( T* )malloc(sizeof(T) * n);
	for (size_t i = 0; i < count; ++i){
		new (&new_elements[i]) T(elements[i]);	
		// Destruct previous vector
		elements[i].~T();
	}
	free(elements);
	//count = count; // count does not change
	reserved = n;
	elements = new_elements;
}

#endif

// Destructor
template <typename T>
Vector<T>::~Vector()
{
	for (size_t i = 0; i < count; ++i){
		elements[i].~T();
	}
	free(elements);

	elements = NULL;
	count = 0;
	reserved = 0;
}

// Add and access front and back
template <typename T>
void Vector<T>::Push_Back(const T& rhs)
{
	if ((count+1) > reserved){
		Reserve(count+1);
	}
	new (&elements[count]) T(rhs);
	count ++;

}

template <typename T>
void Vector<T>::Push_Front(const T& rhs)
{
	if (count < reserved){
		for (size_t i=count; i > 0; --i){
			//new (&elements[i]) T(elements[i-1]);  
			// Use in-place new?? 
			elements[i] = elements[i-1];
		}
		new (&elements[0]) T(rhs);	
		count++;
		//reserved does not change
	}
	else{
		// Similar to Reserve(), but reduces copy_count.
		// Do not use Reserve(), it will do unnecessary copies
		T* new_elements = ( T* )malloc(sizeof(T) * (count+1));

		for (size_t i = count; i > 0 ; --i){
			new (&new_elements[i]) T(elements[i-1]);	
			elements[i-1].~T();
		}
		free(elements);
		elements = new_elements;
		new (&elements[0]) T(rhs);
		//new (&new_elements[0]) T(rhs);
		//elements = new_elements;
		count++; 
		reserved = count;
	}
}

template <typename T>
void Vector<T>::Pop_Back()
{ // Remove last element
	elements[(count-1)].~T();
	count --;
	// count--;
	// elements[count].~T();
}

template <typename T>
void Vector<T>::Pop_Front()
{ // Remove first element
	
	for (size_t i = 0; i < count-1; ++i){
		elements[i].~T();
		new (&elements[i]) T(elements[i+1]); // Use in-place new??
	}
	elements[count-1].~T();
	count --;
}

// Element Access
template <typename T>
T& Vector<T>::Front() const
{
	return elements[0];
}

// Element Access
template <typename T>
T& Vector<T>::Back() const
{
	return elements[count-1];
}

template <typename T>
const T& Vector<T>::operator[](size_t i) const
{ // const element access
	return (*this).elements[i];
}

template <typename T>
T& Vector<T>::operator[](size_t i)
{//nonconst element access
	return (*this).elements[i];
}

template <typename T>
size_t Vector<T>::Size() const
{
	return count;
}

template <typename T>
bool Vector<T>::Empty() const
{
	return (count == 0);
}

// Implement clear
template <typename T>
void Vector<T>::Clear()
{
	for(size_t i = 0; i < count; ++i ){
		elements[i].~T();
	}
	count = 0;
}

// Iterator access functions
template <typename T>
VectorIterator<T> Vector<T>::Begin() const
{
  return VectorIterator<T>(elements);
}

template <typename T>
VectorIterator<T> Vector<T>::End() const
{
	return VectorIterator<T>(elements+count);
}

#ifdef GRAD_STUDENT
// Erase and insert
template <typename T>
void Vector<T>::Erase(const VectorIterator<T>& it)
{
	// Find it.current's index in the vector
	size_t ind;
	for (size_t i = 0; i < count; ++i)
	{
		if (it.current == &elements[i]){
			ind = i;
			break;
		}
	}
	// Regenerate elements after "ind"
	// Similar to pop_front()
	for (size_t i = ind; i < count-1; ++i){
		elements[i].~T();
		new (&elements[i]) T(elements[i+1]); 
	}
	
	elements[count-1].~T();
	count --;
}

template <typename T>
void Vector<T>::Insert(const T& rhs, const VectorIterator<T>& it)
{
	// Find it.current's index in the vector
	size_t ind;
	for (size_t i = 0; i < count; ++i)
	{
		if (it.current == &elements[i]){
			ind = i;
			break;
		}
	}

	// Regenerate elements after "ind"
	// Similar to push_front()
	if (count < reserved){
		for (size_t i = count; i > ind; --i){
			elements[i] = elements[i-1];
		}
		new (&elements[ind]) T(rhs);	
		count++;
		}
	else{
		T* new_elements = ( T* )malloc(sizeof(T) * (count+1));

		for(size_t i = 0; i < ind; ++i){
			new (&new_elements[i]) T(elements[i]);	
		 	elements[i].~T();
		}
		// for(size_t i = (ind+1); i < count+1; ++i){
		// 	new (&new_elements[i]) T(elements[i-1]);
		// 	elements[i-1].~T();
		// }
		for (size_t i = count; i > ind; --i ){
			new (&new_elements[i]) T(elements[i-1]);
			elements[i-1].~T();
		}
		free(elements);
		new (&new_elements[ind]) T(rhs);
		elements = new_elements;
		count++; 
		reserved = count;
	}

}
#endif

// Implement the iterators

// Constructors
template <typename T>
VectorIterator<T>::VectorIterator()
{
	current = NULL;
}

template <typename T>
VectorIterator<T>::VectorIterator(T* c)
{
	current = c;
}

// Copy constructor
template <typename T>
VectorIterator<T>::VectorIterator(const VectorIterator<T>& rhs)
{
	current = rhs.current;
}

// Iterator defeferencing operator
template <typename T>
T& VectorIterator<T>::operator*() const
{
	return *current;
}

// Prefix increment
template <typename T>
VectorIterator<T>  VectorIterator<T>::operator++()
{
	current ++;
	return *this;
}

// Postfix increment
template <typename T>
VectorIterator<T> VectorIterator<T>::operator++(int)
{	
	VectorIterator<T> before_incre(*this);
	current ++;
	return before_incre;
}

// Comparison operators
template <typename T>
bool VectorIterator<T>::operator !=(const VectorIterator<T>& rhs) const
{
	return (current != rhs.current);
}

template <typename T>
bool VectorIterator<T>::operator ==(const VectorIterator<T>& rhs) const
{
	return (current == rhs.current);
}




