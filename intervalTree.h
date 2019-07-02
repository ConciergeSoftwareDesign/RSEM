#include <iostream> 
#include "Transcript.h"

using namespace std; 
  
struct ITInterval
{ 
    int low, high; 
}; 
  
struct ITNode
{ 
    ITInterval *interval; 
    int max; 
    ITNode *left, *right;
    Transcript *transcript;
}; 
  
ITNode* newNode(ITInterval i, Transcript *t) 
{ 
    ITNode *temp = new ITNode; 
    temp->interval = new ITInterval(i); 
    temp->max = i.high; 
    temp->transcript = t;
    temp->left = temp->right = NULL; 
}; 
  
ITNode* insert(ITNode* root, ITInterval i, Transcript *t) 
{ 
    // If the tree is empty, insert the new node as the root
    if (root == NULL) 
        return newNode(i, t); 
  
    // Get low value of interval at root 
    int lowValue = root->interval->low; 
  
    // If root's low value is smaller, then new interval goes to 
    // left subtree 
    if (i.low < lowValue) 
        root->left = insert(root->left, i, t); 
  
    // Else, new node goes to right subtree. 
    else
        root->right = insert(root->right, i, t); 
  
    // Update the max value of this ancestor if needed 
    if (root->max < i.high) 
        root->max = i.high; 
  
    return root; 
} 
  
// A utility function to check if given two intervals overlap 
bool intervalsOverlap(ITInterval i1, ITInterval i2) 
{ 
    if (i1.low <= i2.high && i2.low <= i1.high) 
        return true; 
    return false; 
} 
  
// The main function that searches a given interval i in a given 
// Interval Tree. 
ITNode* overlapSearch(ITNode *root, ITInterval i) 
{ 
    // Base Case, tree is empty 
    if (root == NULL) return NULL; 
  
    // If given interval overlaps with root 
    if (intervalsOverlap(*(root->interval), i)) 
        return root; 
  
    // If left child of root is present and max of left child is 
    // greater than or equal to given interval, then i may 
    // overlap with an interval is left subtree 
    if (root->left != NULL && root->left->max >= i.low) 
        return overlapSearch(root->left, i); 
  
    // Else interval can only overlap with right subtree 
    return overlapSearch(root->right, i); 
} 