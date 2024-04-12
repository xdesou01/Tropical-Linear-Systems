#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <list>
#include <queue>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <set>
#include <climits>
#include <assert.h>

using namespace std;

vector<vector<int>> retrieveMatrix() {
    //vector<vector<int>> x(4, vector<int>(5, 1)); //change for input
    // vector<vector<int>> x = {{1, 0}}; //works
    // vector<vector<int>> x = {{0, 1}}; //works
    // vector<vector<int>> x = {{1, 2, 3}, {3, 2, 1}}; //works, solution found
    vector<vector<int>> x = {{0, 1, 2, 3}, 
                             {0, 0, 1, 3}, 
                             {1, 0, 0, 2}, 
                             {0, 1, 0, 2}}; //works, solution found
    // vector<vector<int>> x = {{1, 2, 3, 4}};
    //vector<vector<int>> x = {{0, 1, 2, 3}, {0, 0, 1, 3}, {1, 0, 0, 2}, {0, 1, 0, 2}};
    // vector<vector<int>> x = {{1, 1, 0, 3}, {1, 2, 1, 3}};
    //vector<vector<int>> x = {{1, 2}, {3, 2}}; //works, no solution detected
    return x;
}


void insertIfNotInJ(int j, unordered_set<int> &J, unordered_set<int> &curMinColsNotInJ) {
    if (J.find(j) == J.end()) curMinColsNotInJ.insert(j);
}


class Matrix {
    private:
    vector<vector<int>> mat;

    public:
    Matrix(vector<vector<int>> &v) {this->mat = v;}
    int get(int row, int col) {
        return mat[row - 1][col - 1];
    }
    void set(int row, int col, int val) {
        mat[row - 1][col - 1] = val;
    }
    int rows() {return mat.size();}
    int cols() {
        if (mat.size() == 0) return 0; //empty case
        else                 return mat[0].size();
    }
    void swapCols(int c1, int c2) {
        for (int i = 0; i < mat.size(); i++) {
            std::swap(mat[i][c1 - 1], mat[i][c2 - 1]);
        }
    }
    void swapRows(int r1, int r2) {
        std::swap(mat[r1 - 1], mat[r2 - 1]);
    }
    void addRow(vector<int> &v) {
        mat.insert(mat.begin(), v);
    }
};

class SolnVector {
    private:
    vector<int> vec;

    public:
    SolnVector(int n) {
        vector<int> v(n, 0);
        vec = v;
    }
    void add(int index, int val) {vec[index - 1] += val;}
    int get(int index) {return vec[index - 1];}
    void printSolution() {
        cout << "SOLUTION:" << endl;
        for (int j = 0; j < vec.size(); j++) {
            cout << vec[j] << " ";
        }
        cout << endl;
    }
};

//pseudocode: m = numRows, n = numCols, x = solnVector

//do until solution is valid or invalid:
    //initialize J as empty set
    //go through each row 1..m
        //if row has one with exactly one minimum outside J
        //add row number to J
    //if |J| == 1, success (do some steps and return x)
    //if |J| == n, fail (no solution as demonstrated by [dis]provable by contradiction)
    //otherwise, increment values in x as appropriate
SolnVector Grigoriev_Algorithm(vector<vector<int>> &matrix);
SolnVector AGG_Algorithm(vector<vector<int>> &matrix);
void solve_Grigoriev(Matrix &matrix, SolnVector &x);
void solve_AGG(Matrix &matrix, SolnVector &x);


int main(int argc, char *argv[]) {
    vector<vector<int>> input_matrix = retrieveMatrix(); //input
    //SolnVector x1 = Grigoriev_Algorithm(input_matrix);
    SolnVector x2 = AGG_Algorithm(input_matrix);
//     vector<vector<int>> vec_matrix = retrieveMatrix(); //input
//     Matrix matrix(vec_matrix);

//     int m = matrix.rows(), n = matrix.cols();
//     SolnVector x(n); //solution vector initialized with zeroes

//     while (true) {
//         unordered_set<int> J;
//         for (int i = 1; i <= m; i++) {
//             //cout << "J SIZE2 = " << J.size() << endl;
//             int curMin = INT_MAX;
//             unordered_set<int> curMinColsNotInJ;
//             for (int j = 1; j <= n; j++) {
//                 int curVal = matrix.get(i, j) + x.get(j);
//                 if (curVal < curMin) {
//                     curMin = curVal;
//                     curMinColsNotInJ.clear();
//                     insertIfNotInJ(j, J, curMinColsNotInJ);

//                 }
//                 else if (curVal == curMin) {
//                     insertIfNotInJ(j, J, curMinColsNotInJ);
//                     //cout << curMinColsNotInJ.size() << ", " << J.size() << endl;
//                 }
//             }
//             if (curMinColsNotInJ.size() == 1) {
//                 for (int j : curMinColsNotInJ) {
//                     //cout << "inserting: " << j << endl;
//                     J.insert(j); //extend J with j, continue
//                 }
//             }
//         }

//         int k = J.size();
//         if (k == n) {
//             cout << "there is no solution for this input" << endl;
//             return EXIT_SUCCESS;
//         }
//         if (k == 0) { //is this an end case?
//             x.printSolution();
//             return EXIT_SUCCESS;
//         }
//         if (k == 1) {

// //ASK ABOUT THIS
// // Now, let k = |J| < n. If k = 1, we add to x1 the number
// // min2≤j≤n{a1,j +xj}−(a1,1 +x1) ≥ 1 and obtain a solution of (0.1).
// // Thereupon, we apply Lemma 1.2 to the obtained solution; as a
// // result, the algorithm terminates the inductive step


// //J contains the only column that has a row (single) with a minimum. We just need to find the row and  
// //subtract index j's value from the second minimum (not including index j in the row). We then add that to x[j].
//             int finalJ;
//             for (int row : J) finalJ = row; //weird syntax, but must happen exactly once because k == 1
//             for (int i = 1; i <= m; i++) {
//                 int jVal = matrix.get(i, finalJ) + x.get(finalJ);
//                 bool singleMin = true; //signifies whether the row has a single minimum
//                 for (int j = 1; singleMin && j <= n; j++) {
//                     if (j != finalJ && jVal >= (matrix.get(i, j)) + x.get(j)) singleMin = false;       
//                 }
//                 if (singleMin) {
//                     int secondMin = INT_MAX;
//                     for (int j = 1; j <= n; j++) {
//                         if (j != finalJ) secondMin = min(secondMin, matrix.get(i, j) + x.get(j));      
//                     }
//                     x.add(finalJ, secondMin - jVal); //achieve final equality
//                     x.printSolution();
//                     return EXIT_SUCCESS;
//                 }
//             }
//             cout << "Error: no solution found in solution branch" << endl;
//             return EXIT_FAILURE;
//         }

//         //get attracted rows so far
//         unordered_set<int> attractedRows;
//         for (int i = 1; i <= m; i++) {
//             unordered_set<int> minCols;
//             int rowMinElem = INT_MAX;
//             for (int j = 1; j <= n; j++) {
//                 int elem = matrix.get(i, j) + x.get(j);
//                 if (elem < rowMinElem) {
//                     minCols = {j};
//                     rowMinElem = elem;
//                 }
//                 else if (elem == rowMinElem) {
//                     minCols.insert(j);
//                 }
//             }
//             bool valid = true;
//             for (int minCol : minCols) { //this is inefficient, but it's okay for now
//                 if (J.find(minCol) == J.end()) {
//                     valid = false;
//                     break;
//                 }
//             }
//             if (valid) attractedRows.insert(i);
//         }
//         //cout << "aRows size: " << attractedRows.size() << endl;

//         //update x as necessary before continuing
//         int a = INT_MAX, leftMin = INT_MAX, rightMin = INT_MAX;
//         for (int i : attractedRows) {
//             for (int j = 1; j <= n; j++) {
//                 if (J.find(j) == J.end()) { //cols not in J, for leftMin calculation
//                     leftMin = min(leftMin, matrix.get(i, j) + x.get(j));
//                 }
//                 rightMin = min(rightMin, matrix.get(i, j) + x.get(j));
//             }
//             assert(leftMin - rightMin >= 1); //according to paper this must be true
//             a = min(a, leftMin - rightMin);
//         }
//         for (int j : J) x.add(j, a); //add 'a' to x at indices in j
//     }

    return 0;
}

void solve_Grigoriev(Matrix &matrix, SolnVector &x) {
    cout << "before: ";
    x.printSolution();
    cout << endl;
    int m = matrix.rows(), n = matrix.cols();
    while (true) {
        unordered_set<int> J;
        for (int i = 1; i <= m; i++) {
            //cout << "J SIZE2 = " << J.size() << endl;
            int curMin = INT_MAX;
            unordered_set<int> curMinColsNotInJ;
            for (int j = 1; j <= n; j++) {
                int curVal = matrix.get(i, j) + x.get(j);
                if (curVal < curMin) {
                    curMin = curVal;
                    curMinColsNotInJ.clear();
                    insertIfNotInJ(j, J, curMinColsNotInJ);

                }
                else if (curVal == curMin) {
                    insertIfNotInJ(j, J, curMinColsNotInJ);
                    //cout << curMinColsNotInJ.size() << ", " << J.size() << endl;
                }
            }
            if (curMinColsNotInJ.size() == 1) {
                for (int j : curMinColsNotInJ) {
                    //cout << "inserting: " << j << endl;
                    J.insert(j); //extend J with j, continue
                }
            }
        }

        int k = J.size();
        if (k == n) {
            cout << "there is no solution for this input" << endl;
            return;
        }
        if (k == 0) { //solved already
            return;
        }
        if (k == 1) {

//ASK ABOUT THIS
// Now, let k = |J| < n. If k = 1, we add to x1 the number
// min2≤j≤n{a1,j +xj}−(a1,1 +x1) ≥ 1 and obtain a solution of (0.1).
// Thereupon, we apply Lemma 1.2 to the obtained solution; as a
// result, the algorithm terminates the inductive step


//J contains the only column that has a row (single) with a minimum. We just need to find the row and  
//subtract index j's value from the second minimum (not including index j in the row). We then add that to x[j].
            int finalJ;
            for (int row : J) finalJ = row; //weird syntax, but must happen exactly once because k == 1            //x.add(finalRow) += minRange(i, 2, n) - (matrix.get(1, 1) + 0); //TODO fix

            for (int i = 1; i <= m; i++) {
                int jVal = matrix.get(i, finalJ) + x.get(finalJ);
                bool singleMin = true; //signifies whether the row has a single minimum
                for (int j = 1; singleMin && j <= n; j++) {
                    if (j != finalJ && jVal >= (matrix.get(i, j)) + x.get(j)) singleMin = false;       
                }
                if (singleMin) {
                    int secondMin = INT_MAX;
                    for (int j = 1; j <= n; j++) {
                        if (j != finalJ) secondMin = min(secondMin, matrix.get(i, j) + x.get(j));      
                    }
                    x.add(finalJ, secondMin - jVal); //achieve final equality
                    cout << "after: ";
                    x.printSolution();
                    cout << endl;
                    return;
                }
            }
            cout << "Error: no solution found in solution branch" << endl;
            return;
        }

        //get attracted rows so far
        unordered_set<int> attractedRows;
        for (int i = 1; i <= m; i++) {
            unordered_set<int> minCols;
            int rowMinElem = INT_MAX;
            for (int j = 1; j <= n; j++) {
                int elem = matrix.get(i, j) + x.get(j);
                if (elem < rowMinElem) {
                    minCols = {j};
                    rowMinElem = elem;
                }
                else if (elem == rowMinElem) {
                    minCols.insert(j);
                }
            }
            bool valid = true;
            for (int minCol : minCols) { //this is inefficient, but it's okay for now
                if (J.find(minCol) == J.end()) {
                    valid = false;
                    break;
                }
            }
            if (valid) attractedRows.insert(i);
        }
        //cout << "aRows size: " << attractedRows.size() << endl;

        //update x as necessary before continuing
        int a = INT_MAX, leftMin = INT_MAX, rightMin = INT_MAX;
        for (int i : attractedRows) {
            for (int j = 1; j <= n; j++) {
                if (J.find(j) == J.end()) { //cols not in J, for leftMin calculation
                    leftMin = min(leftMin, matrix.get(i, j) + x.get(j));
                }
                rightMin = min(rightMin, matrix.get(i, j) + x.get(j));
            }
            assert(leftMin - rightMin >= 1); //according to paper this must be true
            a = min(a, leftMin - rightMin);
        }
        for (int j : J) x.add(j, a); //add 'a' to x at indices in j
    }
}

void solve_AGG(Matrix &matrix, SolnVector &x) {
    int m = matrix.rows(), n = matrix.cols();
    while (true) {
        unordered_map<int, int> toAdd; //maps rows to lifting amount
        for (int i = 1; i <= m; i++) {
            pair<int, int> minimum_idx = {INT_MAX, -1};
            bool unique = true;
            for (int j = 1; j <= n; j++) {
                int cur = matrix.get(i, j) + x.get(j);
                if (cur == minimum_idx.first) {
                    unique = false;
                }
                else if (cur < minimum_idx.first) {
                    unique = true;
                    minimum_idx = {cur, j};
                }
            }
            if (unique) {
                int nextMinimum = INT_MAX;
                for (int j = 1; j <= n; j++) {
                    if (j == minimum_idx.second) continue;
                    int cur = matrix.get(i, j) + x.get(j);
                    nextMinimum = min(nextMinimum, cur);
                }
                int amtToLift = nextMinimum - minimum_idx.first; //must be positive
                int colToLift = minimum_idx.second;
                if (toAdd.find(colToLift) == toAdd.end()) toAdd[colToLift] = amtToLift;
                else toAdd[colToLift] = min(toAdd[colToLift], amtToLift);
            }
        }

        if (toAdd.empty()) {
            cout << "printing AGG soln: ";
            x.printSolution();
            return;
        } //end condition

        for (const auto &it : toAdd) {
            int idx = it.first, amt = it.second;
            x.add(idx, amt);
        }
    }
}

SolnVector Grigoriev_Algorithm(vector<vector<int>> &matrix) {
    SolnVector s(matrix[0].size());
    vector<vector<int>> emptyMatrix = {};
    Matrix m(emptyMatrix);

    for (vector<int> r : matrix) {
        m.addRow(r);
        cout << "before0, with cursize " << m.rows() << ": ";
        s.printSolution();
        cout << endl;
        solve_Grigoriev(m, s);
        cout << "after0: ";
        s.printSolution();
        cout << endl;
        // s.printSolution();
    }
    // printSolution(s);
    return s;
}

SolnVector AGG_Algorithm(vector<vector<int>> &matrix) {
    SolnVector s(matrix[0].size());
    Matrix m(matrix);
    solve_AGG(m, s);
    return s;
}






//TODO define big intereger class to get ints of size 2^n  (or size b^n specifially)
//also define randomizer for matrix of size m * n

// Akian, Gaubert, Guterman
// Davydowv
