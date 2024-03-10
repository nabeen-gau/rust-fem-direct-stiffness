//! Personal Matrix Calulation Library
//!

use std::fmt;

/// Calculates the determinant of given matrix by gauss jordan elimination
/// TODO: add pivoting
fn determinant(matrix: Matrix) -> f64 {
    let mut det = 0.0;

    if matrix.columncount == 1 {
        return matrix.items[0];
    }

    for i in 0..matrix.columncount {
        let new_size = matrix.items.len() - matrix.columncount * 2 + 1;

        if new_size == 1 {
            return matrix.items[0] * matrix.items[3] - matrix.items[1] * matrix.items[2];
        }

        let mut submatrix = Matrix::new(vec![0.0; new_size]);
        let mut index = 0;
        for j in 0..matrix.columncount {
            for k in 0..matrix.columncount {
                if k != i && j != 0 {
                    // println!("{}{}", j,k);
                    submatrix.items[index] = matrix.items[j * matrix.columncount + k];
                    index += 1;
                }
            }
        }

        let col = i % matrix.columncount;
        let row = (i as f64 / matrix.columncount as f64) as usize;
        if (row + col) % 2 == 0 {
            det += matrix.items[i] * determinant(submatrix);
        } else {
            det -= matrix.items[i] * determinant(submatrix);
        }
    }
    det
}

/// Adds identity matrix to given matrix as augmented matrix
/// Example:
/// 1   2
/// 3   4
/// changes to
/// 1   2   1   0
/// 3   4   0   1
fn add_identity(mut matrix: Matrix) -> Matrix {
    let size = matrix.columncount;

    for i in (0..size).rev() {
        let mut unity = vec![0.0; size];
        unity[i] = 1.0;
        matrix.items.splice((i+1)*matrix.columncount..(i+1)*matrix.columncount, unity);

    }
    matrix.columncount *= 2;
    matrix
}

// fn get_rc(matrix: Matrix, index: usize) -> (usize, usize) {
//     let col = index % matrix.columncount;
//     let row = (index as f64 / matrix.columncount as f64) as usize;
//     (row, col)
// }

/// performs gauss jordan elimination on given matrix
/// todo: add pivoting
fn gauss_jordan(mut matrix: Matrix) -> Matrix {
    for i in 0..matrix.rowcount(){
        let diag = matrix.get_by_rc(i, i);
        for j in 0..matrix.columncount {
            let val = matrix.get_by_rc(i, j);
            matrix.set_by_rc(i, j, val/diag);
        }

        for k in 0..matrix.rowcount() {
            if k!= i {

                let factor = matrix.get_by_rc(k, i);
                for l in 0..matrix.columncount {
                    if i!=k{
                        let val = matrix.get_by_rc(k, l);
                        let mainval = matrix.get_by_rc(i, l);
                        matrix.set_by_rc(k, l, val-factor*mainval);
                    }
                }
            }
        }
    }

    
    matrix
}

/// Matrix struct
/// Example
/// Matrix::new(vec![1, 2, 3, 4]) to create this struct with following elements
/// 1   2
/// 3   4
#[derive(Debug)]
pub struct Matrix {
    pub columncount: usize,
    pub items: Vec<f64>,
}

impl fmt::Display for Matrix{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // write!(f, "")?;
        for (i, s) in self.items.iter().enumerate() {
            if i > 0 {
                write!(f, "\t")?;
            }
            if i%self.columncount == 0 {
                write!(f, "\n")?;
            }
            // let n = format!("{:.2e}", s);
            if s < &1000.0 && s > &-1000.0{
                write!(f, "{}", format!("{:.2}", s))?;
            }
            else {
                write!(f, "{}", format!("{:.2e}", s))?;
            }
        }
        write!(f, "")
    }
}

impl Matrix {
    /// Create a square matrix with given data as elements row-wise
    /// # Example
    /// use ds::Matrix;
    /// let c = Matrix::new(vec![1.0, 2.0, 0.0, 4.0, 0.0, 1.0, 3.0, 5.0, 2.0]);
    /// creates a 3x3 matrix with elements as follows
    /// 1.0     2.0     0.0
    /// 4.0     0.0     1.0
    /// 3.0     5.0     2.0
    ///
    pub fn new(data: Vec<f64>) -> Matrix {
        Matrix {
            columncount: (data.len() as f64).sqrt() as usize,
            items: data,
        }
    }

    /// size is the order of matrix example 6x6 matrix size=6
    pub fn zeros(size: usize) -> Matrix {
        Matrix {
            columncount: size,
            items: vec![0.0;size*size],
        }
    }

    pub fn clone(&self) -> Matrix {
        Matrix {
            columncount: self.columncount,
            items: self.items.clone(),
        }
    }

    /// return the no. of rows of the matrix
    fn rowcount(&self) -> usize {
        self.items.len()/self.columncount
    }

    /// get the element of matrix by row, column index
    fn get_by_rc(&self, row:usize, col:usize) -> f64 {
        self.items[row*self.columncount+col]
    }

    /// set the element of matrix by row, column index
    fn set_by_rc(&mut self, row:usize, col:usize, value: f64) {
        self.items[row*self.columncount+col] = value;
    }


    /// multipy the matrix with a given matrix of correct size
    pub fn multipy(&self, matrix: Matrix) -> Matrix {
        let mut result = Matrix::new(vec![0.0; self.items.len()]);
        for i in 0..self.columncount {
            for j in 0..self.columncount {
                result.items[i * self.columncount + j] = 0.0;
                for k in 0..self.columncount {
                    result.items[i * self.columncount + j] += self.items[i * self.columncount + k]
                        * matrix.items[k * self.columncount + j];
                }
            }
        }
        result
    }

    /// multipy the matrix with a appropriate sized vector
    pub fn vec_multipy(&self, vector: Vec<f64>) -> Vec<f64> {
        let mut product = vec![0.0;vector.len()];
        for vec_index in 0..product.len() {
            for mat_index in 0..self.columncount {
                product[vec_index] += self.get_by_rc(vec_index, mat_index)*vector[mat_index]; 
            }
        }
        product
    }

    /// Returns transpose of the matrix
    pub fn transpose(&self) -> Matrix {
        let mut transpose = Matrix::new(vec![0.0; self.items.len()]);
        
        for i in 0..self.columncount {
            for j in 0..self.columncount {
                transpose.set_by_rc(i, j, self.get_by_rc(j, i))
            }
        }

        transpose
    }

    /// returns the inverse of matrix if possible
    pub fn inverse(&self) -> Matrix {
        let mut matrix = gauss_jordan(add_identity(self.clone()));
        matrix.columncount /= 2;

        let mut count = 0;
        let mut removeflag = false;
        for i in (0..matrix.items.len()).rev() {
            if removeflag {
                matrix.items.remove(i);
            }
            count += 1;
            if count == matrix.columncount {
                count = 0;
                removeflag = !removeflag;
            }
        }
        matrix
    }


    /// returns the determinant if possible
    pub fn determinant(&self) -> f64 {
        determinant(self.clone())
    }

    /// swap two rows of the matrix
    /// todo: for implementation of pivoting in gauss elimination
    /// not sure it is needed here
    fn swap_row(&mut self, row1: usize, row2: usize) {
        let mat_copy = self.items.clone();

        let start_id = row1 * self.columncount;
        let end_id = start_id + self.columncount;

        let start_id2 = row2 * self.columncount;
        let end_id2 = start_id2 + self.columncount;

        for (i, j) in (start_id2..end_id2).zip(start_id..end_id) {
            self.items[i] = mat_copy[j];
        }

        for (i, j) in (start_id2..end_id2).zip(start_id..end_id) {
            self.items[j] = mat_copy[i];
        }
    }

    pub fn remove_rowcol(&mut self, index: usize) {

        let col = self.columncount;
        for j in (index+col+index*col..index+col*(col-1)+1).step_by(col).rev() {
            self.items.remove(j);
        }
        for i in (col*index..col*index+col).rev() {
            self.items.remove(i);
        }
        for j in (index..index+index*col).step_by(col).rev() {
            self.items.remove(j);
        }

        self.columncount -= 1;

    }

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matrix_multiplication() {
        let b = Matrix {
            columncount: 3,
            items: vec![1.0, 3.0, 3.0, 4.0, 5.0, 1.0, 2.0, 3.0, 4.0],
        };

        let c = Matrix {
            columncount: 3,
            items: vec![1.0, 2.0, 0.0, 4.0, 0.0, 1.0, 3.0, 5.0, 2.0],
        };

        let mul_matrix = b.multipy(c);
        assert_eq!(
            mul_matrix.items,
            vec![22.0, 17.0, 9.0, 27.0, 13.0, 7.0, 26.0, 24.0, 11.0]
        )
    }

    #[test]
    fn determinant_test() {
        let b = Matrix {
            columncount: 3,
            items: vec![1.0, 3.0, 3.0, 4.0, 5.0, 1.0, 2.0, 3.0, 4.0],
        };
        assert_eq!(b.determinant(), -19.0);
    }

    #[test]
    fn swap_row_test() {
        let mut b = Matrix {
            columncount: 3,
            items: vec![1.0, 3.0, 3.0, 4.0, 5.0, 1.0, 2.0, 3.0, 4.0],
        };
        b.swap_row(1, 2);
        assert_eq!(
            b.items,
            vec![
                1.0, 3.0, 3.0,
                2.0, 3.0, 4.0,
                4.0, 5.0, 1.0,
            ]
        )
    }

    #[test]
    fn inverse_test() {

        let b = Matrix {
            columncount: 3,
            items: vec![1.0, 3.0, 3.0, 4.0, 5.0, 1.0, 2.0, 3.0, 4.0],
        };
        let inv = vec![-0.8947368421052632, 0.15789473684210525, 0.631578947368421, 0.736842105263158, 0.10526315789473684, -0.5789473684210525, -0.10526315789473688, -0.15789473684210525, 0.3684210526315789];

        assert_eq!(b.inverse().items, inv);
    }

    #[test]
    fn vec_multipy() {
        let d = Matrix::new(vec![1.0, 2.0, 3.0, 4.0]);

        assert_eq!(d.vec_multipy(vec![3.0, -6.0]), vec![-9.0, -15.0]);
    }

    #[test]
    fn transpose() {
        let d = Matrix::new(vec![1.0, 2.0, 3.0, 4.0]);

        assert_eq!(d.transpose().items, vec![1.0, 3.0, 2.0, 4.0]);
    }

    #[test]
    fn remove_rowcol_test() {
        let mut d = Matrix::new(vec![1.0, 2.0, 3.0, 4.0, 
                                    5.0, 6.0, 7.0, 8.0,
                                    9.0, 10.0, 11.0, 12.0,
                                    13.0, 14.0, 15.0, 16.0]);
        let index = 2;
        d.remove_rowcol(index);

        assert_eq!(d.items, vec![1.0, 2.0, 4.0,
                                5.0, 6.0, 8.0,
                                13.0, 14.0, 16.0]);
    }
}
