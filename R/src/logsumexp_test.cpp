#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

VectorXd logsumexp(const MatrixXd x, int axis) {
  if (axis == 0) { // Compute log-sum-exp along rows
    VectorXd maxCoeffs = x.rowwise().maxCoeff();
    
    return maxCoeffs.array() + (x.colwise() - maxCoeffs).array().exp().rowwise().sum().log();
  } else if (axis == 1) { // Compute log-sum-exp along columns
    VectorXd maxCoeffs = x.colwise().maxCoeff();
    
    return maxCoeffs.transpose().array() + (x.rowwise() - maxCoeffs.transpose()).array().exp().colwise().sum().log();
  } else {
    // Invalid axis value
    std::cerr << "Invalid axis value. Use 0 for column-wise or 1 for row-wise." << std::endl;
    return VectorXd::Zero(x.rows());
  }
  
  return VectorXd::Zero(x.rows());
}

// Expose the Sinkhorn function to R
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
int main() {
  // Example usage
  MatrixXd inputMatrix(3, 4);
  inputMatrix << 1.0, 2.0, 3.0, 4.0,
                 5.0, 6.0, 7.0, 8.0,
                 9.0, 10.0, 11.0, 12.0;
  
  VectorXd resultRowWise = logsumexp(inputMatrix, 0); // Row-wise
  VectorXd resultColWise = logsumexp(inputMatrix, 1); // Column-wise
  
  std::cout << "Matrix:\n" << inputMatrix << std::endl;
  std::cout << "Row-wise result:\n" << resultRowWise << std::endl;
  std::cout << "Column-wise result:\n" << resultColWise << std::endl;
  
  return 0;
}

/***R
main()
*/