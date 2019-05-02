#ifndef MMTBX_TLS_UTILS_H
#define MMTBX_TLS_UTILS_H

#include <string>

// Basic data types
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <scitbx/sym_mat3.h>

// Allow arrays of the above
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/accessors/c_grid.h>

#include <mmtbx/tls/decompose.h>

namespace mmtbx { namespace tls { namespace utils {

// These are used a lot...
namespace af = scitbx::af;
namespace dec = mmtbx::tls::decompose;

// Commonly used types
typedef scitbx::vec3<double> vec;
typedef scitbx::mat3<double> mat;
typedef scitbx::sym_mat3<double> sym;
// 1-dimensional arrays
typedef af::shared<double> dblArr1d;
typedef af::shared<size_t> selArr1d;
typedef af::shared<vec>    vecArr1d;
typedef af::shared<sym>    symArr1d;
// N-dimensional arrays
typedef af::versa<double, af::flex_grid<> > dblArrNd;
typedef af::versa<vec,    af::flex_grid<> > vecArrNd;
typedef af::versa<sym,    af::flex_grid<> > symArrNd;

typedef enum
{
  TLSBlankString = 0,
  TLSTranslation = 1 << 0,
  TLSLibration   = 1 << 1,
  TLSScrew       = 1 << 2,
} TLSComponent;

vecArrNd uij_eigenvalues(const symArrNd& uijs);

//! Convenience class for the manipulation, validation and usage of TLS matrices
/*! Convenience class containing functions for:
 *    Adding TLS matrices
 *    Multiplying TLS matrices
 *    Running TLS-validation routines
 *    Generating Uij from TLS Matrices
 */
class TLSMatrices {

  public:
    //! Constructor from nothing
    TLSMatrices();

    //! Constructor from separate matrices
    TLSMatrices(const sym &T, const sym &L, const mat &S);

    //! Constructor from another TLSMatrices object
    TLSMatrices(const TLSMatrices &other);

    //! Constructor from array of 21 values
    TLSMatrices(const dblArr1d &values);

    //! Get number of decimal places for rounding values
    static int getPrecision();

    //! Set number of decimal places for rounding values
    static void setPrecision(int decimals);

    //! Get default tolerance for determining validity of matrices
    static double getTolerance();

    //! Set default tolerance for determining validity of matrices
    static void setTolerance(double tolerance);

    //! Override + operator for adding two classes together.
    TLSMatrices& operator+(const TLSMatrices &other) const;

    //! Override * operator for multiplying by a scalar
    TLSMatrices& operator*(double scalar) const;

    //! Add another TLSMatrices object to this one, matrix by matrix.
    void add(const TLSMatrices &other);

    //! Test if the matrix values of any of the selected components are nonzero
    bool any(const std::string & component_string = "TLS", double tolerance = -1);

    //! Create a copy of the class
    TLSMatrices* copy() const;

    //! Return a class that determines the physical motions associated with the TLS matrices
    dec::decompose_tls_matrices decompose(double tolerance = -1) const;

    //! Extract matrix values as an array, given an enum 0-7 which selects combinations of T(1), L(2) & S(4) by bitwise operations.
    dblArr1d getValuesByInt(const TLSComponent &components, bool include_szz = true) const;

    //! Extract matrix values as an array, given a selection string containing combinations of T,L or S.
    dblArr1d getValuesByString(const std::string & component_string = "TLS", bool include_szz = true) const;

    //! Determine whether the TLS matrices are physically valid by decomposition
    bool isValid(double tolerance = -1) const;

    //! Scale all matrix values by a scalar
    void multiply(double scalar);

    //! Scale the TLS matrices so that the generated Uij are on a particular scale.
    /*!
     * On failure, return negative value
     */
    double normalise(const vecArr1d &sites_cart, const vec &origin, double target = 1.0, double tolerance = -1);

    //! Determine the number of parameters (total or free)
    int paramCount(bool free = true, bool non_zero = false);

    //! Set all matrix elements to 0
    void reset();

    //! Set matrices values from an array and a matrix selection enum (0-7)
    void setValuesByInt(const dblArr1d &values, const TLSComponent &components, bool include_szz = true);

    //! Set matrices vaulues from an array and a matrix selection string (containing combinations of T, L or S)
    void setValuesByString(const dblArr1d &values, const std::string &component_string = "TLS", bool include_szz = true);

    //! Generates summary string for the class
    std::string summary();

    //! Convert TLS matrices to Uijs from coordinates and an origin
    symArr1d uijs(const vecArr1d &sites_cart, const vec &origin) const;

    //! Returns T-matrix
    sym getT();
    //! Returns L-matrix
    sym getL();
    //! Returns S-matrix
    mat getS();

  private:
    //! T-, L- & S-matrix
    sym T;
    sym L;
    mat S;

    static double tol; //!< Tolerance determining when a number is considered non-zero
    static double rnd; //!< Number of decimals that the matrices will be rounded to

    //! Test if the matrix values of any of the selected components are nonzero (using enum 0-7)
    bool anyByInt(const TLSComponent &components, double tolerance = -1);

    //! Format a sym_mat3 matrix neatly to buffer
    std::string matrix_to_string(const sym &M);

    //! Format a mat3 matrix neatly to buffer
    std::string matrix_to_string(const mat &M);

    //! Round matrix values to the set precision
    void round();

    //! Set a tolerance value to the default value if less than zero
    void sanitiseTolerance(double *tolerance) const;

    //! Scale all matrices by a scalar
    void scale(double mult);

    //! Scale a block of values by a scalar
    void scaleComponent(double *comp, int size, double mult);

    //! Set the trace of the S-matrix by changing Szz
    void setSzzValueFromSxxAndSyy(double target_trace = 0);

    //! Convert a string to TLSComponents enum based on whether it contains the letters T, L or S
    TLSComponent stringToComponents(const std::string &component_string) const;

};

//! Overload that allows left-multiplication of TLSMatrices/TLSAmplitudes
TLSMatrices operator*(const double &v, const TLSMatrices &x);

//! Class containing a list of amplitudes for multiplying onto TLSMatrices
/*! Takes a TLSMatrices object and returns a series of TLSMatrices objects,
 * where each of the TLSMatrices is multiplied by a different amplitude.
 */
class TLSAmplitudes {

  public:
    //! Constructor from number of amplitudes
    TLSAmplitudes(size_t n);

    //! Constructor from another TLSAmplitudes object
    TLSAmplitudes(const TLSAmplitudes &other);

    //! Constructor from array of values
    TLSAmplitudes(const dblArr1d &values);

    //! Get the number of decimal places to be used for rounding
    static int getPrecision();

    //! Set the number of decimal places to be used for rounding
    static void setPrecision(int decimals);

    //! Get default tolerance for determining validity of matrices
    static double getTolerance();

    //! Set default tolerance for determining validity of matrices
    static void setTolerance(double tolerance);

    //! Indexing operator returns the i'th element of the amplitudes
    double operator[] (const int index);

    //! Override + operator for adding two classes together.
    TLSAmplitudes& operator+(const TLSAmplitudes &other) const;

    //! Override * operator for multiplying by a scalar
    TLSAmplitudes& operator*(double scalar) const;

    //! Addition of another class
    void add(const TLSAmplitudes &other);

    //! Test if any amplitudes values are nonzero
    bool any(double tolerance = -1);

    //! Create a copy of the class
    TLSAmplitudes* copy() const;

    //! Description strings of the functions
    std::string getDescription();

    //! Get the amplitude values (by value!)
    dblArr1d getValues() const;

    //! Get the amplitude values for a selection of indices (by value!)
    dblArr1d getValuesBySelection(const selArr1d &selection);

    //! Multiplication by a scaler
    void multiply(double scalar);

    //! Divide the amplitudes by their average, and return multiplier to be applied to the TLSMatrices
    double normalise(double target = 1.0);

    //! Count the number of (non-zero) parameters
    int paramCount(bool non_zero = false);

    //! Set all amplitudes to 1.0
    void reset();

    //! Round all amplitudes to a number of decimal places
    void round();

    //! Set matrices vals from an array of the same size as the object
    void setValues(const dblArr1d &values);

    //! Set matrices vals from an array and a selection of indices
    void setValuesBySelection(const dblArr1d &values, const selArr1d &selection);

    //! Return the number of amplitudes values
    int size() const;

    //! Get summary string for the amplitudes in the object
    std::string summary();

    void zeroValues();

    void zeroNegativeValues();

  private:
    //! Array of amplitudes
    dblArr1d vals;

    static double tol; //!< Tolerance determining when a number is considered non-zero
    static double rnd; //!< Number of decimals that the matrices will be rounded to

    // Description strings
    static std::string short_description;
    static std::string description;

    //! Set a tolerance value to the default value if less than zero
    void sanitiseTolerance(double *tolerance);

    void scale(double multiplier);

    void validateSelection(const selArr1d &selection);

};

//! Overload that allows left-multiplication of TLSAmplitudes
TLSAmplitudes operator*(const double &v, const TLSAmplitudes &x);

class TLSMatricesAndAmplitudes {

  public :
    //! Constructor from number of amplitudes
    TLSMatricesAndAmplitudes(size_t n);

    //! Constructor from other class
    TLSMatricesAndAmplitudes(const TLSMatricesAndAmplitudes& other);

    //! Constructor from matrix-amplitude pair
    TLSMatricesAndAmplitudes(TLSMatrices& matrices, TLSAmplitudes& amplitudes);

    //! Constructor from matrix and amplitude arrays
    TLSMatricesAndAmplitudes(const dblArr1d& matrix_values, const dblArr1d& amplitude_values);

    //! Create a copy of the class
    TLSMatricesAndAmplitudes* copy() const;

    //! Return a reference to the owned matrices class
    TLSMatrices* getMatrices();
    const TLSMatrices* getMatricesConst() const;

    //! Return a reference to the owned amplitudes class
    TLSAmplitudes* getAmplitudes();
    const TLSAmplitudes* getAmplitudesConst() const;

    //! Return a list of amplitude-multiplied TLSMatrices
    af::shared<TLSMatrices> expand();

    //! Return a list of amplitude-multiplied TLSMatrices by selection
    af::shared<TLSMatrices> expand(const selArr1d &selection);

    bool isValid(double tolerance = -1);

    bool isValid(const selArr1d &selection, double tolerance = -1);

    bool isNull(double matricesTolerance = -1, double amplitudesTolerance = -1);

    int paramCount(bool free = true, bool non_zero = false);

    double normaliseByAmplitudes(double target = 1.0);

    double normaliseByMatrices(const vecArrNd &sites_carts, const vecArr1d &origins, double target = 1.0);

    void reset();

    void setLabel(int label);

    std::string summary();

    symArrNd uijs(const vecArrNd &sites_carts, const vecArr1d &origins);

    symArrNd uijs(const vecArrNd &sites_carts, const vecArr1d &origins, const selArr1d &selection);

    void resetIfNull(double matricesTolerance = -1, double amplitudesTolerance = -1);

  private:
    TLSMatrices* mats;
    TLSAmplitudes* amps;

    int lbl;

    symArrNd uijs(const vecArrNd &sites_carts, const vecArr1d &origins, const af::shared<TLSMatrices> &tls_matrices);

};

class TLSMatricesAndAmplitudesList {

  public:
    //! Constructor from length and number of amplitudes
    TLSMatricesAndAmplitudesList(size_t length, size_t n_amplitudes);

    //! Constructor for existing TLSMatricesAndAmplitudesList object
    TLSMatricesAndAmplitudesList(const TLSMatricesAndAmplitudesList &other);

    //! Constructor from two arrays of matrices/amplitudes
    TLSMatricesAndAmplitudesList(const dblArrNd& matrix_values, const dblArrNd& amplitude_values);

    TLSMatricesAndAmplitudes* operator[] (const int index);

    TLSMatricesAndAmplitudesList* copy() const;

    TLSMatricesAndAmplitudes* get(const int index);

    TLSMatricesAndAmplitudes* getConst(const int index) const;

    bool isNull(double matricesTolerance = -1, double amplitudesTolerance = -1);

    int paramCount(bool free = true, bool non_zero = false);

    void normaliseByAmplitudes(double target = 1.0);

    void normaliseByMatrices(const vecArrNd &sites_carts, const vecArr1d &origins, double target = 1.0);

    void reset();

    void resetMatrices();

    void resetNullModes(double matricesTolerance = -1, double amplitudesTolerance = -1);

    size_t size() const;

    std::string summary();

    symArrNd uijs(const vecArrNd &sites_carts, const vecArr1d &origins);

    symArrNd uijs(const vecArrNd &sites_carts, const vecArr1d &origins, const selArr1d &selection);

    void zeroAmplitudes(const selArr1d &selection);

    void zeroNegativeAmplitudes();

  private:
    af::shared<TLSMatricesAndAmplitudes*> list;

    void initialiseList(size_t length, size_t n_amplitudes);

    void validateIndex(size_t index) const;

    void validateSelection(const selArr1d &selection);

};

} } } // close namepsace mmtbx/tls/utils

#endif
