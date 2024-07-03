%let pgm=utl-challenges-associated-with-hilbert-matrices-in-regression-design;

Challenges associated with hilbert matrices in regression design


  Two Sections

      1 hilbert regression
      2 independence

github
https://tinyurl.com/mtxavtpk
https://github.com/rogerjdeangelis/utl-challenges-associated-with-hilbert-matrices-in-regression-design/tree/main

/*               _     _
 _ __  _ __ ___ | |__ | | ___ _ __ ___
| `_ \| `__/ _ \| `_ \| |/ _ \ `_ ` _ \
| |_) | | | (_) | |_) | |  __/ | | | | |
| .__/|_|  \___/|_.__/|_|\___|_| |_| |_|
|_|
*/

/**************************************************************************************************************************/
/*                                                                                                                        */
/*  Reproduce the known linear relationship below using regression                                                        */
/*                                                                                                                        */
/*  Y = 2*X1 + 3*X2 -1.5*X3 + 0.5*X4 -X5  (note equation goes through the origin)                                         */
/*                                                                                                                        */
/*  INPUT                                                                                                                 */
/*                                            Use this arbitrary equation                                                 */
/*                   Hilbert Rows             to create the Y                Regr Err      Y+Err                          */
/*        X1      X2      X3      X4     X5   Y=2*X1+3*X2-1.5*X3+0.5*X4-X5   normal(0,.1)  Yhat                           */
/*   ---------------------------------------  ----------------------------   ------------  --------                       */
/*   1.00000 0.50000 0.33333 0.25000 0.20000         2.92500                 -0.02677       3.04907                       */
/*   0.50000 0.33333 0.25000 0.20000 0.16667         1.55833                  0.03448       1.45694                       */
/*   0.33333 0.25000 0.20000 0.16667 0.14286         1.05714                 -0.08589       1.02468                       */
/*   0.25000 0.20000 0.16667 0.14286 0.12500         0.79643                 -0.10936       0.72875                       */
/*   0.20000 0.16667 0.14286 0.12500 0.11111         0.63710                 -0.01505       0.67609                       */
/*                                                                                                                        */
/*  In this case we know the linear 'best' coeficients before regression                                                  */
/*                                                                                                                        */
/*  Known best coeficients                                                                                                */
/*                                                                                                                        */
/*   X1            2                                                                                                      */
/*   X2            3                                                                                                      */
/*   X3          1.5                                                                                                      */
/*   X4          0.5                                                                                                      */
/*   X5          1.0                                                                                                      */
/*                                                                                                                        */
/*  Equation                                                                                                              */
/*                                                                                                                        */
/*  Y=2*X1+3*X2-1.5*X3+0.5*X4-X5                                                                                          */
/*                                                                                                                        */
/*  The Regression estimates are very different from the known values.                                                    */
/*                                                                                                                        */
/*                         Regression                                                                                     */
/*             Known       Coeficient                                                                                     */
/*  Variable  Coeficients   Estimates                                                                                     */
/*                                                                                                                        */
/*  X1             2         2.68105                                                                                      */
/*  X2             3        10.92363                                                                                      */
/*  X3           1.5       -35.13521                                                                                      */
/*  X4           0.5        26.47417                                                                                      */
/*  X5           1.0               0                                                                                      */
/*                                                                                                                        */
/*  Properties of Hilbert rows.                                   Banded structure and symetry of Hilbert rows.           */
/*                                                                                                                        */
/*    1. No row is a linear combination of other rows.             1 .50 .33 .25 .20 .                                    */
/*       So row relationships are independent and uncorelated,      .50 .33 .25 .20 .17                                   */
/*       even though regression may find a linear                      .33 .25 .20 .17 .14                                */
/*       dependendence do to the instability of row                        .25 .20 .17 .14 .13                             */
/*       relationships.                                                      .20 .17 .14 .13 .11                          */
/*                                                                                                                        */
/*    2. Hilbert matrices(rows) are invertible, positive                                                                  */
/*       definate and have a determinant                                                                                  */
/*                                                                                                                        */
/*  Now lets run a regression of Yhat on X1-X5                                                                            */
/*                                                                                                                        */
/*  proc reg data=have;                                                                                                   */
/*   model yhat = h1-h5 / noint;                                                                                          */
/*  run;quit;                                                                                                             */
/*                                                                                                                        */
/*                                                                                                                        */
/**************************************************************************************************************************/

/*   _     _ _ _               _                                      _
/ | | |__ (_) | |__   ___ _ __| |_   _ __ ___  __ _ _ __ ___  ___ ___(_) ___  _ __
| | | `_ \| | | `_ \ / _ \ `__| __| | `__/ _ \/ _` | `__/ _ \/ __/ __| |/ _ \| `_ \
| | | | | | | | |_) |  __/ |  | |_  | | |  __/ (_| | | |  __/\__ \__ \ | (_) | | | |
|_| |_| |_|_|_|_.__/ \___|_|   \__| |_|  \___|\__, |_|  \___||___/___/_|\___/|_| |_|
 _                    _                         |___/
(_)_ __  _ __  _   _| |_
| | `_ \| `_ \| | | | __|
| | | | | |_) | |_| | |_
|_|_| |_| .__/ \__,_|\__|
        |_|
*/

data
     sd1.hmatrix (keep=h:)
     sd1.beta   (keep=beta)
     sd1.yhat   (keep=yhat YXact Err);

  %let n=5;

  call streaminit(4321);

  array xs[&n,&n] x1-x%eval(&n*&n);
  array h[&n] h1-h&n;
  array betas[&n] b1-b5 (2 3 -1.5 + 0.5 -1);

  n=&n;

  /*----                                                      ----*/
  /*---- Create Hilbert rows                                  ----*/
  /*----                                                      ----*/

  do i=1 to n;
    beta=betas[i];
    output sd1.beta;
  end;

  /*----                                                      ----*/
  /*---- Create Hilbert rows                                  ----*/
  /*----                                                      ----*/

  do i=1 to n;
    do j=1 to n;
       xs[i, j] = 1 / (i + j - 1);
    end;
  end;

  /*----                                                      ----*/
  /*---- Output Hilbert rows                                  ----*/
  /*----                                                      ----*/

  do j=1 to n;
    do i=1 to n;
      h[i] = xs[i,j];
    end;
    output sd1.hmatrix;
  end;

  /*----                                                      ----*/
  /*---- Use arbitrary betas(coeficients) and the equation    ----*/
  /*---- YXact = 2*X1 + 3*X2 -1.5*X3 + 0.5*X4 -X              ----*/
  /*----                                                      ----*/

  do i=1 to n;
    yXact=0;
    do j=1 to n;
      yXact=yXact+betas[j]*xs[i,j];
    end;

  /*----                                                      ----*/
  /*---- Use arbitrary betas(coeficients) and the equation    ----*/
  /*---- YXact = 2*X1 + 3*X2 -1.5*X3 + 0.5*X4 -X              ----*/
  /*----                                                      ----*/
    err = rand('normal',0,.1);
    yhat=yXact + err;
    output sd1.yhat;
  end;

  format h: 4.2;

  stop;
run;quit;

data have ;
 merge sd1.hmatrix sd1.yhat sd1.beta;
run;quit;

/**************************************************************************************************************************/
/*                                                                                                                        */
/* WORK.HAVE total obs=5 03JUL2024:11:43:39                                                                               */
/*                                                                                                                        */
/* Obs       H1         H2         H3         H4         H5    BETA   YXACT        ERR        YHAT                        */
/*                                                                                                                        */
/*  1     1.00000    0.50000    0.33333    0.25000    0.20000   2.0  2.92500     0.12407    3.04907                       */
/*  2     0.50000    0.33333    0.25000    0.20000    0.16667   3.0  1.55833    -0.05353    1.50480                       */
/*  3     0.33333    0.25000    0.20000    0.16667    0.14286  -1.5  1.05714    -0.10139    0.95575                       */
/*  4     0.25000    0.20000    0.16667    0.14286    0.12500   0.5  0.79643     0.06896    0.86539                       */
/*  5     0.20000    0.16667    0.14286    0.12500    0.11111  -1.0  0.63710    -0.03246    0.60464                       */
/*                                                                                                                        */
/**************************************************************************************************************************/

/*
 _ __  _ __ ___   ___ ___  ___ ___
| `_ \| `__/ _ \ / __/ _ \/ __/ __|
| |_) | | | (_) | (_|  __/\__ \__ \
| .__/|_|  \___/ \___\___||___/___/
|_|
*/

proc reg data=have;
 model yhat = h1-h5 / noint;
run;quit;

/**************************************************************************************************************************/
/*                                                                                                                        */
/*                                                                                                                        */
/*  The Regression estimates are very different from the known values.                                                    */
/*                                                                                                                        */
/*                         Regression                                                                                     */
/*             Known       Coeficient                                                                                     */
/*  Variable  Coeficients   Estimates                                                                                     */
/*                                                                                                                        */
/*  X1             2         2.68105                                                                                      */
/*  X2             3        10.92363                                                                                      */
/*  X3           1.5       -35.13521                                                                                      */
/*  X4           0.5        26.47417                                                                                      */
/*  X5           1.0               0                                                                                      */
/*                                                                                                                        */
/*                                                                                                                        */
/* The REG Procedure                                                                                                      */
/* Model: MODEL1                                                                                                          */
/* Dependent Variable: YHAT                                                                                               */
/*                                                                                                                        */
/* Number of Observations Read           5                                                                                */
/* Number of Observations Used           5                                                                                */
/*                                                                                                                        */
/*                                                                                                                        */
/* NOTE: No intercept in model. R-Square is redefined.                                                                    */
/*                                                                                                                        */
/*                              Analysis of Variance                                                                      */
/*                                                                                                                        */
/*                                     Sum of           Mean                                                              */
/* Source                   DF        Squares         Square    F Value    Pr > F                                         */
/*                                                                                                                        */
/* Model                     4       13.57724        3.39431     284.02    0.0445                                         */
/* Error                     1        0.01195        0.01195                                                              */
/* Uncorrected Total         5       13.58919                                                                             */
/*                                                                                                                        */
/*                                                                                                                        */
/* Root MSE              0.10932    R-Square     0.9991                                                                   */
/* Dependent Mean        1.39593    Adj R-Sq     0.9956                                                                   */
/* Coeff Var             7.83142                                                                                          */
/*                                                                                                                        */
/* NOTE: Model is not full rank. Least-squares solutions for the parameters are not                                       */
/*       unique. Some statistics will be misleading. A reported DF of                                                     */
/*       0 or B means that the estimate is biased.                                                                        */
/*                                                                                                                        */
/* NOTE: The following parameters have been set to 0, since the variables are a                                           */
/*       linear combination of other variables as shown.                                                                  */
/*                                                                                                                        */
/*                                                                                                                        */
/* H5 =  -0.01638 * H1 + 0.30988 * H2 - 1.34466 * H3 + 2.03866 * H4                                                       */
/*                                                                                                                        */
/*                                                                                                                        */
/*                         Parameter Estimates                                                                            */
/*                                                                                                                        */
/*                      Parameter       Standard                                                                          */
/* Variable     DF       Estimate          Error    t Value    Pr > |t|                                                   */
/*                                                                                                                        */
/* H1            B        2.68105       20.22017       0.13      0.9161                                                   */
/* H2            B       10.92363      216.46719       0.05      0.9679                                                   */
/* H3            B      -35.13521      506.58567      -0.07      0.9559                                                   */
/* H4            B       26.47417      323.31905       0.08      0.9480                                                   */
/* H5            0              0              .        .         .                                                       */
/*                                                                                                                        */
/**************************************************************************************************************************/

/*___    _           _                           _
|___ \  (_)_ __   __| | ___ _ __   ___ _ __   __| | ___ _ __   ___ ___
  __) | | | `_ \ / _` |/ _ \ `_ \ / _ \ `_ \ / _` |/ _ \ `_ \ / __/ _ \
 / __/  | | | | | (_| |  __/ |_) |  __/ | | | (_| |  __/ | | | (_|  __/
|_____| |_|_| |_|\__,_|\___| .__/ \___|_| |_|\__,_|\___|_| |_|\___\___|
                           |_|
*/


The apparent contradiction between the independence of rows in a Hilbert matrix and the
linear dependence found in regression can be explained by the matrix's extreme
ill-conditioning. Here's a detailed explanation:

Row independence:
  The rows (and columns) of a Hilbert matrix are indeed linearly independent. This means
  that the matrix is invertible and has full rank. Mathematically, there is no exact
  linear combination of rows that can produce another row.

Ill-conditioning:
  Despite being invertible, the Hilbert matrix is notoriously ill-conditioned.
  The condition number of an n x n Hilbert matrix grows exponentially with n, approximately
  as O((1+v2)^(4n)/vn). This means that even for relatively small sizes, the matrix
  becomes extremely sensitive to small perturbations.
  Near linear dependence:
  Due to its ill-conditioning, the columns (or rows) of a Hilbert matrix are nearly
  linearly dependent. Each column is very close to being a multiple of the other columns,
  even though they are not exactly linearly dependent.

Numerical instability:
  In practical computations using finite-precision arithmetic (as in most computer systems),
  the near-linear dependence of the Hilbert matrix's columns can lead to numerical instability.
  Rounding errors and limitations in floating-point precision can make the matrix appear
  singular or very close to singular in numerical computations.

Regression and least squares:
  When using the Hilbert matrix in regression or least squares problems,
  the near-linear dependence of its columns can cause significant issues.
  The least squares solution becomes highly sensitive to small changes in
  the input data, and the coefficients can be wildly inaccurate due to the
  amplification of numerical errors.
Practical implications:
  In practice, regression algorithms may struggle to distinguish between the
  nearly linearly dependent columns, leading to results that suggest linear
  dependence even though the matrix is mathematically of full rank.

In summary, while the Hilbert matrix is theoretically invertible with
independent rows, its extreme ill-conditioning makes it behave as if it
were singular or had linearly dependent columns in numerical computations.
This is why regression algorithms, which are subject to numerical limitations,
may find apparent linear dependence despite the theoretical independence of the rows.

/*              _
  ___ _ __   __| |
 / _ \ `_ \ / _` |
|  __/ | | | (_| |
 \___|_| |_|\__,_|

*/

