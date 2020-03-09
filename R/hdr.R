#' Hubbard Decision Research Pseudo-Random Number Generator
#'
#' Hubbard Decision Research (HDR) Pseudo-Random Number Generator (PRNG)
#' @param x a vector (or a `n * m` matrix) of seeds, where `m`` corresponds to dimensions of random numbers.
#' Default `m<=5`. See **HDR Dimensions** section below.
#' @param t1,t2 T constants (prime numbers) for 1st and 2nd term, respectively.
#' The length of these vectors determine maximum number of dimensions for HDR PRNG.
#' Default values are t1=c(2499997, 1800451, 2000371, 1796777, 2299603)
#'                and t2=c(2246527, 2399993, 2100869, 1918303, 1624729)
#' @param a1,a2 A constants. Default values are a1=7450589, a2=7450987
#' @param b1,b2 B constants. Default values are b1=4658, b2=7580
#' @param c1,c2 C constants. Default values are c1=7450581, c2=7560584
#' @param d1,d2 D constants. Default values are d1=383, d2=17669
#' @param e1,e2 E constants. Default values are e1=99991, e2=7440893
#' @param f1,f2 F constants. Default values are f1=7440893, f2=1343
#'
#' @details
#'
#' HDR PRNG is given by the formula:
#'
#' \ifelse{html}{\out{R(X)=mod(mod(mod(10<sup>15</sup>-11,mod(X x T<sub>1</sub>,A<sub>1</sub>)*B<sub>1</sub>+C<sub>1</sub>)*D<sub>1</sub>,E<sub>1</sub>)*F<sub>1</sub>+
#' mod(mod(10<sup>15</sup>-11, mod(X x T<sub>2</sub>,A<sub>2</sub>)*B<sub>2</sub>+C<sub>2</sub>)*D<sub>2</sub>,E<sub>2</sub>)*F<sub>2</sub>,2<sup>32</sup>)}}{\deqn{R(x)=mod(mod(mod(10^15-11,mod(x*T1,A1)*B1+C1)*D1,E1)*F1+mod(mod(10^15-11, mod(x*T2,A2)*B2+C2)*D2,E2)*F2,2^32)}}
#'
#' Further details on each of the dimensions
#'
#' | Term | Dimension | Description |
#' | ---- | --------- | ----------- |
#' | 1 | Trial ID | This represents a unique identifier for a given scenario in a simulation. This 8 decimal digit identifier allows for up to 100 million unique trials for each variable in a model|
#' | 2 | Variable ID | This is a unique identifier for a variable. It would be an 8-digit variable ID allowing for up to 100 million unique variables. For example, if “Monthly Demand for Product X” and “average time spent in activity Y” were variables in a model, they would each be given unique variable IDs. Organizations may structure their Variables IDs so that related variables are in groups. For example, perhaps all marketing and sales related variables have “11” for the first two digits and all cybersecurity related variables have “73” for the first two digits, and so on. Variable IDs could be assigned on an ad hoc basis but a large organization making many models with a lot of shared variables would  want  to  develop  an  internal  library  of  assigned  variable  IDs  similar  to  an accountant’s “chart of accounts.”|
#' | 3 | Entity ID | This  identifies an organization  or some other category  of  users.  A corporation or government agency may be assigned a unique 8 decimal digit Entity ID. Since this provides for 100 million potential entities, that should be enough for every business, not for profit and government agency that wants one on the planet. This is useful if there  are  models  using  random  variables  from  many  organizations  do  not  have variables that produce the same random sequences. For example, many banks may use variables defined by the FDIC for “stress testing” to ensure banks are financially stable even during times of economic stress. The bank would want to ensure that internally defined variables with the same Variable ID are not correlated to the shared variables. The FDIC would supply the variable ID along with the Entity ID of the FDIC so that every  bank  using  those  variables  produces  the  same  sequence  while  avoiding duplicating the sequence of internally defined variables. A default Entity ID of 0 can be used by anyone as long as sharing variables would not be an issue.|
#' | 4 | Time ID | This identifies a particular time unit for a given variable/trial/entity combination. This allows one scenario for a given variable to contain an entire unique time series. A 7-digit time series ID would allow for time series containing 115 days of seconds, 19 years of minutes, or 27,397 years of days. This is an optional dimension. Variables that do not represent a time series use the default Time ID of 0.|
#' | 5 | Agent ID | This provides a fifth optional dimension for the counter-based PRNG. One possible use is as an identify for agents in agent-based modeling. If this ID is not used, the default value is 0.|
#'
#' @return vector or pseudo-random numbers related for every one of (combination of) seeds.
#' @references D. W. Hubbard, "A Multi-Dimensional, Counter-Based Pseudo Random Number Generator
#' as a Standard for Monte Carlo Simulations," 2019 Winter Simulation Conference (WSC),
#' National Harbor, MD, USA, 2019, pp. 3064-3073. DOI: 10.1109/WSC40007.2019.9004773
#'
#' @examples
#' rHDR(c(1:10))
#' rHDR(matrix(c(1:10), byrow=TRUE, nrow=5))
#' @export
rHDR <- function(x, t1=c(2499997, 1800451, 2000371, 1796777, 2299603), a1=7450589, b1=4658, c1=7450581, d1=383, e1=99991, f1=7440893,
                    t2=c(2246527, 2399993, 2100869, 1918303, 1624729), a2=7450987, b2=7580, c2=7560584, d2=17669, e2=7440893, f2=1343){
   if(is.null(dim(x)[2])){
    dim(x) <- c(length(x),1)
   }
  # limit it to dimensions supplied in t1 and t2 vectors
  w <- seq.int(min(dim(x)[2], length(t1)))
  s1 <- as.vector(x %*% t1[w])
  s2 <- as.vector(x %*% t2[w])

  rx <- `%%`(
    (`%%`(
      `%%`(1e15-11L,
           `%%`(s1,a1)*b1 + c1)*d1, e1)*f1 +
       `%%`(
         `%%`(1e15-11L,
              `%%`(s2,a2)*b2 + c2)*d2, e2))*f2, 2^32)

  (rx +.5)/2^32

}

