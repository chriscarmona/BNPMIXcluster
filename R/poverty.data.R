#'
#' @docType data
#'
#' @name
#'     poverty.data
#'
#' @title
#'     Poverty data for testing the \emph{BNPMIXcluster} package
#'
#' @description
#'
#'     Poverty indicators observed in Mexico for 2014.
#'
#'     The original data is available in the file "\code{R_2014.zip}" from CONEVAL's website: \url{http://www.coneval.org.mx/Medicion/MP/Paginas/Programas_BD_10_12_14.aspx}
#'
#'     (download zip file directly from:
#'     \url{http://www.coneval.org.mx/Medicion/MP/Documents/Programas_calculo_pobreza_10_12_14/R_2014.zip}
#'     )
#'
#'     This data frame presents indicators aggregated by household. The aggregation was done by the authors according with code in section \code{Examples}.
#'
#' @usage
#'     poverty.data
#'
#' @details
#'
#'     \code{poverty.data} is a data frame with 58121 rows and 13 variables, with the following columns:
#' \describe{
#'     \item{\code{proyecto}}{ Data source identifier (1=MCS, 2=ENIGH) }
#'     \item{\code{folioviv}}{ Household identifier level 1 }
#'     \item{\code{foliohog}}{ Household identifier level 2 }
#'     \item{\code{ict_norm}}{(continuous) Total income in the household (in Mexican pesos).}
#'     \item{\code{ic_ali}}{(binary) Indicator for deprivation to feeding: 1-yes,0-no }
#'     \item{\code{ic_asalud}}{(binary) Indicator for deprivation to health services: 1-yes,0-no }
#'     \item{\code{ic_cv}}{(binary) Housing quality: 1-bad, 0-good }
#'     \item{\code{ic_rezedu}}{(binary) Indicator for education backwardness: 1-yes,0-no }
#'     \item{\code{ic_sbv}}{(binary) Indicator for deprivation to basic public services: 1-yes,0-no }
#'     \item{\code{ic_segsoc}}{(binary) Indicator for deprivation to social security: 1-yes,0-no }
#'     \item{\code{niv_ed}}{(categorical, ordered) Maximum educational level in the household: 0-incomplete primary: 1-incomplete secondary, 2-complete secondary or more }
#'     \item{\code{tam_loc}}{(categorical, nominal) Size of locality according to the number of people living there: 1-(100000, \eqn{\infty}), 2-[15000, 100000), 3-[2500, 15000), 4-[0, 2500)}
#'     \item{\code{factor_hog}}{Expansion factor for the household, according to the complex survey design.}
#' }
#'
#' @examples
#'
#' ##### Generates poverty.data using the original data from CONEVAL's website #####
#'
#' \dontrun{
#' # step 1:
#' #     Download and unzip the file "R_2014.zip"
#' #     available in:
#' #     http://www.coneval.org.mx/Medicion/MP/Documents/Programas_calculo_pobreza_10_12_14/R_2014.zip
#'
#' # step 2:
#' #     extract and read the csv file "pobreza_14.csv"
#'
#' coneval.poverty.data <- read.csv("pobreza_14.csv", na.strings=c("NA",""))
#'
#' # step 3:
#' #     Execute the following code...
#'
#' var_id <- c("proyecto","folioviv","foliohog","numren")
#' for(i in match(var_id,colnames(coneval.poverty.data)) ){
#'      coneval.poverty.data[,i] <- formatC( x=as.numeric(coneval.poverty.data[,i]),
#'                                           width=max(nchar(coneval.poverty.data[,i])),
#'                                           format="f",flag="0",digits=0
#'                                           )
#' }
#'
#' # normalizing the continuous variable for income #
#' b <- quantile(coneval.poverty.data$ict,probs=0.01)
#' coneval.poverty.data$ict_norm <- log(coneval.poverty.data$ict+b)
#'
#' # Aggregating data at household level
#' Y_names <- c("ict_norm",
#'              "ic_ali","ic_asalud","ic_cv",
#'              "ic_rezedu","ic_sbv","ic_segsoc",
#'              "niv_ed","tam_loc")
#' agg_form <- as.formula( paste( "cbind(",paste(c(Y_names,"factor_hog"),collapse=",") ,")",
#'                                "~proyecto+folioviv+foliohog"
#'                              )
#'                       )
#'
#' poverty.data <- aggregate(agg_form,FUN="max",data=coneval.poverty.data)
#'
#' }
#'
#' @seealso
#' \code{\link{MIXclustering}}
#'
NULL
