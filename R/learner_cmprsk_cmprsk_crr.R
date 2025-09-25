#' @title Fine-Gray Competing Risks Regression Learner
#'
#' @name mlr_learners_cmprsk.crr
#'
#' @description
#' A learner for Fine-Gray competing risks regression to estimate cumulative incidence functions (CIFs)
#' for multiple mutually exclusive events. Calls [cmprsk::crr()] from \CRANpkg{cmprsk} for training and [cmprsk::predict.crr()]
#' for prediction. Supports fixed and time-varying covariates.
#'
#' @template learner
#' @templateVar id cmprsk.crr
#'
#' @section Parameter Details:
#' \describe{
#'   \item{maxiter}{(`integer(1)`)\cr Maximum number of iterations for convergence.}
#'   \item{gtol}{(`numeric(1)`)\cr Convergence tolerance for gradient.}
# `   \item{na.action){`logical`)\cr ??? .)
#'   \item{init_list}{(`list()`)\cr Initial values for model coefficients, named by event type.}
#'   \item{censor_group}{(`character(0)`)\cr Censoring group column name, if any.}
#'   \item{cov2_info}{(`list()`)\cr Information for time-varying covariates, including `feat2` (names) and `tf` (transformation function).}
#' }
#'
#' @section Custom mlr3 parameters:
#' * `cov2_info`: Configuration for time-varying features (default: `NULL` for model with fixed covariates).
#' * `init_list`: Initial regression parameters per cause (default: `NULL` for zeros).
#' * `censor_group`: Column name for groups with different censoring distributions (default: `character(0)`).
#' @param cov2_info 'list()'\cr
#' Optional configuration for time-varying covariates, enabling the learner to
#' model covariate effects that change over time. This list must contain:
#'
#'
#' \describe{
#'   \item{feat2}{'character()'\cr
#'     A vector of feature names from the task's feature set that are
#'     treated as time-varying. }
#'   \item{tf}{'function(uft)'\cr
#'     A user-defined function specifying how the features in 'feat2' vary
#'     over time. It takes one argument: 'uft' (a numeric vector of unique
#'     failure times from the training data). The function's behavior is
#'     described in the cmprsk::crr() documentation and must return a matrix
#'     with:\cr
#'     - 'nrow = length(uft)' (matching the number of unique failure times).\cr
#'     - 'ncol' equal to the number of columns in the 'cov2' matrix (derived
#'       from 'feat2' vector via binding columnwise `model.matrix'), where each column corresponds to a
#'       time-varying effect for each column in 'cov2' matrix.\cr
#'     Example: 'function(uft) log(uft)' applies a logarithmic transformation
#'     to the time points.}
#'   \item{feat2only}{'character()' or 'NULL'\cr
#'     A vector of feature names that are used solely to build the
#'     time-varying covariate matrix ('cov2') and are excluded from the fixed
#'     covariate matrix ('cov1'). Must be a subset of 'feat2'. If 'NULL'
#'     (default), all features in the task contribute to 'cov1', and 'feat2'
#'     defines 'cov2' matrix.}
#' }
#' If 'cov2_info' is 'NULL' (default), the learner treats all covariates as
#' fixed.
#'
#'
#' @section Initial parameter values:
#' - `maxiter`: Set to 100 (package default is 10).
#' - `gtol`: Set to 1e-6 (package default is 1e-3).
#' - `na.action`: Set to na.omit.
#' - `init_list`: Set to NULL (no corresponding parameter in cmprsk
#' - `censor_group`: Set to NULL,
#' - `cov2_info`: Set to NULL

#'
#' @references
#' `r mlr3misc::format_bib("finegray1999crr")`
#'
#' @template seealso_learner
#' @examplesIf learner_is_runnable("cmprsk.crr")
#' library(mlr3)
#' library(mlr3proba)
#'
#' # Define the learner
#' learner <- lrn("cmprsk.crr")
#' print(learner)
#'
#' # Define a task
#' task <- tsk("pbc")
#' task$select(c("age", "bili", "sex"))
#' task$set_col_roles(cols = "status", add_to = "stratum")
#'
#' # Create train and test sets
#' ids <- partition(task)
#'
#' # Train the learner
#' learner$train(task, row_ids = ids$train)
#'
#' # Print model, convergence, and importance
#' print(learner$model)
#' print(learner$convergence())
#'
#' # Make predictions
#' predictions <- learner$predict(task, ids$test)
#'
#' # Score the predictions
#' predictions$score()
#' @export
LearnerCompRisksFineGrayCRR <- R6::R6Class(

  "LearnerCompRisksFineGrayCRR",
  inherit = mlr3proba::LearnerCompRisks,
  public = list(

    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    initialize = function() {
      param_set <- ps(
        maxiter = p_int(lower = 1L, upper = 1000L, default = 10L, tags = c("train")),
        gtol = p_dbl(lower = 1e-9, upper = 1e-3, default = 1e-6, tags = c("train")),
        variance = p_lgl(default = TRUE, tags = c("train")),
        na.action = p_uty(default = stats::na.omit, tags = "train"),
        init_list = p_uty(default = NULL, tags = c("train")),
        censor_group = p_uty(default = NULL, tags = c("train")),
        cov2_info = p_uty(default = NULL, tags = c("train", "predict"))
      )

      param_set$values <- list(
        variance = TRUE,
        na.action = stats::na.omit
       )

      super$initialize(
        id = "cmprsk.crr",
        predict_types = "cif",
        feature_types = c("logical", "integer", "numeric", "factor"),
        packages = c("mlr3extralearners", "cmprsk"),
        param_set = param_set,
        properties = c("missings"), # "convergence",
        man = "mlr3extralearners::mlr_learners_cmprsk.crr",
        label = "Competing Risks Regression: Fine-Gray model"
      )
    }
  ) # end of public list
  ,
  private = list(
    .create_cov = function(data, feat) {
      cmtx_fun <- function(data, cx) {
        # cx is a covariate name,
        # logical, numeric, factor, integer vars accepted.
        cff <- paste0("~", cx)
        ff <- as.formula(cff)
        m <- model.matrix(ff, data = data)[, -1]
        if (is.vector(m, mode = "numeric")) {
          m <- matrix(m, ncol = 1)
          colnames(m) <- cx
        }
        m # Returns model matrix (wout intercept) for one covariate
      }

      logger <- lgr::get_logger("mlr3")
      func <- "[cmprsl.crr] private$.create_cov "
      logger$debug("%s STARTS. ", func)
      res <- NULL
      assertDataFrame(data, min.rows = 1, min.cols = 1)
      assertCharacter(feat, any.missing = FALSE)
      assertSubset(unique(feat), colnames(data))

      if (length(feat)) {
        cmtx_list <- lapply(feat, function(cx) cmtx_fun(data, cx))
        # Bind matrices by columns
        res <- do.call(cbind, cmtx_list)
      }
      logger$debug("%s ENDS. ", func)
      res
    },

    .create_xcov = function(task, cov2_info) {
      logger <- lgr::get_logger("mlr3")
      func <- "[cmprsk.crr] private$.create_xcov: "
      logger$debug("%s STARTS", func)

      # Get data and feature names
      data <- task$data()
      feature_names <- task$feature_names
      uft <- task$unique_event_times()

      # Check cov2_info validity
      check_cov2_info <- private$.check_cov2_info(cov2_info, feature_names, uft)
      if (!check_cov2_info) {
        err <- sprintf("%s Invalid cov2_info parameter", func)
        logger$error(err)
        stop(err)
      }

      # Determine feat1 and feat2
      if (!is.null(cov2_info)) {
        feat2 <- cov2_info$feat2
        cov2only <- cov2_info$cov2only
        feat1 <- setdiff(unique(c(feature_names, feat2)), unique(cov2only))
      } else {
        feat1 <- feature_names
        feat2 <- NULL
      }

      # Create cov1 and cov2 matrices
      xcov_list <- vector("list", length = 2)
      names(xcov_list) <- c("cov1", "cov2")

      xcov_list$cov1 <- if (!is.null(feat1)) private$.create_cov(data, feat1) else NULL
      xcov_list$cov2 <- if (!is.null(feat2)) private$.create_cov(data, feat2) else NULL
      logger$debug("%s ENDS", func)

      return(xcov_list)
    },
    
    .check_cov2_info = function(cov2_info, feature_names, uft) {
      logger <- lgr::get_logger("mlr3")
      func <- "[cmprsk:crr] private$.check_cov2_info"

      chk_err <- function(chk) {
        if (is.character(chk)) {
          err <- sprintf("%s %s", func, chk)
          logger$error(err)
          stop(err)
        }
        return()
      }

      # Check if cov2_info is NULL (valid for fixed covariates)
      if (is.null(cov2_info)) {
        return(TRUE)
      }
      # Additional checks go here
      if (!typeof(cov2_info) == "list") {
        logger$error("%s cov2_info is not a list", func)
        stop("ERROR: cov2_info is not a list")
      }

      t1 <- c("feat2", "tf")
      t2 <- names(cov2_info)
      tt <- intersect(t1, t2)
      if (!identical(tt, t1)) {
        err <- sprintf("%s cov2_info should contain at least two components: feat2, and uft", func)
        logger$error(err)
        stop(err)
      }

      chk <- check_character(cov2_info$feat2, min.len = 1)
      chk_err(chk)

      chk <- check_subset(cov2_info$feat2, feature_names, empty.ok = FALSE)
      chk_err(chk)
      return(TRUE)
    },
    
  .add_to_args = function(args = list(), x = NULL) {
    if (is.null(x)) {
      return(args)
    }
      logger <- lgr::get_logger("mlr3")
      func <- "[cmprsk.crr] private$.train: "
     
       logger$debug("%s Execution STARTS", func)
  
      print("----add_to_args1 ---")
       print(typeof(x))
        print(str(x))
        print(names(x))

    olist =  if (typeof(x) != "list") setNames(list(x), deparse(substitute(x))) else xlist
     print("----add_to_args2 ---")
        print(names(olist))
       print(str(olist))
      print("----add_to_arg3 ---")
      
    mlr3misc::insert_named(args, olist)
  },
  
    .train = function(task) {
      logger <- lgr::get_logger("mlr3")
      func <- "[cmprsk.crr] private$.train: "

      logger$debug("%s STARTS", func)
      # Parameter Extraction

      pv <- self$param_set$get_values(tags = "train")

      message("pv1---msg")
      print(pv)
      
      target_names <- task$target_names
      target_df = task$data(cols = target_names)
      
      cov2_info = pv$cov2_info

      # Create list with cov1 and/or cov2
      args = list()
      xcov_args = private$.create_xcov(task, cov2_info)
      #message("===xcov_args")
      #print(names(xcov_args))
      #print(str(xcov_args))
      #cov1 = xcov_args$cov1
      #names(cov1) = "cov1"
      #print("==== cov1 ===")

      #print(str(cov1))
      args = private$.add_to_args(args, xcov_args) 
      #cov2 = xcov_args$cov2
      #print("==== cov2 ===")
      #print(str(cov2))
      #args = private$.add_to_args(args, cov2) 
      tf = cov2_info$tf
      args = private$.add_to_args(args, tf) 

      if (!is.null(pv$censor_group)) {
        logger$debug("%s censor_group name in backend data is %s", func, pv$censor_group)
        cengroup = task$backend$data(rows = task$row_ids, cols = pv$censor_group)[[1]]
      } else {
        cengroup = NULL
      }
      args = private$.add_to_args(args, cengroup) 

      subset  = NULL  #  a logical vector extracted subset
      args = private$.add_to_args(args, subset) 

      print("----args before the loop")
      print(names(args))
      print(str(args))

      unique_events = task$unique_events()
      model = lapply(seq_along(unique_events), function(i) {
        uei = unique_events[i]
        logger$debug("%s Training for cause = %s", func, uei)
        init = if (!is.null(pv$init_list)) pv$init_list[[uei]] else NULL
        args = private$.add_to_args(args, cov2) 
        print("----args inside the loop")
        print(names(args))
        print(str(args))

        rlang::exec(cmprsk::crr,
          ftime = target_df[[1]],
	  fstatus = target_df[[2]],
	  cencode = 0L,
	  failcode = uei,
	  na.action = pv$na.action,
	  maxiter <- pv$maxiter %??% 10L,
	  gtol = pv$gtol %??% 1e-6,
	  variance = pv$variance,
	  !!!args)
      })


      logger$debug("%s Training completed, model has %d components", func, length(model))
      return(model)
    },
    .predict = function(task) {
      logger <- lgr::get_logger("mlr3")
      func <- "[cmprsk.crr] predict$.predict "
      logger$debug("%s Starts ", func)
      pv <- self$param_set$get_values(tags = "predict")
      cmp_events <- task$cmp_events
      cov2_info <- pv$cov2_info
      xcov_args <- private$.create_xcov(task, cov2_info)
      xcov_nms <- names(xcov_args)

      uftimes <- task$unique_event_times()

      cif_list <- lapply(
        seq_along(cmp_events),
        function(cmp_event) {
          logger$debug("%s -- cmp_event/cause = %d", func, cmp_event)
          fit <- self$model[[cmp_event]]
          obj <- list(object = fit)
          uft <- fit$uftime
          args <- insert_named(obj, xcov_args)
          pred <- rlang::exec(cmprsk::predict.crr, !!!args)
          cif <- t(pred[, -1])
          colnames(cif) <- as.character(uft)
          cif
        }
      )
      names(cif_list) <- cmp_events

      return(list(cif = cif_list))
    }
  ) # end of private list
) # end Rclass definition

.extralrns_dict$add("cmprsk.crr", LearnerCompRisksFineGrayCRR)
