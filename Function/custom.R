#-------------------------------------------------------------------------------------------------
# custom_meta_aggregate was inspired by https://github.com/meghapsimatrix/shiny_egm/blob/main/tidy_meta.R
# but the computational of mean effects is changed. We used robust point estimation to estimate mean effects. 
# Details see Yang Y, Lagisz M, Williams C, et al. Robust point and variance estimation for ecological and evolutionary meta-analyses with selective reporting and dependent effect sizes[J]. 2023. https://ecoevorxiv.org/repository/view/6018/ 
#-------------------------------------------------------------------------------------------------

#' @title custom_meta_aggregate
#' @description custom_meta_aggregate is used to aggregate effect sizes
#'  and their variances from multiple studies.
#'  The function can also handle the situation when there is only one effect size in the dataset.
#'  In this case, the function will return the effect size and its variance.
#' @param data A data frame containing effect sizes and their variances.
#'  The data frame should contain at least three columns: study_id, yi, and vi.
#'  study_id is the unique identifier for each study. 
#'  yi is the effect size and vi is the variance of the effect size.  
#' @param rho is the intra-class correlation coefficient (ICC) for the random effects model.
#'  The default value of rho is 0.5.
#' @return A data frame containing the aggregated effect size, the number of studies, and the number of effect sizes.

custom_meta_aggregate <- function(data, rho = 0.5){
  
  n_studies <- length(unique(data$study_id))
  n_es <- nrow(data)
  n_info <- tibble(n_studies = n_studies, n_es = n_es)
  data$study_id <- as.factor(data$study_id)
  
  if(n_studies >= 2){
    data$es_id <- 1:nrow(data)
    data$es_id <- as.factor(data$es_id)
    VCV <- vcalc(vi = vi,
                 cluster = study_id, 
                 rho = rho, 
                 obs = es_id, 
                 data = data) 
    
    suppressWarnings(mod <- rma.mv(yi = yi, 
                                   V = VCV, 
                                   method = "REML", 
                                   test = "t", 
                                   dfs = "contain", 
                                   data = data)) # robust point estimation
    
    # Equal-effect model or weighted least square also can be used to obtain robust point estimate
    #suppressWarnings(mod <- rma(yi = yi, 
    #                            vi = vi, 
    #                            method = "EE",
    #                            data = data)) 
    
    #suppressWarnings(mod <- lm(yi ~ 1, weights = 1/vi, data=data))
    mod_rob <- robust(mod, 
                      cluster = study_id, 
                      adjust = TRUE, 
                      clubSandwich = TRUE) 
    
    res <- tibble(estimate = mod_rob$beta[[1]],
                  SE = mod_rob$se,
                  CI_L = mod_rob$ci.lb,
                  CI_U = mod_rob$ci.ub)

    
  } else if(n_es == 1){
    
    res <- tibble(estimate = data$yi,
                  SE = sqrt(data$vi), 
                  CI_L = data$yi - qnorm(1-0.05/2,lower.tail = T) * sqrt(data$vi), 
                  CI_U = data$yi + qnorm(1-0.05/2,lower.tail = T) * sqrt(data$vi))
    
  } else{
    
    suppressWarnings(mod <- rma(yi = yi,
                                vi = vi,
                                method = "EE",
                                data = data)) # equal-effect model
    
    res <- tibble(estimate = mod$beta[[1]],
                  SE = mod$se,
                  CI_L = mod$ci.lb,
                  CI_U = mod$ci.ub)
    
  } 
  
  
  summary_data <- bind_cols(res, n_info) %>%
    mutate_if(is.numeric, round, 3) %>%
    select(estimate, n_studies, n_es)
  
  return(summary_data)
}



#-------------------------------------------------------------------------------------------------
# custom_altmetric_aggregate is used to aggregate altmetric scores
#-------------------------------------------------------------------------------------------------

#' @title custom_altmetric_aggregate
#' @description custom_altmetric_aggregate is used to aggregate altmetric scores.
#' @param  data A data frame containing altmetric scores. 
#' I needs to have study_id column with study IDs and score column with altmetric scores.
#' @return A data frame containing the aggregated altmetric score, the number of studies, and the number of effect sizes.

custom_altmetric_aggregate <- function(data){
  
  n_studies <- length(unique(data$study_id))
  n_es <- nrow(data)
  n_info <- tibble(n_studies = n_studies, n_es = n_es)
  data$study_id <- as.factor(data$study_id)
  
  if(n_studies >= 2){
    suppressMessages(suppressWarnings(mod <- lme4::glmer(score ~ 1 + (1 | study_id), 
                                        data = data, 
                                        family = poisson(link="log"))))
    
    res <- suppressMessages(tibble(estimate = exp(summary(mod)$coefficients[1]), # medium
                  # mean - exp((summary(mod)$coefficients[1] + 
                  # 0.5*(summary(mod)$varcor[[1]][[1]] + # sigma^2 for study level
                  #       summary(mod)$sigma^2) # residual level - we cannot do like what we did above as we added 0.025
                  # ))
                  SE = exp(summary(mod)$coefficients[2]),
                  CI_L = exp(confint(mod)[2,1]),
                  CI_U = exp(confint(mod)[2,2])))
    
  } else if(n_es == 1){
    
    res <- tibble(estimate = data$score,
                  SE = NA, 
                  df = NA, 
                  CI_L = NA, 
                  CI_U = NA)
    
  } else{
    
    suppressWarnings(mod <- glm(score ~ 1, 
                                data = data, 
                                family = poisson(link="log")))
    
    res <- suppressMessages(suppressMessages(tibble(estimate = exp(summary(mod)$coefficients[1]), # medium
                                   SE = exp(summary(mod)$coefficients[2]),
                                   CI_L = exp(confint(mod)[1]),
                                   CI_U = exp(confint(mod)[2]))))
    
  } 
  
  
  summary_data <- bind_cols(res, n_info) %>%
    mutate_if(is.numeric, round, 3) %>%
    select(estimate, n_studies, n_es)
  
  return(summary_data)
}


#-------------------------------------------------------------------------------------------------
# custom_translation_aggregate is used to aggregate policy and patent citation counts
#-------------------------------------------------------------------------------------------------

#' @title custom_translation_aggregate
#' @description custom_translation_aggregate is used to aggregate policy and patent citation counts.
#' @param  data A data frame containing policy and patent citation counts.
#' I needs to have study_id column with study IDs and score column with altmetric scores.
#' @return A data frame containing the aggregated policy and patent citation counts,
#'   the number of studies, and the number of effect sizes.

custom_translation_aggregate <- function(data){
  
  n_studies <- length(unique(data$study_id))
  n_es <- nrow(data)
  n_info <- tibble(n_studies = n_studies, n_es = n_es)
  data$study_id <- as.factor(data$study_id)
  
  if(n_es == 1){
    res <- tibble(estimate = data$count,
                  SE = NA, 
                  df = NA, 
                  CI_L = NA, 
                  CI_U = NA)
    
  }  else{
    
    #suppressWarnings(mod <- glm(count ~ 1, data = data, family = poisson(link="log")))
    
    #res <- suppressMessages(suppressMessages(tibble(estimate = exp(summary(mod)$coefficients[1]), # medium
    #                                                SE = exp(summary(mod)$coefficients[2]),
    #                                                CI_L = exp(confint(mod)[1]),
    #                                                CI_U = exp(confint(mod)[2]))))
    res <- tibble(estimate = sum(data$count/n_es),
                  SE = NA, 
                  df = NA, 
                  CI_L = NA, 
                  CI_U = NA)
  } 
  
  
  summary_data <- bind_cols(res, n_info) %>%
    mutate_if(is.numeric, round, 3) %>%
    select(estimate, n_studies, n_es)
  
  return(summary_data)
}



#-------------------------------------------------------------------------------------------------
# # function custom_getAltmetrics()
#-------------------------------------------------------------------------------------------------

custom_getAltmetrics <- function(doi = NULL,
                          foptions = list(),
                          ...) {
  if (!is.null(doi)) doi <- stringr::str_c("doi/", doi)
  identifiers <- purrr::compact(list(doi))
  if (!is.null(identifiers)) {
    ids <- identifiers[[1]]
  }
  base_url <- "http://api.altmetric.com/v1/"
  #request <- httr::GET(paste0(base_url, ids), httr::add_headers("user-agent" = "#rstats rAltmertic package https://github.com/ropensci/rAltmetric"))
  request <- httr::GET(paste0(base_url, ids))
  results <-
    jsonlite::fromJSON(httr::content(request, as = "text"), flatten = TRUE)
  results <- rlist::list.flatten(results)
  class(results) <- "altmetric"
  return(results)
}

# helper function for Altmetric analysis
# format altmetric object
format.Altmetric <- function(altmetric.object) {
  stats <- altmetric.object[grep("^cited", names(altmetric.object))]
  stats <- data.frame(stats, stringsAsFactors = FALSE)
  data.frame(paper_title = altmetric.object$title,
             journal = altmetric.object$journal,
             doi = altmetric.object$doi,
             #subject = altmetric.object$subjects,
             Altmetric.score = altmetric.object$score,
             stats = stats)
}

# create a dataframe function
altmetric_df <- function(altmetric.object) {
  df <- data.frame(t(unlist(altmetric.object)), stringsAsFactors = FALSE)
}
#altmetric.crawler[[n]]  <-  try(list(altmetric_df(custom_getAltmetrics(doi = DOIs[n]))))
# create a function to summarize Altmetric object
summary.altmetric <- function(x, ...) {
  if (inherits(x, "altmetric"))  {
    string <- "Altmetrics on: \"%s\" with altmetric_id: %s published in %s."
    vals   <- c(x$title,  x$altmetric_id, x$journal)
    if("journal" %in% names(x)) {
      cat(do.call(sprintf, as.list(c(string, vals))))
    } else {
      string <- "Altmetrics on: \"%s\" with altmetric_id: %s"
      cat(do.call(sprintf, as.list(c(string, vals))))
    }
    cat("\n")
    stats <- x[grep("^cited", names(x))]
    stats <- data.frame(stats, stringsAsFactors = FALSE)
    print(data.frame(stats = t(stats)))
  }
}

altmetric_summary <- function(object) {
  altmetric.summary <- data.frame(paper = sapply(object, function(x)  ifelse(class(x) == "data.frame",x$paper_title,NA)),
                                  Journal = sapply(object, function(x) ifelse(class(x) == "data.frame",x$journal,NA)),
                                  DOI = sapply(object, function(x) ifelse(class(x) == "data.frame",x$doi,NA)),
                                  #subject = sapply(altmetric.crawler2, function(x) ifelse(class(x) == "data.frame",x$subject,NA)),
                                  `Altmetric score` = sapply(object, function(x) ifelse(class(x) == "data.frame",x$Altmetric.score,0)),
                                  Policy = sapply(object, function(x) ifelse(class(x) == "data.frame",ifelse(!is.null(x$stats.cited_by_policies_count),x$stats.cited_by_policies_count,0),0)),
                                  Patent = sapply(object, function(x) ifelse(class(x) == "data.frame",ifelse(!is.null(x$stats.cited_by_patents_count),x$stats.cited_by_patents_count,0),0))
  )
  return(altmetric.summary)
}





#-------------------------------------------------------------------------------------------------
# geom_sankey comes from ihttps://github.com/davidsjoberg/ggsankey/blob/main/R/sankey.R
#-------------------------------------------------------------------------------------------------
utils::globalVariables(c(".", ".data", "x", "node", "next_node", "next_x", "..r"))
# importFrom(ggplot2, "%+replace%")
#' @importFrom ggplot2 %+replace%

# ** Support functions ----------
prepare_params <- function(...) {
  # Prepare aesthics for flow lines
  flow.aes <- list(...)
  removes <- names(flow.aes) %>%
    stringr::str_extract_all(., "(?<=flow.).*") %>% unlist()
  removes2 <- names(flow.aes) %>%
    stringr::str_subset(., "node") %>% unlist()
  flow.aes[c(removes, removes2)] <- NULL
  names(flow.aes) <- names(flow.aes) %>%
    stringr::str_replace_all("flow.", "")
  
  # Prepare aesthics for node boxes
  node.aes <- list(...)
  removes <- names(node.aes) %>%
    stringr::str_extract_all(., "(?<=node.).*") %>% unlist()
  removes2 <- names(node.aes) %>%
    stringr::str_subset(., "flow") %>% unlist()
  node.aes[c(removes, removes2)] <- NULL
  names(node.aes) <- names(node.aes) %>%
    stringr::str_replace_all(., "node.", "")
  
  return(list(flow.aes, node.aes))
}

find_default_space <- function(.df) {
  .df %>%
    dplyr::group_by(.data$n_x) %>%
    dplyr::summarise(n_groups = dplyr::n_distinct(.data$node),
                     freq = sum(.data$freq, na.rm = TRUE)) %>%
    dplyr::mutate(v = .data$freq / .data$n_groups / 4) %>%
    dplyr::pull(.data$v) %>%
    max()
}

sigmoid <- function(x_from, x_to, y_from, y_to, smooth = 5, n = 300) {
  x <- seq(-smooth, smooth, length = n)
  y <- exp(x) / (exp(x) + 1)
  out <- data.frame(x = (x + smooth) / (smooth * 2) * (x_to - x_from) + x_from,
                    y = y * (y_to - y_from) + y_from)
}


#long format -----------------------------------------------------------------
dlong <- function(.df, ..., value = NULL) {
  if("..r" %in% names(.df)) stop("The column name '..r' is not allowed")
  .vars <- dplyr::quos(...)
  
  if(!missing(value)) {
    value_var <- dplyr::enquo(value)
    out <- .df %>%
      dplyr::select(!!!.vars, value = !!value_var) %>%
      dplyr::mutate(..r = dplyr::row_number()) %>%
      tidyr::gather(x, node, -..r, -value) %>%
      dplyr::arrange(.data$..r) %>%
      dplyr::group_by(.data$..r) %>%
      dplyr::mutate(next_x = dplyr::lead(.data$x),
                    next_node = dplyr::lead(.data$node)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(-..r) %>%
      dplyr::relocate(value, .after = dplyr::last_col())
  } else {
    out <- .df %>%
      dplyr::select(!!!.vars) %>%
      dplyr::mutate(..r = dplyr::row_number()) %>%
      tidyr::gather(x, node, -..r) %>%
      dplyr::arrange(.data$..r) %>%
      dplyr::group_by(.data$..r) %>%
      dplyr::mutate(next_x = dplyr::lead(.data$x),
                    next_node = dplyr::lead(.data$node)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(-..r)
  }
  
  levels <- unique(out$x)
  
  out %>%
    dplyr::mutate(dplyr::across(c(x, next_x), ~factor(., levels = levels)))
}


#' @title sankey_themes
#' @name theme_sankey
#' @aliases theme_alluvial
#' @aliases theme_sankey_bump
#'
#' @description Minimal themes for sankey, alluvial and sankey bump plots
#'
#' @param base_size base font size, given in pts.
#' @param base_family base font family
#' @param base_line_size base size for line elements
#' @param base_rect_size base size for rect elements
#'
#' @export
theme_sankey <-
  function(base_size = 11,
           base_family = "",
           base_line_size = base_size / 22,
           base_rect_size = base_size / 22) {
    {
      ggplot2::theme_bw(
        base_size = base_size,
        base_family = base_family,
        base_line_size = base_line_size,
        base_rect_size = base_rect_size
      ) %+replace%
        ggplot2::theme(
          panel.border = ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          axis.line = ggplot2::element_line(colour = "black",
                                            size = ggplot2::rel(1)),
          legend.key = ggplot2::element_blank(),
          strip.background = ggplot2::element_rect(
            fill = "white",
            colour = "transparent",
            size = ggplot2::rel(2)
          ),
          complete = TRUE,
          axis.line.y = ggplot2::element_blank(),
          axis.line.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank()
        )
    }
  }

#' @rdname theme_sankey
#' @export
theme_alluvial <-
  function(base_size = 11,
           base_family = "",
           base_line_size = base_size / 22,
           base_rect_size = base_size / 22) {
    {
      ggplot2::theme_bw(
        base_size = base_size,
        base_family = base_family,
        base_line_size = base_line_size,
        base_rect_size = base_rect_size
      ) %+replace%
        ggplot2::theme(
          panel.border = ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          axis.line = ggplot2::element_blank(),
          legend.key = ggplot2::element_blank(),
          strip.background = ggplot2::element_rect(
            fill = "white",
            colour = "transparent",
            size = ggplot2::rel(2)
          ),
          complete = TRUE,
          axis.line.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank()
        )
    }
  }

#' @rdname theme_sankey
#' @export
theme_sankey_bump <-
  function(base_size = 11,
           base_family = "",
           base_line_size = base_size / 22,
           base_rect_size = base_size / 22) {
    {
      ggplot2::theme_bw(
        base_size = base_size,
        base_family = base_family,
        base_line_size = base_line_size,
        base_rect_size = base_rect_size
      ) %+replace%
        ggplot2::theme(
          panel.border = ggplot2::element_blank(),
          panel.grid.major.x = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          axis.line = ggplot2::element_blank(),
          legend.key = ggplot2::element_blank(),
          strip.background = ggplot2::element_rect(
            fill = "white",
            colour = "transparent",
            size = ggplot2::rel(2)
          ),
          complete = TRUE,
          axis.line.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          panel.grid.major.y = ggplot2::element_line("gray90")
        )
    }
  }


# FLOW LAYER ---------
StatSankeyFlow <- ggplot2::ggproto("StatSankeyFlow", ggplot2::Stat,
                                   extra_params = c("n_grid", "na.rm", "type", "width", "space", "smooth"),
                                   
                                   setup_data = function(data, params) {
                                     purrr::map_dfr(unique(data$PANEL),
                                                    ~{
                                                      data <- data %>% dplyr::filter(PANEL == .x)
                                                      
                                                      data <- data %>%
                                                        dplyr::mutate(dplyr::across(c(x, next_x), ~as.numeric(.), .names = ("n_{.col}")))
                                                      
                                                      if(!("value" %in% names(data))) {
                                                        flow_data <- data %>%
                                                          dplyr::mutate(group = 1) %>%
                                                          dplyr::group_by(n_x, node, n_next_x, next_node) %>%
                                                          dplyr::summarise(flow_freq = dplyr::n(), .groups = "keep") %>%
                                                          dplyr::ungroup()
                                                        
                                                        data <- data %>%
                                                          dplyr::mutate(group = 1) %>%
                                                          dplyr::select(-n_next_x, -next_node, -next_x) %>%
                                                          dplyr::group_by_all() %>%
                                                          dplyr::summarise(freq = dplyr::n(), .groups = "keep") %>%
                                                          dplyr::ungroup()
                                                      } else {
                                                        flow_data <- data %>%
                                                          dplyr::mutate(group = 1) %>%
                                                          dplyr::group_by(n_x, node, n_next_x, next_node) %>%
                                                          dplyr::summarise(flow_freq = sum(value, na.rm = TRUE), .groups = "keep") %>%
                                                          dplyr::ungroup()
                                                        
                                                        data <- data %>%
                                                          dplyr::mutate(group = 1) %>%
                                                          dplyr::select(-n_next_x, -next_node, -next_x) %>%
                                                          dplyr::group_by_at(dplyr::vars(dplyr::everything(), -value)) %>%
                                                          dplyr::summarise(freq = sum(value, na.rm = TRUE),, .groups = "keep") %>%
                                                          dplyr::ungroup()
                                                      }
                                                      
                                                      if(is.null(params$space)) {
                                                        params$space <- find_default_space(data)
                                                      }
                                                      
                                                      data <- data %>%
                                                        dplyr::group_by(n_x) %>%
                                                        dplyr::mutate(ymax = cumsum(freq) + (dplyr::row_number() - 1)*params$space,
                                                                      ymin = ymax - freq) %>%
                                                        dplyr::ungroup()
                                                      
                                                      if(params$type == "sankey") {
                                                        data <- data %>%
                                                          dplyr::group_by(n_x) %>%
                                                          dplyr::mutate(ymin = ymin - max(ymax)/2,
                                                                        ymax = ymax - max(ymax)/2) %>%
                                                          dplyr::ungroup()
                                                      } else if (params$type == "alluvial"){
                                                        data <- data
                                                      }
                                                      
                                                      data <- data %>%
                                                        dplyr::mutate(xmin = n_x - params$width/2,
                                                                      xmax = n_x + params$width/2)
                                                      
                                                      if("shift" %in% names(data)) {
                                                        data <- data %>%
                                                          dplyr::mutate(dplyr::across(dplyr::contains("y"), ~ . + shift))
                                                      }
                                                      
                                                      df <- data %>%
                                                        dplyr::left_join(flow_data, by = c("n_x", "node"))
                                                      
                                                      
                                                      
                                                      flows <- df %>%
                                                        dplyr::left_join(df %>%
                                                                           dplyr::select(n_x, node, ymin_end = ymin, ymax_end = ymax, xmin_end = xmin, xmax_end = xmax) %>%
                                                                           dplyr::distinct(),
                                                                         by = c("n_next_x" = "n_x", "next_node" = "node")) %>%
                                                        tidyr::drop_na(n_x, node, next_node, n_next_x, ymax_end, ymin_end, xmax_end, xmin_end) %>%
                                                        dplyr::mutate(r = dplyr::row_number()) %>%
                                                        dplyr::arrange(n_x, -r) %>%
                                                        dplyr::select(-r) %>%
                                                        dplyr::group_by(n_x, node) %>%
                                                        dplyr::mutate(cum_flow_freq = cumsum(flow_freq) - flow_freq) %>%
                                                        dplyr::ungroup() %>%
                                                        dplyr::group_by(n_x, n_next_x, node, next_node) %>%
                                                        dplyr::mutate(flow_start_ymax = ymax - cum_flow_freq,
                                                                      flow_start_ymin = flow_start_ymax - flow_freq)
                                                      
                                                      flows <- flows %>%
                                                        dplyr::arrange(n_x, n_next_x, next_node) %>%
                                                        dplyr::group_by(n_next_x, next_node) %>%
                                                        dplyr::mutate(cum_flow_freq_end = cumsum(flow_freq) - flow_freq) %>%
                                                        dplyr::mutate(flow_end_ymax = ymax_end - cum_flow_freq_end,
                                                                      flow_end_ymin = flow_end_ymax - flow_freq) %>%
                                                        dplyr::ungroup()
                                                      
                                                      flows <- flows %>%
                                                        dplyr::select(-n_x, -node, -freq, -ymax, -ymin, -xmin, -n_next_x, -next_node, -flow_freq, -ymin_end, -ymax_end, -xmax_end, -cum_flow_freq, -cum_flow_freq_end) %>%
                                                        dplyr::mutate(group = dplyr::row_number())
                                                      
                                                      flows %>%
                                                        dplyr::mutate(smooth = params$smooth) %>%
                                                        as.data.frame()
                                                    })
                                     
                                     
                                     
                                   },
                                   
                                   compute_group = function(data, scales) {
                                     
                                     out1 <- sigmoid(data$xmax, data$xmin_end, data$flow_start_ymax, data$flow_end_ymax,
                                                     smooth = data$smooth)
                                     out2 <- sigmoid(data$xmin_end, data$xmax, data$flow_end_ymin, data$flow_start_ymin,
                                                     smooth = data$smooth)
                                     dplyr::bind_rows(out1, out2)
                                   }
)


# FLOW SANKEYBUMP LAYER ---------
StatSankeyBumpFlow <- ggplot2::ggproto("StatSankeyBumpFlow", ggplot2::Stat,
                                       extra_params = c("na.rm", "type", "space", "smooth"),
                                       
                                       setup_data = function(data, params) {
                                         
                                         purrr::map_dfr(unique(data$PANEL),
                                                        ~{
                                                          data <- data %>% dplyr::filter(PANEL == .x)
                                                          
                                                          data <- data %>%
                                                            dplyr::mutate(nodes = paste(node, x)) %>%
                                                            dplyr::arrange(x, -value) %>%
                                                            dplyr::mutate(bbb = dplyr::row_number()) %>%
                                                            dplyr::arrange(bbb) %>%
                                                            dplyr::mutate(nodes = fct_reorder(nodes, value, mean)) %>%
                                                            dplyr::arrange(node, x) %>%
                                                            dplyr::group_by(node) %>%
                                                            dplyr::mutate(next_x = dplyr::lead(x),
                                                                          node = nodes,
                                                                          next_node = dplyr::lead(nodes)) %>%
                                                            dplyr::ungroup() %>%
                                                            dplyr::arrange(x, node)
                                                          
                                                          data <- data %>%
                                                            dplyr::mutate(dplyr::across(c(x, next_x), ~as.numeric(.), .names = ("n_{.col}")))
                                                          
                                                          if(!("value" %in% names(data))) {
                                                            flow_data <- data %>%
                                                              dplyr::mutate(group = 1) %>%
                                                              dplyr::group_by(n_x, node, n_next_x, next_node) %>%
                                                              dplyr::summarise(flow_freq = dplyr::n(), .groups = "keep") %>%
                                                              dplyr::ungroup()
                                                            
                                                            data <- data %>%
                                                              dplyr::mutate(group = 1) %>%
                                                              dplyr::select(-n_next_x, -next_node) %>%
                                                              dplyr::group_by_all() %>%
                                                              dplyr::summarise(freq = dplyr::n(), .groups = "keep") %>%
                                                              dplyr::ungroup()
                                                          } else {
                                                            flow_data <- data %>%
                                                              dplyr::mutate(group = 1) %>%
                                                              dplyr::group_by(n_x, node, n_next_x, next_node) %>%
                                                              dplyr::summarise(flow_freq = sum(value, na.rm = TRUE), .groups = "keep") %>%
                                                              dplyr::ungroup()
                                                            
                                                            data <- data %>%
                                                              dplyr::mutate(group = 1) %>%
                                                              dplyr::select(-n_next_x, -next_node) %>%
                                                              dplyr::group_by_at(dplyr::vars(dplyr::everything(), -value)) %>%
                                                              dplyr::summarise(freq = sum(value, na.rm = TRUE), .groups = "keep") %>%
                                                              dplyr::ungroup()
                                                          }
                                                          
                                                          if(is.null(params$space)) {
                                                            params$space <- find_default_space(data)
                                                          }
                                                          
                                                          data <- data %>%
                                                            dplyr::group_by(n_x) %>%
                                                            dplyr::arrange(node) %>%
                                                            dplyr::mutate(ymax = cumsum(freq) + (dplyr::row_number() - 1)*params$space,
                                                                          ymin = ymax - freq) %>%
                                                            dplyr::ungroup()
                                                          
                                                          if(params$type == "sankey") {
                                                            data <- data %>%
                                                              dplyr::group_by(n_x) %>%
                                                              dplyr::mutate(ymin = ymin - max(ymax)/2,
                                                                            ymax = ymax - max(ymax)/2) %>%
                                                              dplyr::ungroup()
                                                          } else if (params$type == "alluvial"){
                                                            data <- data
                                                          }
                                                          
                                                          data <- data %>%
                                                            dplyr::mutate(xmin = n_x,
                                                                          xmax = n_x)
                                                          
                                                          df <- data %>%
                                                            dplyr::left_join(flow_data, by = c("n_x", "node"))
                                                          
                                                          flows <- df %>%
                                                            dplyr::left_join(df %>%
                                                                               dplyr::select(n_x, node, ymin_end = ymin, ymax_end = ymax, xmin_end = xmin, xmax_end = xmax, flow_freq_end = flow_freq) %>%
                                                                               dplyr::distinct(),
                                                                             by = c("n_next_x" = "n_x", "next_node" = "node")) %>%
                                                            tidyr::drop_na(n_x, node, next_node, n_next_x, ymax_end, ymin_end, xmax_end, xmin_end) %>%
                                                            dplyr::mutate(r = dplyr::row_number()) %>%
                                                            dplyr::arrange(n_x, -r) %>%
                                                            dplyr::select(-r) %>%
                                                            dplyr::group_by(n_x, node) %>%
                                                            dplyr::mutate(cum_flow_freq = cumsum(flow_freq) - flow_freq) %>%
                                                            dplyr::ungroup() %>%
                                                            dplyr::group_by(n_x, n_next_x, node, next_node) %>%
                                                            dplyr::mutate(flow_start_ymax = ymax - cum_flow_freq,
                                                                          flow_start_ymin = flow_start_ymax - flow_freq)
                                                          
                                                          flows <- flows %>%
                                                            dplyr::arrange(n_x, n_next_x, next_node) %>%
                                                            dplyr::group_by(n_next_x, next_node) %>%
                                                            dplyr::mutate(cum_flow_freq_end = cumsum(flow_freq_end) - flow_freq_end) %>%
                                                            dplyr::mutate(flow_end_ymax = ymax_end - cum_flow_freq_end,
                                                                          flow_end_ymin = flow_end_ymax - flow_freq_end) %>%
                                                            dplyr::ungroup()
                                                          
                                                          flows <- flows %>%
                                                            dplyr::select(-n_x, -node, -freq, -ymax, -ymin, -xmin, -n_next_x, -next_node, -flow_freq, -ymin_end, -ymax_end, -xmax_end, -cum_flow_freq, -cum_flow_freq_end) %>%
                                                            dplyr::mutate(group = dplyr::row_number())
                                                          
                                                          flows %>%
                                                            rowwise() %>%
                                                            dplyr::mutate(..groupqq = stringr::str_remove(nodes, as.character(x))) %>%
                                                            dplyr::ungroup() %>%
                                                            dplyr::group_by(..groupqq) %>%
                                                            dplyr::mutate(group = dplyr::cur_group_id()) %>%
                                                            dplyr::ungroup() %>%
                                                            dplyr::select(-..groupqq) %>%
                                                            dplyr::mutate(smooth = params$smooth) %>%
                                                            as.data.frame()
                                                        })
                                       },
                                       
                                       compute_group = function(data, scales) {
                                         
                                         out1 <- purrr::map_dfr(1:nrow(data), ~{
                                           datat <- data %>% dplyr::slice(.x)
                                           sigmoid(datat$xmax, datat$xmin_end, datat$flow_start_ymax, datat$flow_end_ymax,
                                                   smooth = datat$smooth)
                                         }) %>%
                                           dplyr::arrange(x)
                                         out2 <- purrr::map_dfr(1:nrow(data), ~{
                                           datat <- data %>% dplyr::slice(.x)
                                           sigmoid(datat$xmin_end, datat$xmax, datat$flow_end_ymin, datat$flow_start_ymin,
                                                   smooth = datat$smooth)
                                         }) %>%
                                           dplyr::arrange(-x)
                                         
                                         dplyr::bind_rows(out1, out2)
                                       }
)

# TEXT LAYER -------
StatSankeyText <- ggplot2::ggproto("StatSankeyText", ggplot2::Stat,
                                   extra_params = c("n_grid", "na.rm", "type", "width", "space"),
                                   
                                   setup_data = function(data, params) {
                                     
                                     purrr::map_dfr(unique(data$PANEL),
                                                    ~{
                                                      data <- data %>% dplyr::filter(PANEL == .x)
                                                      
                                                      data <- data %>%
                                                        dplyr::mutate(dplyr::across(c(x, next_x), ~as.numeric(.), .names = ("n_{.col}")))
                                                      
                                                      if(!("value" %in% names(data))) {
                                                        data <- data %>%
                                                          dplyr::mutate(group = 1) %>%
                                                          dplyr::select(-n_next_x, -next_node) %>%
                                                          dplyr::group_by_all() %>%
                                                          dplyr::summarise(freq = dplyr::n(), .groups = "keep") %>%
                                                          dplyr::ungroup()
                                                      } else {
                                                        data <- data %>%
                                                          dplyr::mutate(group = 1) %>%
                                                          dplyr::select(-n_next_x, -next_node) %>%
                                                          dplyr::group_by_at(dplyr::vars(dplyr::everything(), -value)) %>%
                                                          dplyr::summarise(freq = sum(value, na.rm = TRUE), .groups = "keep") %>%
                                                          dplyr::ungroup()
                                                      }
                                                      
                                                      if(is.null(params$space)) {
                                                        params$space <- find_default_space(data)
                                                      }
                                                      
                                                      data <- data %>%
                                                        dplyr::group_by(n_x) %>%
                                                        dplyr::mutate(ymax = cumsum(freq) + (dplyr::row_number() - 1)*params$space,
                                                                      ymin = ymax - freq) %>%
                                                        dplyr::ungroup()
                                                      
                                                      if(params$type == "sankey") {
                                                        data <- data %>%
                                                          dplyr::group_by(n_x) %>%
                                                          dplyr::mutate(ymin = ymin - max(ymax)/2,
                                                                        ymax = ymax - max(ymax)/2) %>%
                                                          dplyr::ungroup()
                                                      } else if (params$type == "alluvial"){
                                                        data <- data
                                                      }
                                                      
                                                      data <- data %>%
                                                        dplyr::mutate(xmin = n_x - params$width/2,
                                                                      xmax = n_x + params$width/2)
                                                      
                                                      data <- data %>%
                                                        dplyr::mutate(x = n_x,
                                                                      y = ymin + (ymax - ymin)/2)
                                                      
                                                      if("shift" %in% names(data)) {
                                                        data <- data %>%
                                                          dplyr::mutate(dplyr::across(dplyr::contains("y"), ~ . + shift))
                                                      }
                                                      
                                                      
                                                      return(as.data.frame(data))
                                                    })
                                   },
                                   
                                   compute_group = function(data, scales) {
                                     data
                                   }
)


# NODE LAYER -------
StatSankeyNode <- ggplot2::ggproto("StatSankeyNode", ggplot2::Stat,
                                   extra_params = c("n_grid", "na.rm", "type", "width", "space", "smooth"),
                                   
                                   setup_data = function(data, params) {
                                     
                                     purrr::map_dfr(unique(data$PANEL),
                                                    ~{
                                                      
                                                      data <- data %>% dplyr::filter(PANEL == .x)
                                                      data <- data %>%
                                                        dplyr::mutate(dplyr::across(c(x, next_x), ~as.numeric(.), .names = ("n_{.col}")))
                                                      
                                                      if(!("value" %in% names(data))) {
                                                        data <- data %>%
                                                          dplyr::mutate(group = 1) %>%
                                                          dplyr::select(-n_next_x, -next_node, -next_x) %>%
                                                          dplyr::group_by_all() %>%
                                                          dplyr::summarise(freq = dplyr::n(), .groups = "keep") %>%
                                                          dplyr::ungroup()
                                                      } else {
                                                        data <- data %>%
                                                          dplyr::mutate(group = 1) %>%
                                                          dplyr::select(-n_next_x, -next_node, -next_x) %>%
                                                          dplyr::group_by_at(dplyr::vars(dplyr::everything(), -value)) %>%
                                                          dplyr::summarise(freq = sum(value, na.rm = TRUE), .groups = "keep") %>%
                                                          dplyr::ungroup()
                                                      }
                                                      
                                                      if(is.null(params$space)) {
                                                        params$space <- find_default_space(data)
                                                      }
                                                      
                                                      data <- data %>%
                                                        dplyr::group_by(n_x) %>%
                                                        dplyr::mutate(ymax = cumsum(freq) + (dplyr::row_number() - 1)*params$space,
                                                                      ymin = ymax - freq) %>%
                                                        dplyr::ungroup()
                                                      
                                                      if(params$type == "sankey") {
                                                        data <- data %>%
                                                          dplyr::group_by(n_x) %>%
                                                          dplyr::mutate(ymin = ymin - max(ymax)/2,
                                                                        ymax = ymax - max(ymax)/2) %>%
                                                          dplyr::ungroup()
                                                      } else if (params$type == "alluvial"){
                                                        data <- data
                                                      }
                                                      
                                                      data <- data %>%
                                                        dplyr::mutate(xmin = n_x - params$width/2,
                                                                      xmax = n_x + params$width/2)
                                                      
                                                      if("shift" %in% names(data)) {
                                                        data <- data %>%
                                                          dplyr::mutate(dplyr::across(dplyr::contains("y"), ~ . + shift))
                                                      }
                                                      
                                                      return(as.data.frame(data))
                                                    })
                                     
                                   },
                                   
                                   compute_group = function(data, scales) {
                                     data
                                   }
)


sankey_p <- function(mapping = NULL,
                        data = NULL,
                        position = "identity",
                        na.rm = FALSE,
                        show.legend = NA,
                        space = NULL,
                        type = "sankey",
                        width = .1,
                        smooth = 8,
                        inherit.aes = TRUE,
                        ...
) {
  params_list <- prepare_params(...)
  
  list(
    flow = ggplot2::layer(
      stat = StatSankeyFlow,
      data = data,
      mapping = mapping,
      geom = "polygon",
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = purrr::flatten(
        list(
          na.rm = na.rm,
          width = width,
          space = space,
          smooth = smooth,
          type = type,
          params_list[[1]]
        )
      )
    ),
    
    node = ggplot2::layer(
      stat = StatSankeyNode,
      data = data,
      mapping = mapping,
      geom = ggplot2::GeomRect,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = purrr::flatten(
        list(
          na.rm = na.rm,
          width = width,
          space = space,
          smooth = smooth,
          type = type,
          params_list[[2]]
        )
      )
    )
  )
  
  
}


sankey_p_label <- function(mapping = NULL,
                              data = NULL,
                              position = "identity",
                              na.rm = FALSE,
                              show.legend = NA,
                              space = NULL,
                              type = "sankey",
                              width = .1,
                              inherit.aes = TRUE,
                              ...) {
  # Prepare aesthics for label
  label.aes <- list(...)
  
  list(
    label = ggplot2::layer(
      stat = StatSankeyText,
      data = data,
      mapping = mapping,
      geom = "label",
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = purrr::flatten(
        list(
          na.rm = na.rm,
          width = width,
          space = space,
          type = type,
          label.aes
        )
      )
    )
  )
}
