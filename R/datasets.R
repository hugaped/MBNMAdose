#' Studies of triptans for headache pain relief
#'
#' A dataset from a systematic review of interventions for pain relief in migraine \insertCite{thorlund2014}{MBNMAdose}.
#' The outcome is binary, and represents (as aggregate data) the proportion of participants who were
#' headache-free at 2 hours. Data are from patients who had had at least one migraine attack, who were
#' not lost to follow-up, and who did not violate the trial protocol. The dataset includes 70 Randomised-Controlled
#' Trials (RCTs), comparing 7 triptans with placebo. Doses are standardised as relative to a "common" dose,
#' and in total there are 23 different treatments (combination of dose and agent).
#'
#' @format A data frame in long format (one row per arm and study), with with 181 rows and 6 variables:
#' * `studyID` Study identifiers
#' * `AuthorYear` The author and year published of the study
#' * `N` Numeric data indicating the number of participants in a study arm
#' * `r` Numeric data indicating the number of responders (headache free at 2 hours) in a study arm
#' * `dose` Numeric data indicating the standardised dose received
#' * `agent` Factor data indicating the agent to which participants were randomised
#'
#' @source
#' \insertAllCited{}
"HF2PPITT"




#' Studies of treatments for Serum Uric Acid reduction in patients with gout
#'
#' A dataset from a systematic review of interventions for lowering Serum Uric Acid (SUA) concentration in
#' patients with gout **(not published previously)**. The outcome is continuous, and aggregate data responses
#' correspond to the mean change from baseline in SUA in mg/dL at 2 weeks follow-up. The dataset includes 10
#' Randomised-Controlled Trials (RCTs), comparing 5 different agents, and placebo. Data for one agent (RDEA)
#' arises from an RCT that is not placebo-controlled, and so is not connected to the network directly. In total
#' there were 19 different treatments (combination of dose and agent).
#'
#' @format A data frame in long format (one row per arm and study), with 27 rows and 5 variables:
#' * `studyID` Study identifiers
#' * `y` Numeric data indicating the mean change from baseline in SUA in a study arm
#' * `se` Numeric data indicating the standard error for the mean change from baseline in SUA in a study arm
#' * `agent` Character data indicating the agent to which participants were randomised
#' * `dose` Numeric data indicating the standardised dose received
#'
#' @source Pfizer Ltd.
"GoutSUA_2wkCFB"




#' Studies of treatments for pain relief in patients with osteoarthritis
#'
#' A dataset from a systematic review of interventions for pain relief in osteoarthritis, used previously
#' in \insertCite{pedder2019;textual}{MBNMAdose}. The outcome is continuous, and aggregate data responses correspond to
#' the mean WOMAC pain score at 2 weeks follow-up. The dataset includes 18 Randomised-Controlled Trials (RCTs),
#' comparing 8 different agents with placebo. In total there were 26 different treatments (combination of dose and
#' agent). The active treatments can also be grouped into 3 different classes, within which they have similar
#' mechanisms of action.
#'
#'
#' @format A data frame in long format (one row per arm and study), with 74 rows and 7 variables:
#' * `studyID` Study identifiers
#' * `agent` Character data indicating the agent to which participants were randomised
#' * `dose` Numeric data indicating the standardised dose received
#' * `class` Character data indicating the drug class to which the agent belongs to
#' * `y` Numeric data indicating the mean pain score on the WOMAC scale in a study arm
#' * `se` Numeric data indicating the standard error for the mean pain score on the WOMAC scale in a study arm
#'
#' @references
#' \insertAllCited{}
#'
#' @source Pfizer Ltd.
"osteopain_2wkabs"






#' Studies of alogliptin for lowering blood glucose concentration in patients with type II diabetes
#'
#' A dataset from a systematic review of Randomised-Controlled Trials (RCTs) comparing different doses of
#' alogliptin with placebo \insertCite{langford2016}{MBNMAdose}. The systematic review was simply performed and was intended to
#' provide data to illustrate a statistical methodology rather than for clinical inference. Alogliptin is
#' a treatment aimed at reducing blood glucose concentration in type II diabetes. The outcome is continuous,
#' and aggregate data responses correspond to the mean change in HbA1c from baseline to follow-up in studies
#' of at least 12 weeks follow-up. The dataset includes 14 Randomised-Controlled Trials (RCTs), comparing 5
#' different doses of alogliptin with placebo, leading to 6 different treatments (combination of dose and agent)
#' within the network.
#'
#' `alog_pcfb` is a data frame in long format (one row per arm and study), with the variables `studyID`, `agent`, `dose`, `y`, `se`, and `N`.
#'
#' @format A data frame in long format (one row per arm and study), with 46 rows and 6 variables:
#' * `studyID` Study identifiers
#' * `agent` Character data indicating the agent to which participants were randomised
#' * `dose` Numeric data indicating the standardised dose received
#' * `y` Numeric data indicating the mean change from baseline in blood glucose concentration (mg/dL) in a study arm
#' * `se` Numeric data indicating the standard error for the mean change from baseline in blood glucose concentration (mg/dL) in a study arm
#'
#' @references
#' \insertAllCited{}
#'
"alog_pcfb"

