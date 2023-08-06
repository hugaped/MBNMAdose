# Global variables for MBNMAdose
# Author: Hugo Pedder
# Date created: 2022-02-22

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c(".", "studyID", "agent", "dose", "Var1", "value",
                                                        "Parameter", "fupdose", "groupvar", "y", "sd", "stansd",
                                                        "network", "a", "param", "med", "l95", "u95", "value",
                                                        "Estimate", "2.5%", "50%", "97.5%", "treatment",
                                                        "study", "arm", "mod1.mean", "mod2.mean", "argcheck",
                                                        "Comparison", "Evidence", "splinefun", "level",
                                                        "user.str"))
