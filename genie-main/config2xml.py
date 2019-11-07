#!/usr/bin/python

import sys
import os.path
import re

# initialisation...
vars = {}
varlist = ["EXPID","RESTARTREAD","CHECKFLUXES"]

# NB "ma_" -> <control> is special
# GEM common is included in control section
# initialise with default science module choices
configflags = {
    # flag name in control,shortname,defual config, prefix 
    # atmos
    "ma_flag_ebatmos" : ["embm",".FALSE.","ea"],
    "ma_flag_goldsteinocean" : ["goldstein",".FALSE.","go"],
    "ma_flag_goldsteinseaice" : ["goldsteinseaice",".FALSE.","gs"],
    "ma_flag_ents" : ["ents",".FALSE.","el"],
    "ma_flag_biogem" : ["biogem",".FALSE.","bg"],
    "ma_flag_atchem" : ["atchem",".FALSE.","ac"],
    "ma_flag_sedgem" : ["sedgem",".FALSE.","sg"],
    "ma_flag_rokgem" : ["rokgem",".FALSE.","rg"],
    "ma_flag_ichem" : ["ichem",".FALSE.","ci"]
    "gl_flag_gemlite" : ["gemlite",".TRUE.","gl"],
    }

control = {} 

# paramarrays!
# gem
atm_select = {
    "gm_atm_select_1" : ["1",""],
    "gm_atm_select_2" : ["2",""],
    "gm_atm_select_3" : ["3",""],
    "gm_atm_select_4" : ["4",""],
    "gm_atm_select_5" : ["5",""],
    "gm_atm_select_6" : ["6",""],
    "gm_atm_select_7" : ["7",""],
    "gm_atm_select_8" : ["8",""],
    "gm_atm_select_9" : ["9",""],
    "gm_atm_select_10" : ["10",""],
    "gm_atm_select_11" : ["11",""],
    "gm_atm_select_12" : ["12",""],
    "gm_atm_select_13" : ["13",""],
    "gm_atm_select_14" : ["14",""],
    "gm_atm_select_15" : ["15",""],
    "gm_atm_select_16" : ["16",""],
    "gm_atm_select_17" : ["17",""],
    "gm_atm_select_18" : ["18",""],
    "gm_atm_select_19" : ["19",""]
    }
ocn_select = {
    "gm_ocn_select_1" : ["1",""],
    "gm_ocn_select_2" : ["2",""],
    "gm_ocn_select_3" : ["3",""],
    "gm_ocn_select_4" : ["4",""],
    "gm_ocn_select_5" : ["5",""],
    "gm_ocn_select_6" : ["6",""],
    "gm_ocn_select_7" : ["7",""],
    "gm_ocn_select_8" : ["8",""],
    "gm_ocn_select_9" : ["9",""],
    "gm_ocn_select_10" : ["10",""],
    "gm_ocn_select_11" : ["11",""],
    "gm_ocn_select_12" : ["12",""],
    "gm_ocn_select_13" : ["13",""],
    "gm_ocn_select_14" : ["14",""],
    "gm_ocn_select_15" : ["15",""],
    "gm_ocn_select_16" : ["16",""],
    "gm_ocn_select_17" : ["17",""],
    "gm_ocn_select_18" : ["18",""],
    "gm_ocn_select_19" : ["19",""],
    "gm_ocn_select_20" : ["20",""],
    "gm_ocn_select_21" : ["21",""],
    "gm_ocn_select_22" : ["22",""],
    "gm_ocn_select_23" : ["23",""],
    "gm_ocn_select_24" : ["24",""],
    "gm_ocn_select_25" : ["25",""],
    "gm_ocn_select_26" : ["26",""],
    "gm_ocn_select_27" : ["27",""],
    "gm_ocn_select_28" : ["28",""],
    "gm_ocn_select_29" : ["29",""],
    "gm_ocn_select_30" : ["30",""],
    "gm_ocn_select_31" : ["31",""],
    "gm_ocn_select_32" : ["32",""],
    "gm_ocn_select_33" : ["33",""],
    "gm_ocn_select_34" : ["34",""],
    "gm_ocn_select_35" : ["35",""],
    "gm_ocn_select_36" : ["36",""],
    "gm_ocn_select_37" : ["37",""],
    "gm_ocn_select_38" : ["38",""],
    "gm_ocn_select_39" : ["39",""],
    "gm_ocn_select_40" : ["40",""],
    "gm_ocn_select_41" : ["41",""],
    "gm_ocn_select_42" : ["42",""],
    "gm_ocn_select_43" : ["43",""],
    "gm_ocn_select_44" : ["44",""],
    "gm_ocn_select_45" : ["45",""],
    "gm_ocn_select_46" : ["46",""],
    "gm_ocn_select_47" : ["47",""],
    "gm_ocn_select_48" : ["48",""],
    "gm_ocn_select_49" : ["49",""],
    "gm_ocn_select_50" : ["50",""],
    "gm_ocn_select_51" : ["51",""],
    "gm_ocn_select_52" : ["52",""],
    "gm_ocn_select_53" : ["53",""],
    "gm_ocn_select_54" : ["54",""],
    "gm_ocn_select_55" : ["55",""],
    "gm_ocn_select_56" : ["56",""],
# NH4_new (Fanny, June 2015)
    "gm_ocn_select_76" : ["76",""],
# NO3_new (Fanny, Sep 2015)
    "gm_ocn_select_77" : ["77",""]
    }
sed_select = {
    "gm_sed_select_1" : ["1",""],
    "gm_sed_select_2" : ["2",""],
    "gm_sed_select_3" : ["3",""],
    "gm_sed_select_4" : ["4",""],
    "gm_sed_select_5" : ["5",""],
    "gm_sed_select_6" : ["6",""],
    "gm_sed_select_7" : ["7",""],
    "gm_sed_select_8" : ["8",""],
    "gm_sed_select_9" : ["9",""],
    "gm_sed_select_10" : ["10",""],
    "gm_sed_select_11" : ["11",""],
    "gm_sed_select_12" : ["12",""],
    "gm_sed_select_13" : ["13",""],
    "gm_sed_select_14" : ["14",""],
    "gm_sed_select_15" : ["15",""],
    "gm_sed_select_16" : ["16",""],
    "gm_sed_select_17" : ["17",""],
    "gm_sed_select_18" : ["18",""],
    "gm_sed_select_19" : ["19",""],
    "gm_sed_select_20" : ["20",""],
    "gm_sed_select_21" : ["21",""],
    "gm_sed_select_22" : ["22",""],
    "gm_sed_select_23" : ["23",""],
    "gm_sed_select_24" : ["24",""],
    "gm_sed_select_25" : ["25",""],
    "gm_sed_select_26" : ["26",""],
    "gm_sed_select_27" : ["27",""],
    "gm_sed_select_28" : ["28",""],
    "gm_sed_select_29" : ["29",""],
    "gm_sed_select_30" : ["30",""],
    "gm_sed_select_31" : ["31",""],
    "gm_sed_select_32" : ["32",""],
    "gm_sed_select_33" : ["33",""],
    "gm_sed_select_34" : ["34",""],
    "gm_sed_select_35" : ["35",""],
    "gm_sed_select_36" : ["36",""],
    "gm_sed_select_37" : ["37",""],
    "gm_sed_select_38" : ["38",""],
    "gm_sed_select_39" : ["39",""],
    "gm_sed_select_40" : ["40",""],
    "gm_sed_select_41" : ["41",""],
    "gm_sed_select_42" : ["42",""],
    "gm_sed_select_43" : ["43",""],
    "gm_sed_select_44" : ["44",""],
    "gm_sed_select_45" : ["45",""],
    "gm_sed_select_46" : ["46",""],
    "gm_sed_select_47" : ["47",""],
    "gm_sed_select_48" : ["48",""],
    "gm_sed_select_49" : ["49",""],
    "gm_sed_select_50" : ["50",""],
    "gm_sed_select_51" : ["51",""],
    "gm_sed_select_52" : ["52",""],
    "gm_sed_select_53" : ["53",""],
    "gm_sed_select_54" : ["54",""],
    "gm_sed_select_55" : ["55",""],
    "gm_sed_select_56" : ["56",""]
    }
# atchem
atm_init = {
    "ac_atm_init_1" : ["1",""],
    "ac_atm_init_2" : ["2",""],
    "ac_atm_init_3" : ["3",""],
    "ac_atm_init_4" : ["4",""],
    "ac_atm_init_5" : ["5",""],
    "ac_atm_init_6" : ["6",""],
    "ac_atm_init_7" : ["7",""],
    "ac_atm_init_8" : ["8",""],
    "ac_atm_init_9" : ["9",""],
    "ac_atm_init_10" : ["10",""],
    "ac_atm_init_11" : ["11",""],
    "ac_atm_init_12" : ["12",""],
    "ac_atm_init_13" : ["13",""],
    "ac_atm_init_14" : ["14",""],
    "ac_atm_init_15" : ["15",""],
    "ac_atm_init_16" : ["16",""],
    "ac_atm_init_17" : ["17",""],
    "ac_atm_init_18" : ["18",""],
    "ac_atm_init_19" : ["19",""]
    }
# biogem
ocn_init = {
    "bg_ocn_init_1" : ["1",""],
    "bg_ocn_init_2" : ["2",""],
    "bg_ocn_init_3" : ["3",""],
    "bg_ocn_init_4" : ["4",""],
    "bg_ocn_init_5" : ["5",""],
    "bg_ocn_init_6" : ["6",""],
    "bg_ocn_init_7" : ["7",""],
    "bg_ocn_init_8" : ["8",""],
    "bg_ocn_init_9" : ["9",""],
    "bg_ocn_init_10" : ["10",""],
    "bg_ocn_init_11" : ["11",""],
    "bg_ocn_init_12" : ["12",""],
    "bg_ocn_init_13" : ["13",""],
    "bg_ocn_init_14" : ["14",""],
    "bg_ocn_init_15" : ["15",""],
    "bg_ocn_init_16" : ["16",""],
    "bg_ocn_init_17" : ["17",""],
    "bg_ocn_init_18" : ["18",""],
    "bg_ocn_init_19" : ["19",""],
    "bg_ocn_init_20" : ["20",""],
    "bg_ocn_init_21" : ["21",""],
    "bg_ocn_init_22" : ["22",""],
    "bg_ocn_init_23" : ["23",""],
    "bg_ocn_init_24" : ["24",""],
    "bg_ocn_init_25" : ["25",""],
    "bg_ocn_init_26" : ["26",""],
    "bg_ocn_init_27" : ["27",""],
    "bg_ocn_init_28" : ["28",""],
    "bg_ocn_init_29" : ["29",""],
    "bg_ocn_init_30" : ["30",""],
    "bg_ocn_init_31" : ["31",""],
    "bg_ocn_init_32" : ["32",""],
    "bg_ocn_init_33" : ["33",""],
    "bg_ocn_init_34" : ["34",""],
    "bg_ocn_init_35" : ["35",""],
    "bg_ocn_init_36" : ["36",""],
    "bg_ocn_init_37" : ["37",""],
    "bg_ocn_init_38" : ["38",""],
    "bg_ocn_init_39" : ["39",""],
    "bg_ocn_init_40" : ["40",""],
    "bg_ocn_init_41" : ["41",""],
    "bg_ocn_init_42" : ["42",""],
    "bg_ocn_init_43" : ["43",""],
    "bg_ocn_init_44" : ["44",""],
    "bg_ocn_init_45" : ["45",""],
    "bg_ocn_init_46" : ["46",""],
    "bg_ocn_init_47" : ["47",""],
    "bg_ocn_init_48" : ["48",""],
    "bg_ocn_init_49" : ["49",""],
    "bg_ocn_init_50" : ["50",""],
    "bg_ocn_init_51" : ["51",""],
    "bg_ocn_init_52" : ["52",""],
    "bg_ocn_init_53" : ["53",""],
    "bg_ocn_init_54" : ["54",""],
    "bg_ocn_init_55" : ["55",""],
    "bg_ocn_init_56" : ["56",""],
# NH4_new (Fanny, June 2015)
    "bg_ocn_init_76" : ["76",""],
# NO3_new (Fanny, Sep 2015)
    "bg_ocn_init_77" : ["77",""]
    }
par_atm_force_scale_val = {
    "bg_par_atm_force_scale_val_1" : ["1",""],
    "bg_par_atm_force_scale_val_2" : ["2",""],
    "bg_par_atm_force_scale_val_3" : ["3",""],
    "bg_par_atm_force_scale_val_4" : ["4",""],
    "bg_par_atm_force_scale_val_5" : ["5",""],
    "bg_par_atm_force_scale_val_6" : ["6",""],
    "bg_par_atm_force_scale_val_7" : ["7",""],
    "bg_par_atm_force_scale_val_8" : ["8",""],
    "bg_par_atm_force_scale_val_9" : ["9",""],
    "bg_par_atm_force_scale_val_10" : ["10",""],
    "bg_par_atm_force_scale_val_11" : ["11",""],
    "bg_par_atm_force_scale_val_12" : ["12",""],
    "bg_par_atm_force_scale_val_13" : ["13",""],
    "bg_par_atm_force_scale_val_14" : ["14",""],
    "bg_par_atm_force_scale_val_15" : ["15",""],
    "bg_par_atm_force_scale_val_16" : ["16",""],
    "bg_par_atm_force_scale_val_17" : ["17",""],
    "bg_par_atm_force_scale_val_18" : ["18",""],
    "bg_par_atm_force_scale_val_19" : ["19",""]
    }
par_atm_force_scale_time = {
    "bg_par_atm_force_scale_time_1" : ["1",""],
    "bg_par_atm_force_scale_time_2" : ["2",""],
    "bg_par_atm_force_scale_time_3" : ["3",""],
    "bg_par_atm_force_scale_time_4" : ["4",""],
    "bg_par_atm_force_scale_time_5" : ["5",""],
    "bg_par_atm_force_scale_time_6" : ["6",""],
    "bg_par_atm_force_scale_time_7" : ["7",""],
    "bg_par_atm_force_scale_time_8" : ["8",""],
    "bg_par_atm_force_scale_time_9" : ["9",""],
    "bg_par_atm_force_scale_time_10" : ["10",""],
    "bg_par_atm_force_scale_time_11" : ["11",""],
    "bg_par_atm_force_scale_time_12" : ["12",""],
    "bg_par_atm_force_scale_time_13" : ["13",""],
    "bg_par_atm_force_scale_time_14" : ["14",""],
    "bg_par_atm_force_scale_time_15" : ["15",""],
    "bg_par_atm_force_scale_time_16" : ["16",""],
    "bg_par_atm_force_scale_time_17" : ["17",""],
    "bg_par_atm_force_scale_time_18" : ["18",""],
    "bg_par_atm_force_scale_time_19" : ["19",""]
    }
par_ocn_force_scale_val = {
    "bg_par_ocn_force_scale_val_1" : ["1",""],
    "bg_par_ocn_force_scale_val_2" : ["2",""],
    "bg_par_ocn_force_scale_val_3" : ["3",""],
    "bg_par_ocn_force_scale_val_4" : ["4",""],
    "bg_par_ocn_force_scale_val_5" : ["5",""],
    "bg_par_ocn_force_scale_val_6" : ["6",""],
    "bg_par_ocn_force_scale_val_7" : ["7",""],
    "bg_par_ocn_force_scale_val_8" : ["8",""],
    "bg_par_ocn_force_scale_val_9" : ["9",""],
    "bg_par_ocn_force_scale_val_10" : ["10",""],
    "bg_par_ocn_force_scale_val_11" : ["11",""],
    "bg_par_ocn_force_scale_val_12" : ["12",""],
    "bg_par_ocn_force_scale_val_13" : ["13",""],
    "bg_par_ocn_force_scale_val_14" : ["14",""],
    "bg_par_ocn_force_scale_val_15" : ["15",""],
    "bg_par_ocn_force_scale_val_16" : ["16",""],
    "bg_par_ocn_force_scale_val_17" : ["17",""],
    "bg_par_ocn_force_scale_val_18" : ["18",""],
    "bg_par_ocn_force_scale_val_19" : ["19",""],
    "bg_par_ocn_force_scale_val_20" : ["20",""],
    "bg_par_ocn_force_scale_val_21" : ["21",""],
    "bg_par_ocn_force_scale_val_22" : ["22",""],
    "bg_par_ocn_force_scale_val_23" : ["23",""],
    "bg_par_ocn_force_scale_val_24" : ["24",""],
    "bg_par_ocn_force_scale_val_25" : ["25",""],
    "bg_par_ocn_force_scale_val_26" : ["26",""],
    "bg_par_ocn_force_scale_val_27" : ["27",""],
    "bg_par_ocn_force_scale_val_28" : ["28",""],
    "bg_par_ocn_force_scale_val_29" : ["29",""],
    "bg_par_ocn_force_scale_val_30" : ["30",""],
    "bg_par_ocn_force_scale_val_31" : ["31",""],
    "bg_par_ocn_force_scale_val_32" : ["32",""],
    "bg_par_ocn_force_scale_val_33" : ["33",""],
    "bg_par_ocn_force_scale_val_34" : ["34",""],
    "bg_par_ocn_force_scale_val_35" : ["35",""],
    "bg_par_ocn_force_scale_val_36" : ["36",""],
    "bg_par_ocn_force_scale_val_37" : ["37",""],
    "bg_par_ocn_force_scale_val_38" : ["38",""],
    "bg_par_ocn_force_scale_val_39" : ["39",""],
    "bg_par_ocn_force_scale_val_40" : ["40",""],
    "bg_par_ocn_force_scale_val_41" : ["41",""],
    "bg_par_ocn_force_scale_val_42" : ["42",""],
    "bg_par_ocn_force_scale_val_43" : ["43",""],
    "bg_par_ocn_force_scale_val_44" : ["44",""],
    "bg_par_ocn_force_scale_val_45" : ["45",""],
    "bg_par_ocn_force_scale_val_46" : ["46",""],
    "bg_par_ocn_force_scale_val_47" : ["47",""],
    "bg_par_ocn_force_scale_val_48" : ["48",""],
    "bg_par_ocn_force_scale_val_49" : ["49",""],
    "bg_par_ocn_force_scale_val_50" : ["50",""],
    "bg_par_ocn_force_scale_val_51" : ["51",""],
    "bg_par_ocn_force_scale_val_52" : ["52",""],
    "bg_par_ocn_force_scale_val_53" : ["53",""],
    "bg_par_ocn_force_scale_val_54" : ["54",""],
    "bg_par_ocn_force_scale_val_55" : ["55",""],
    "bg_par_ocn_force_scale_val_56" : ["56",""]
    }
par_ocn_force_scale_time = {
    "bg_par_ocn_force_scale_time_1" : ["1",""],
    "bg_par_ocn_force_scale_time_2" : ["2",""],
    "bg_par_ocn_force_scale_time_3" : ["3",""],
    "bg_par_ocn_force_scale_time_4" : ["4",""],
    "bg_par_ocn_force_scale_time_5" : ["5",""],
    "bg_par_ocn_force_scale_time_6" : ["6",""],
    "bg_par_ocn_force_scale_time_7" : ["7",""],
    "bg_par_ocn_force_scale_time_8" : ["8",""],
    "bg_par_ocn_force_scale_time_9" : ["9",""],
    "bg_par_ocn_force_scale_time_10" : ["10",""],
    "bg_par_ocn_force_scale_time_11" : ["11",""],
    "bg_par_ocn_force_scale_time_12" : ["12",""],
    "bg_par_ocn_force_scale_time_13" : ["13",""],
    "bg_par_ocn_force_scale_time_14" : ["14",""],
    "bg_par_ocn_force_scale_time_15" : ["15",""],
    "bg_par_ocn_force_scale_time_16" : ["16",""],
    "bg_par_ocn_force_scale_time_17" : ["17",""],
    "bg_par_ocn_force_scale_time_18" : ["18",""],
    "bg_par_ocn_force_scale_time_19" : ["19",""],
    "bg_par_ocn_force_scale_time_20" : ["20",""],
    "bg_par_ocn_force_scale_time_21" : ["21",""],
    "bg_par_ocn_force_scale_time_22" : ["22",""],
    "bg_par_ocn_force_scale_time_23" : ["23",""],
    "bg_par_ocn_force_scale_time_24" : ["24",""],
    "bg_par_ocn_force_scale_time_25" : ["25",""],
    "bg_par_ocn_force_scale_time_26" : ["26",""],
    "bg_par_ocn_force_scale_time_27" : ["27",""],
    "bg_par_ocn_force_scale_time_28" : ["28",""],
    "bg_par_ocn_force_scale_time_29" : ["29",""],
    "bg_par_ocn_force_scale_time_30" : ["30",""],
    "bg_par_ocn_force_scale_time_31" : ["31",""],
    "bg_par_ocn_force_scale_time_32" : ["32",""],
    "bg_par_ocn_force_scale_time_33" : ["33",""],
    "bg_par_ocn_force_scale_time_34" : ["34",""],
    "bg_par_ocn_force_scale_time_35" : ["35",""],
    "bg_par_ocn_force_scale_time_36" : ["36",""],
    "bg_par_ocn_force_scale_time_37" : ["37",""],
    "bg_par_ocn_force_scale_time_38" : ["38",""],
    "bg_par_ocn_force_scale_time_39" : ["39",""],
    "bg_par_ocn_force_scale_time_40" : ["40",""],
    "bg_par_ocn_force_scale_time_41" : ["41",""],
    "bg_par_ocn_force_scale_time_42" : ["42",""],
    "bg_par_ocn_force_scale_time_43" : ["43",""],
    "bg_par_ocn_force_scale_time_44" : ["44",""],
    "bg_par_ocn_force_scale_time_45" : ["45",""],
    "bg_par_ocn_force_scale_time_46" : ["46",""],
    "bg_par_ocn_force_scale_time_47" : ["47",""],
    "bg_par_ocn_force_scale_time_48" : ["48",""],
    "bg_par_ocn_force_scale_time_49" : ["49",""],
    "bg_par_ocn_force_scale_time_50" : ["50",""],
    "bg_par_ocn_force_scale_time_51" : ["51",""],
    "bg_par_ocn_force_scale_time_52" : ["52",""],
    "bg_par_ocn_force_scale_time_53" : ["53",""],
    "bg_par_ocn_force_scale_time_54" : ["54",""],
    "bg_par_ocn_force_scale_time_55" : ["55",""],
    "bg_par_ocn_force_scale_time_56" : ["56",""]
    }

# embm
diffamp = {
    "ea_12" : ["1",""],
    "ea_13" : ["2",""]
    }
betaz = {
    "ea_16" : ["1",""],
    "ea_18" : ["2",""]
    }
betam = {
    "ea_17" : ["1",""],
    "ea_19" : ["2",""]
    }
# goldstein
diff = {
    "go_14" : ["1",""],
    "go_15" : ["2",""]
    }

deletelist = []

makeargs = {}
macros = {}
makearglist = ["GENIEDP","IGCMATMOSDP","FLAG_GLIMMER","FLAG_MOSESTRIFFID"]
macrolist = ["GENIENXOPTS","GENIENYOPTS","GENIENLOPTS"]
macrolist = macrolist + ["GOLDSTEINNLONSOPTS","GOLDSTEINNLATSOPTS","GOLDSTEINNLEVSOPTS","GOLDSTEINMAXISLESOPTS"]
macrolist = macrolist + ["GOLDSTEINNTRACSOPTS"]
macrolist = macrolist + ["GENIELEVRFOPTS","GENIEMXLEVOPTS"]
macrolist = macrolist + ["SEAICEOPTS"]
macrolist = macrolist + ["IGCMNWJ2OPTS","IGCMNNOPTS","IGCMMMOPTS","IGCMPQSATTESTOPTS"]
macrolist = macrolist + ["SEDGEMNLONSOPTS","SEDGEMNLATSOPTS","SEDGEMNLONSOPTS","SEDGEMNLATSOPTS"]
macrolist = macrolist + ["ROKGEMNLONSOPTS","ROKGEMNLATSOPTS"]
testing = {}
testlist = ["TESTFILE","ASSUMEDGOOD_NAME"]
testlist = testlist + ["TEST_NAME","KNOWNGOOD_NAME","CHECKFLUXES"]
# to hold data from namelists.sh
# plus special for gl_config_file
namelists = {
    "gl_config_file" : "config_file"
    }
configfile = {} # to hold data from config file

# Parse command line--only one arg permitted
if len(sys.argv) != 2:
    print "usage: config2xml <config-file-name>"
    sys.exit(2)
configname = sys.argv[1]
# ensure that the arg is a genuine file
if not os.path.isfile(configname):
    print "error: could not find file %s" % configname
    sys.exit(1)

# load genie-main/namelists.sh
# we will need this to give proper names to the vars in the config file

# regex to match the lines we're interested in
p = re.compile(r'\w+="*\$\w+"*', re.IGNORECASE)
# regex to spot when we need <varref>
pp = re.compile(r'\$(\w+)\/', re.IGNORECASE)

fp = open('./namelists.sh')
while 1:
    line = fp.readline()                 # read the next line
    if not line: break                  # break loop on EOF
    line = line.strip()                 # strip '\n'
    if line == '': continue             # ignore blank lines
    if line.startswith('#'): continue   # ignore comment lines
    # try to match our namelist file regex
    m = p.findall(line)
    # store matches in a dictionary
    if m:
        for pair in m:
            [name,nickname]=pair.split('=')
            # we will access proper names, using 'nicknames'
            # found in ASCII config files
            # need to lose the leading '$' off nickname
            nickname=nickname.replace('$','')
            # and also the '"'s
            nickname=nickname.replace(r'"','')
            namelists[nickname]=name
fp.close()

# load the config file
fp = open(configname)
while 1:
    line = fp.readline()                # read the next line
    if not line: break                  # break loop on EOF
    line = line.strip()                 # strip '\n'
    if line == '': continue             # ignore blank lines
    if line.startswith('#'): continue   # ignore comment lines
    [nickname,val]=line.split('=',1)
    configfile[nickname]=val
fp.close()

# sift through contents of configfile for useful stuff..
for item in configfile:
    # find the vars in the config file and store in dict
    for var in varlist:
        if var == item:
            vars[item] = configfile[item]
    # find the makeargs
    for makearg in makearglist:
        if makearg == item:
            makeargs[item] = configfile[item]
    # find the macros
    for macro in macrolist:
        if macro == item:
            macros[item] = configfile[item]
    # testing vars
    for var in testlist:
        if var == item:
            testing[item] = configfile[item]
    
    # sort out the <control> and <config> sections
    if (item.startswith("ma_") or item.startswith("gem_") or item.startswith("gm_")):
        done=0
        for paramarray in [atm_select,ocn_select,sed_select]:
            if item in paramarray.keys():
                paramarray[item][1] = configfile[item]
                deletelist.append(item)
                done=1
        if done == 0:
            if (item.startswith("ma_flag_") and not (item.startswith("ma_flag_checkfluxes") or item.startswith("ma_flag_glim"))):
                configflags[item][1] = configfile[item]
            else:
                control[item] = configfile[item]
    # look for members of a paramarray -- finite number of special cases
    for paramarray in [snolook,watertransport,diffamp,betaz,betam,diff,albsnf_nvg,albsnf_max,hcon_nvg,hcap_nvg,atm_init,ocn_init]:
        if item in paramarray.keys():
            paramarray[item][1] = configfile[item]
            deletelist.append(item)

# must delete param array entries outside of loop
for item in deletelist:
    del configfile[item]

#print betaz
#print deletelist
#sys.exit(0)

# Let's start printing!
# header
print '<?xml version="1.0" encoding="UTF-8"?>'
print '<job author="config2xml.py - automatic conversion of ASCII text config file">'
# <vars> section
print '\t<vars>'
for var in vars:
    print '\t\t<var name="' + str(var) + '">' + str(vars[var]) + '</var>' 
print '\t</vars>'
# <config> section
print '\t<config>'
for model in configflags:
    if (configflags[model][1] == ".TRUE." or configflags[model][1] == ".true."):
        print '\t\t<model name="' + str(configflags[model][0]) + '"/>'
print '\t</config>'
# <parameters> section
print '\t<parameters>'
#    <control> subsection
print '\t\t<control>'
for entry in control:
    myval = control[entry]
    repval = pp.sub(r'<varref>\1</varref>/', myval)
    print '\t\t\t<param name="' +str(namelists[entry]) + '">' + str(repval) +'</param>'
        # --gem--
vals = ""
for list in atm_select.values():
    vals = vals + list[1]
if vals:
    print '\t\t\t<paramArray name="atm_select">'
    for key in atm_select.keys():
        if atm_select[key][1] != "":
            print '\t\t\t\t<param index="'+str(atm_select[key][0])+'">'+str(atm_select[key][1])+'</param>'
    print '\t\t\t</paramArray>'
vals = ""
for list in ocn_select.values():
    vals = vals + list[1]
if vals:
    print '\t\t\t<paramArray name="ocn_select">'
    for key in ocn_select.keys():
        if ocn_select[key][1] != "":
            print '\t\t\t\t<param index="'+str(ocn_select[key][0])+'">'+str(ocn_select[key][1])+'</param>'
    print '\t\t\t</paramArray>'
vals = ""
for list in sed_select.values():
    vals = vals + list[1]
if vals:
    print '\t\t\t<paramArray name="sed_select">'
    for key in sed_select.keys():
        if sed_select[key][1] != "":
            print '\t\t\t\t<param index="'+str(sed_select[key][0])+'">'+str(sed_select[key][1])+'</param>'
    print '\t\t\t</paramArray>'
print '\t\t</control>'
for model in configflags:
    # only need output params for selected models..
    if (configflags[model][1] == ".TRUE." or configflags[model][1] == ".true."):
        print '\t\t<model name="' + str(configflags[model][0]) + '">'
        for entry in configfile:
            if entry.startswith(configflags[model][2]+'_'):
                myval = configfile[entry]
                repval = pp.sub(r'<varref>\1</varref>/', myval)
                print '\t\t\t<param name="' +str(namelists[entry]) + '">' + str(repval) +'</param>'
        # plus paramarray specials
        # --atchem--
        if configflags[model][0] == "atchem":
            vals = ""
            for list in atm_init.values():
                vals = vals + list[1]
            if vals:
                print '\t\t\t<paramArray name="atm_init">'
                for key in atm_init.keys():
                    if atm_init[key][1] != "":
                        print '\t\t\t\t<param index="'+str(atm_init[key][0])+'">'+str(atm_init[key][1])+'</param>'
                print '\t\t\t</paramArray>'
        # --biogem--
        if configflags[model][0] == "biogem":
            vals = ""
            for list in ocn_init.values():
                vals = vals + list[1]
            if vals:
                print '\t\t\t<paramArray name="ocn_init">'
                for key in ocn_init.keys():
                    if ocn_init[key][1] != "":
                        print '\t\t\t\t<param index="'+str(ocn_init[key][0])+'">'+str(ocn_init[key][1])+'</param>'
                print '\t\t\t</paramArray>'
        # --embm--
        if configflags[model][0] == "embm":
            vals = ""
            for list in diffamp.values():
                vals = vals + list[1]
            if vals:
                print '\t\t\t<paramArray name="diffamp">'
                for key in diffamp.keys():
                    if diffamp[key][1] != "":
                        print '\t\t\t\t<param index="'+str(diffamp[key][0])+'">'+str(diffamp[key][1])+'</param>'
                print '\t\t\t</paramArray>'
            vals = ""
            for list in betaz.values():
                vals = vals + list[1]
            if vals:
                print '\t\t\t<paramArray name="betaz">'
                for key in betaz.keys():
                    if betaz[key][1] != "":
                        print '\t\t\t\t<param index="'+str(betaz[key][0])+'">'+str(betaz[key][1])+'</param>'
                print '\t\t\t</paramArray>'
            vals = ""
            for list in betam.values():
                vals = vals + list[1]
            if vals:
                print '\t\t\t<paramArray name="betam">'
                for key in betam.keys():
                    if betam[key][1] != "":
                        print '\t\t\t\t<param index="'+str(betam[key][0])+'">'+str(betam[key][1])+'</param>'
                print '\t\t\t</paramArray>'
        # --goldstein--        
        if configflags[model][0] == "goldstein":
            vals = ""
            for list in diff.values():
                vals = vals + list[1]
            if vals:
                print '\t\t\t<paramArray name="diff">'
                for key in diff.keys():
                    if diff[key][1] != "":
                        print '\t\t\t\t<param index="'+str(diff[key][0])+'">'+str(diff[key][1])+'</param>'
                print '\t\t\t</paramArray>'
        print '\t\t</model>'
print '\t</parameters>'
# <build> section
print '\t<build>'
# make args first
for arg in makeargs:
    print '\t\t<make-arg name="'+str(arg)+'">'+str(makeargs[arg])+'</make-arg>'
# then macros
for macro in macros:
    print '\t\t<macro handle="'+str(macro)+'" status="defined">'
    myval = macros[macro]
    myval = myval.replace('$(DEFINE)','')
    myval = myval.replace("'","")
    # may not contain an '='
    splitret = myval.split('=')
    ident = splitret[0]
    print '\t\t\t<identifier>'+str(ident)+'</identifier>'
    if len(splitret) > 1:
        repl = splitret[1]
        print '\t\t\t<replacement>'+str(repl)+'</replacement>'
    print '\t\t</macro>'
print '\t</build>'
# <testing> section
print '\t<testing>'
for var in testing:
    if testing[var] != "${EXPID}_assumedgood":
        print '\t\t<var name="'+str(var)+'">'+str(testing[var])+'</var>'
print '\t</testing>'
print '</job>'
