#!/usr/bin/awk -f
#
# extract violations from haddock print_xpcs.out files
#
# ana_xpcs_viol.awk format=<> print_xpcs.out 
#
# options:
#   format=cyana [default] : write cyana .pcs style file with pcs_obs
#                            and pcs_calc pcs_delta in comments
#   format=cyanapred       : same but swith pcs_obs and pcs_calc
#                            
#   format=anythingelse    : free format table+header                      
#
#
BEGIN {
    stderr="/dev/stderr"

    read_pcs_violation=0
    fviols=0
          
    format = "cyana"
    weight="1.0"
    
} 

FNR == 1 { 
    
    if (format=="cyana") { 
        print "#  . . Atom         PCS   Error    Weight Sample",

            "# pcs_calc pcs_delta Class pdb"
    } else if (format=="cyanapred") { 
        print "#  . . Atom         PCS_calc Error    Weight Sample",
            "# pcs_obs pcs_delta Class pdb"
    } else {
        print "# Trn Tres_ch Tat",
              "rn res_ch pcs_atom",
              "pcs_calc +/- pcs_error pcs_obs pcs_delta # class pdb"
    }
}

/The following pseudocontact shifts have/ { 
    # start of violation section
    read_pcs_violation=1
    
}

/ASSFIL:/{  pdb=$3 }   


/ FVIOLS = / { 

    fviols = $3   # end of violation section
    
    read_pcs_violation=0   # end of violation section

    pcs_class = ""
    Tchain="" ; Tresnum="" ; Tresnam="" ; Tatnam = ""
    Nchain="" ; Nresnum="" ; Nresnam="" ; Natnam = ""
    pcs_calc="";pcs_obs="";pcs_error="";pcs_delta=""

} 
    
read_pcs_violation { 
    if (/Class/) pcs_class = $2
    if (/==============/) { 
        ## start (or end) of single violation
        Tchain="" ; Tresnum="" ; Tresnam="" ; Tatnam = ""
        Nchain="" ; Nresnum="" ; Nresnam="" ; Natnam = ""
        pcs_calc="";pcs_obs="";pcs_error="";pcs_delta=""
    }
    if (NF==0) { 
        ## end of violations section
    }
    

    if (/Set-M-atoms/) { 
        getline 
        Tchain=$1 ; Tresnum=$2 ; Tresnam=$3 ; Tatnam = $4
        next
    }
    
    if (/Set-N-atoms/) { 
        getline 
        Nchain=$1 ; Nresnum=$2 ; Nresnam=$3 ; Natnam = $4
        next
    }
    
    if (/Calc:/&&/Obs:/) {
        pcs_calc=$2 ; pcs_obs=$4
        next
    }
    
    if (/Error:/&&/Delta:/) {
        pcs_error=$2 ; pcs_delta=$4
        
        if (format == "cyana") { 
            print sprintf(" %3i %-4s %-4s",Nresnum,Nresnam,Natnam),
                  sprintf(" %7.3f %7.3f %9.1f",pcs_obs,pcs_error,weight),
                  sprintf("%6s",Tresnum),
                  sprintf("# %8.3f  %8.3f %5s %s",pcs_calc,pcs_delta,pcs_class,pdb)
        } else if (format == "cyanapred") { 
            print " ",
                  Nresnum,Nresnam,Natnam,
                  pcs_calc,pcs_error,weight,
                  Tresnum,
                  "#",pcs_obs,pcs_delta,
                  pcs_class,pdb
           } else  { 
            print " ",
                  Nresnum,Nresnam"_"Nchain,Natnam,
                  Tresnum,Tresnam"_"Tchain,Tatnam,
                  pcs_obs,"+/-",pcs_error,pcs_calc,pcs_delta,"#",
                  pcs_class,pdb
           } 
 
    }
}
