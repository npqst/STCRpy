BEGIN { FS = " "}
{ 
   if ($1 == "set-i-atoms") 
    { istop = 0;
      printf "(";
      while ( istop == 0)
       { getline inew ;
         split (inew,id);
	 if (id[1] != "set-j-atoms") {printf "%5s %5s %5s %5s",id[1],id[2],id[3],id[4]}
	 if (id[1] == "set-j-atoms") 
	  {istop = 1;
           jstop = 0;
	   printf ") (";
           while ( jstop == 0)
            { getline jnew;
              split (jnew,jd);
	      if (jd[1] != "R<average>=") {printf "%5s %5s %5s %5s",jd[1],jd[2],jd[3],jd[4]}
	      if (jd[1] == "R<average>=") 
	       {jstop = 1; 
	        printf ") Rexp= %7.3f Rave= %7.3f Viol= %7.3f \n",jd[4]+jd[7],jd[2],jd[9]; 
		break
	       }
            }
	  }
       }
    }
} 
