BEGIN { FS = " "}
{ 
  icont = 0;
  inum  = 1;
  iave = NF - 2;
  iexp = NF - 4;
  iviol = NF;
  rexp = $iexp;
  ave = $iave;
  viol = $iviol;
  split ($0 ,istore, "Rexp=");
  while (icont == 0)
  {
   lastline = FNR;
   getline new;
   if (FNR == lastline) {icont == 1; exit};
   split (new, jstore, "Rexp=");
   split (new, tmp, " ");
   if (istore[1] == jstore[1])
     {inum = inum + 1;
      ave = ave + tmp[iave];
      viol = viol + tmp[iviol];
     } 
   else if (istore[1] != jstore[1])
     { ave = ave / inum;
       viol = viol / inum;
       printf "Rexp=%8.3f Rave=%8.3f Viol=%8.3f #viol=%5d %s\n",rexp,ave,viol,inum,istore[1];
       istore[1] = jstore [1];
       split (new, tmp, " ");
       if (tmp[1] == "DONE") {exit};
       found = 0;
       ii = 1;
       while (found != 1)
       { if (tmp[ii] == "Rexp=") {found = 1; nf = ii} 
         ii = ii + 1;
       }
       iave = nf + 3;
       iexp = nf + 1;
       iviol = nf + 5;
       rexp = tmp[iexp];
       ave = tmp[iave];
       viol = tmp[iviol];
       inum = 1;
   }
 }
} 
