BEGIN{
 OFS=ORS=""
 }{
 if(FILENAME==ARGV[1]){
   Nmrk++
   Npop=NF
   for(i=1;i<=Npop;i++){DATA_N[NR,i]=$i}
   }
   
   if(FILENAME==ARGV[2]){
    N_Y++
   for(i=1;i<=Npop;i++){DATA_Y[N_Y,i]=$i}
   }
   }END{
   {print "list(I=",Nmrk,",J=",Npop,",Y=structure(\.Data=c(">"WINBUGS_DATA"}
   for(i=1;i<=Nmrk-1;i++){
     for(j=1;j<=Npop;j++){print DATA_Y[i,j],",">"WINBUGS_DATA"}
     }
   for(j=1;j<=Npop-1;j++){print DATA_Y[Nmrk,j],",">"WINBUGS_DATA"}
   {print DATA_Y[Nmrk,Npop],"),.Dim=c(",Nmrk,",",Npop,")),N=structure(\.Data=c(">"WINBUGS_DATA"} 
   for(i=1;i<=Nmrk-1;i++){
     for(j=1;j<=Npop;j++){print DATA_N[i,j],",">"WINBUGS_DATA"}
     }
   for(j=1;j<=Npop-1;j++){print DATA_N[Nmrk,j],",">"WINBUGS_DATA"}   

   {print DATA_Y[Nmrk,Npop],"),.Dim=c(",Nmrk,",",Npop,")))">"WINBUGS_DATA"} 
   
   
      {print "list(c=c(">"WINBUGS_INITS"}
      for(i=1;i<=Npop-1;i++){print "0.05,">"WINBUGS_INITS"}
      {print "0.05),p=c(">"WINBUGS_INITS"}
      for(i=1;i<=Nmrk-1;i++){print "0.5,">"WINBUGS_INITS"}   
      {print "0.5))">"WINBUGS_INITS"}  
      
      }    
