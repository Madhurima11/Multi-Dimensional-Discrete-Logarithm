load "present_group.m";
g1:=Random(Gq) @ iso_Gq_to_Zq @ iso_unit_group_to_Zp;
g2:=Random(Gq) @ iso_Gq_to_Zq @ iso_unit_group_to_Zp;



h:=(g1^Random(0,bound))*(g2^Random(0,bound));


PrintFile("target_bases.m","g1:=":Overwrite:=true);                              
PrintFile("target_bases.m",g1);             
PrintFile("target_bases.m",";");  

PrintFile("target_bases.m","g2:=");                              
PrintFile("target_bases.m",g2);             
PrintFile("target_bases.m",";");  

PrintFile("target_bases.m","h:=");                              
PrintFile("target_bases.m",h);             
PrintFile("target_bases.m",";");

PrintFile("target_bases.m","g1:=Zp!g1;"); 
PrintFile("target_bases.m","g2:=Zp!g2;"); 
PrintFile("target_bases.m","h:=Zp!h;"); 

