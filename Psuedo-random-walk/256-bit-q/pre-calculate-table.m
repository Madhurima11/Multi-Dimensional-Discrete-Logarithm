
load "present_group.m";
load "target_bases.m";
N:=q;
m:=N div (1000*Floor(Log(2,N)));

printf "\n m=%o",m;

table_as_seq:=[ ];
length:=1024;
for i in [1..length] do
	x:=Random(-m,m); 
	y:=Random(-m,m);
	elm:=Zp!(g1^x)*(g2^y);elm:=Z!elm;
	seq:=[elm,x,y];
	Append(~table_as_seq,[elm,x,y]);
end for;
printf "\n #pre calculated table=%o",#table_as_seq;
PrintFile("table.m","table:=":Overwrite:=true);                              
PrintFile("table.m",table_as_seq);             
PrintFile("table.m",";");  
