
load "present_group.m";
load "target_bases.m";
N:=q;
m:=N div (10*Floor(Log(2,N)));

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
printf "\n finding repeatations if any";
printf "\n #table_as_seq=%o, #SequenceToSet(table_as_seq)=%o",#table_as_seq,#SequenceToSet(table_as_seq);
PrintFile("table.m","table:=":Overwrite:=true);                              
PrintFile("table.m",table_as_seq);             
PrintFile("table.m",";");  
