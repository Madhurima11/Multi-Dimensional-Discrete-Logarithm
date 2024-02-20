load "present_group.m";
load "target_bases.m";
printf "\n r1=%o, l1=%o",r1,l1;

N:=q;
den_m_tt:=l1;
m1:=N div (1000*Floor(Log(2,N)));
m:=Floor(m1/den_m_tt);
printf "\n m=%o",m;
load "sizes.m";
load "fun.m";
mult_seq:=[];
tuple_parent:=car<Zp,Z,Z>;
mult_seq:=[];
exponent_seq_for_first_r_terms:=[];
seq_m_hat:=[];
ml:=m;
mb:=m;

for i in [1..r] do
	x:=Random(-ml,mb);
	y:=Random(-ml,mb);
	elm:=Zp!(g1^x)*(g2^y);
	Append(~mult_seq,elm);
	Append(~exponent_seq_for_first_r_terms,[x,y]);
	Append(~seq_m_hat,get_m_hat_seq(elm)); 
end for;

term:=[r1];
for i in [(r1+1)..(choose(r1+l1,r1)-1)] do
	term:=next_term(term);
	x:=&+[exponent_seq_for_first_r_terms[term[i]][1] : i in [1..#term]];
	y:=&+[exponent_seq_for_first_r_terms[term[i]][2] : i in [1..#term]];
	ind_1:=term[1];
	ind_2:=position(Exclude(term,term[1]));
	elm:=mult_seq[ind_1]*mult_seq[ind_2];
	Append(~mult_seq,elm);
	Append(~seq_m_hat,get_m_hat_seq(elm)); 
	Append(~exponent_seq_for_first_r_terms,[x,y]);
end for;
	
elm:=Zp!1;
Append(~mult_seq,elm);
Append(~seq_m_hat,get_m_hat_seq(elm)); 
Append(~exponent_seq_for_first_r_terms,[0,0]);

printf "\n #mult_seq=%o",#mult_seq;
PrintFile("multiplier_table.m","mult_seq:=":Overwrite:=true);                              
PrintFile("multiplier_table.m",mult_seq);             
PrintFile("multiplier_table.m",";");  

PrintFile("exponent_seq_for_first_r_terms.m","exponent_seq_for_first_r_terms:=":Overwrite:=true);                              
PrintFile("exponent_seq_for_first_r_terms.m",exponent_seq_for_first_r_terms);             
PrintFile("exponent_seq_for_first_r_terms.m",";");  


PrintFile("m_hat_table.m","m_hat_seq:=":Overwrite:=true);                              
PrintFile("m_hat_table.m",seq_m_hat);             
PrintFile("m_hat_table.m",";");  
