load "sizes.m";
pow_2_in_w:=Floor(Log(2,w));
function u_ary_bits(a)
	q:=Z!a;
	bit_seq:=[];
	while(q ne 0) do
		r:=q mod u;
		q:=Floor(q/u);
		Append(~bit_seq,r);
	end while;	
	return (bit_seq cat [0: i in [1..(d-#bit_seq)]]);
end function;


function get_m_hat_seq(any_m) 
	arb_m_hat_seq:=[];
	for j in [1..d] do
		num:=(Z!(u^(j-1))*Z!any_m) mod p;
		den:=w_bar;
		Append(~arb_m_hat_seq,Floor(num/den));
	end for;
	return arb_m_hat_seq;
end function;


function add_binary(b1,b2)// assume b1 and b2 are of same size
	add_seq:=[];
	carry:=0;
	for j in [1..#b1] do
		i:=#b1-j+1;
		tmp_sum:=b1[i]+b2[i]+carry;
		if((tmp_sum) lt 2) then
			carry:=0;
		else
			carry:=tmp_sum div 2;
			tmp_sum:=tmp_sum mod 2;	
		end if;
		Append(~add_seq,tmp_sum);
	end for;
	if(carry eq 1) then
		Append(~add_seq,carry);
	end if;
	return Reverse(add_seq);
end function;


function multiply_bin_seq_1_one_element_of_2(seq_of_1, elm_of_2)
	if(elm_of_2 eq 0) then
		return [0];
	else
		return seq_of_1;	
	end if;	
end function;

function multiply(bin_seq_1,bin_seq_2) 
	sum_seq:=[0: i in [1..#bin_seq_1]];
	for i in [1..#bin_seq_2] do
		j:=#bin_seq_2-i+1;
		elm2:=bin_seq_2[j];
		tmp_mult_seq:=multiply_bin_seq_1_one_element_of_2(bin_seq_1, elm2) cat [0: ind in [1..(i-1)]];
		m:=Max(#tmp_mult_seq,#sum_seq);	tmp_mult_seq:=[0: i in [1..(m-#tmp_mult_seq)]] cat 	tmp_mult_seq;
		
sum_seq:=[0: i in [1..(m-#sum_seq)]] cat sum_seq;
		sum_seq:=add_binary(sum_seq,tmp_mult_seq);
		//sum_seq;
	end for;
	return sum_seq;
end function;


function prod_mod_p(z,m)
	bin_z:=Reverse(Intseq(Z!z,2));bin_m:=Reverse(Intseq(Z!m,2));
	mult_prod_in_binary:=multiply(bin_z,bin_m);
	sum:=0; 
	t:=mult_prod_in_binary;
	for ind in [1..#t] do   
		j:=#t-ind+1;            
		sum:=sum+(t[j]*(2^(ind-1)));
	end for;
	return (sum mod p);
end function;

function part_tau_bar(u_bits_z2,m2_m_hat_seq)
	sum_part_tau:=0;	
	for i in [1..d] do
		i_th_u_bit_z:=u_bits_z2[i];
		i_th_m_hat:=m2_m_hat_seq[i];
		bin_z:=Reverse(Intseq(Z!i_th_u_bit_z,2));
		bin_m:=Reverse(Intseq(Z!i_th_m_hat,2)); 
		mult_prod_in_binary:=multiply(bin_z,bin_m);
		sum:=0; 
		t:=mult_prod_in_binary;
		for ind in [1..#t] do   
			j:=#t-ind+1;            
			sum:=sum+(t[j]*(2^(ind-1)));
		end for;
		sum_part_tau:=sum_part_tau+sum;
	end for;
	prod:=ModByPowerOf2(sum,pow_2_in_w);
	return prod;
end function;


function s_bar(u_bits_z2,m2_m_hat_seq)
	value_tau_bar:=Floor(part_tau_bar(u_bits_z2,m2_m_hat_seq)/t_bar);
	if(((value_tau_bar+1) mod r_bar) eq 0)then
		return (-9999999999999);
	else
		return Floor(value_tau_bar/r_bar);
	end if;
end function;

function choose(a,b)
d:=Factorial(b)*Factorial(a-b);
return Factorial(a) div d;
end function;

function g(a,b)
return choose(r-a+b-1,b-1);
end function;


function position(seq)
	k:=#seq;
	sum:=0;  
	for j in [2..k] do                             
		if((seq[j]+1) le r) then
			sum:=sum+&+([g(i,k-j+1):i in [seq[j]+1..r]]);
		end if;
	end for;
	return choose(r+k-1,k-1)+(&+[g(i,k):i in [1..seq[1]]])-sum-1;
end function;


function arrange_in_lex(seq)
	total:=#seq;
	tmp:=[];
	while(#tmp ne total) do
		elm, pos:=Min(seq);
		Append(~tmp,elm );
		Remove(~seq,pos);		
	end while;
	return tmp;
end function;





function next_term(seq)
	
	if(seq[#seq] lt r1) then
		seq[#seq]:=seq[#seq]+1;
	else
		i:=#seq;
		while(seq[i] eq r1) do
			i:=i-1;
			if(i eq 0) then
				break;
			end if;
		end while;
		if(i eq 0) then
			return [1: j in [1..(#seq+1)]];
		else
			tmp:=seq[i]+1;
			for j in [i..#seq] do
				seq[j]:=tmp;
			end for;	
		end if;
	end if;		
	return seq;
end function;

