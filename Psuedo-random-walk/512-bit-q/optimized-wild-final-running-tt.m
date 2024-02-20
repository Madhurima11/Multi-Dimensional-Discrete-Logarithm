load "present_group.m";
printf "\n r1=%o, l1=%o",r1,l1;

printf "\n q=%o",q;
load "target_bases.m";
printf "\n g1=%o",g1;
printf "\n g2=%o",g2;
printf "\n h=%o",h;
printf "\n dis_pt_length=%o",dis_pt_length;


N:=q;

alpha:=0.1;
wl:=Ceiling(-(alpha*N)/2);
wb:=Floor((alpha*N)/2);

load "multiplier_table.m";
total_pre_calculated:=#mult_seq;
printf "\n total_pre_calculated=%o",total_pre_calculated;

load "exponent_seq_for_first_r_terms.m";
load "sizes.m";

load "fun.m";
load "m_hat_table.m";

wild_seq:=[];
exp_wild_seq:=[];



w_dis:=0; 
wild_initiations:=0;

sum_wild_walk_length:=0;


//ii:=0;


d2:=dis_pt_length_tt;num:=r*(r^d2-1); den:=(r-1); tmp_itr_tt:=r_bar+(num/den);
printf "\n tmp_itr_tt=%o",tmp_itr_tt;


points_walked_through:=0;
for_dp_seq_s_bar:=[0: i in [1..dis_pt_length_tt]];
printf "\n time for %o wild points in TT is",length_for_tame_or_wild;
time while((#wild_seq lt length_for_tame_or_wild) ) do
	y1:=Random(wl,wb);
	y2:=Random(wl,wb);
	z:=Zp!(h*((g1^y1)*(g2^y2)));
	wild_initiations:=wild_initiations+1;	
	event_wild:=true;	
	dp:=false;
	seq_s_bar:=[99: i in [1..dis_pt_length_tt]];
	length_of_wild_walk_now:=0;
	while((not dp) and ((event_wild))) do
		u_bits_z:=u_ary_bits(z);
		pos:=total_pre_calculated;
		m:=Zp!mult_seq[pos];
		points_walked_through:=points_walked_through+1;
		present_m_m_hat_seq:=m_hat_seq[pos];
		value_s_bar:=s_bar(u_bits_z,present_m_m_hat_seq);	
		seq_s_bar:=[seq_s_bar[i] : i in [2..#seq_s_bar]] cat [value_s_bar];
		dp:=seq_s_bar eq for_dp_seq_s_bar;
		j:=0;
		index_from_mult_seq:=[];
		while((j lt l1) and (not dp) and (value_s_bar ne -9999999999999) and (event_wild)) do
			random_index_from_table:=value_s_bar +1;
			Append(~index_from_mult_seq,random_index_from_table);
			pos:=position(arrange_in_lex(index_from_mult_seq));
			present_m_m_hat_seq:=m_hat_seq[pos];
			exp_now:=exponent_seq_for_first_r_terms[pos];
			tmp_y1:=y1+exp_now[1];
			tmp_y2:=y2+exp_now[2];
			event_1:=(wl le tmp_y1) and (tmp_y1 le wb); 
			event_2:=(wl le tmp_y2) and (tmp_y2 le wb);  
			event_wild:=event_1 and event_2;
			if(not event_wild) then
				sum_wild_walk_length:=sum_wild_walk_length+length_of_wild_walk_now;
				break;
			end if;
			value_s_bar:=s_bar(u_bits_z,present_m_m_hat_seq);
			if((value_s_bar eq -9999999999999)) then
				break;
			end if;
			m:=Zp!mult_seq[pos];
			
			m_y1:=exp_now[1];
			m_y2:=exp_now[2];			
			seq_s_bar:=[seq_s_bar[i] : i in [2..#seq_s_bar]] cat [value_s_bar];
			dp:=seq_s_bar eq for_dp_seq_s_bar;
			length_of_wild_walk_now:=length_of_wild_walk_now+1;
			j:=j+1;
		end while;
		points_walked_through:=points_walked_through+j;
		if((value_s_bar eq -9999999999999)) then
			break;
		end if;
		z:=Zp!prod_mod_p(z,m);	
		y1:=tmp_y1;
		y2:=tmp_y2;			
		if(dp) then
			w_dis:=w_dis+1;
			Append(~wild_seq,z);
			Append(~exp_wild_seq,[y1,y2]);
			dp:=true;
			seq_s_bar:=[99: i in [1..dis_pt_length_tt]];
		end if;
	end while;	
end while;


printf "\n points_walked_through=%o, points_walked_through/(tmp_itr_tt)=%o",points_walked_through,RealField(5)!(points_walked_through/(tmp_itr_tt));
printf "\n points_walked_through/(wild_length*tmp_itr_tt)=%o",RealField(5)!(points_walked_through/(wild_length*tmp_itr_tt));
printf "\n w_dis=%o, wild_initiations=%o, sum_wild_walk_length=%o",w_dis,wild_initiations,sum_wild_walk_length;
printf "\n number of distinct elements in wild seq=%o",#SequenceToSet(wild_seq);
