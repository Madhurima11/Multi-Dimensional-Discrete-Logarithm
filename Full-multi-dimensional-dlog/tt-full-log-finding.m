load "present_group.m";
printf "\n r1=%o, l1=%o",r,l;

printf "\n q=%o",q;
load "target_bases.m";
printf "\n g1=%o",g1;
printf "\n g2=%o",g2;
printf "\n h=%o",h;
printf "\n dis_pt_length_tt=%o",dis_pt_length_tt;


N:=q;
tl:=0;
tb:=N;
alpha:=0.1;
wl:=Ceiling(-(alpha*N)/2);
wb:=Floor((alpha*N)/2);

load "multiplier_table.m";
total_pre_calculated:=#mult_seq;

load "exponent_seq_for_first_r_terms.m";
load "sizes.m";

load "fun.m";
load "m_hat_table.m";

tame_seq:=[];
exp_tame_seq:=[];


wild_seq:=[];
exp_wild_seq:=[];

t_dis:=0; 
tame_initiations:=0;


w_dis:=0; 
wild_initiations:=0;

sum_tame_walk_length:=0;
sum_wild_walk_length:=0;


//ii:=0;


d2:=dis_pt_length_tt;num:=r*(r^d2-1); den:=(r-1); tmp_itr_tt:=r_bar+(num/den);
printf "\n tmp_itr_tt=%o",tmp_itr_tt;


tame_points_walked_through:=0;
wild_points_walked_through:=0;

for_dp_seq_s_bar:=[0: i in [1..dis_pt_length_tt]];
discrete_log_solved:=0;

round:=0;
printf "\n time for full log finding in tt is";
time while(discrete_log_solved eq 0) do
	round:=round+1;
	printf "\n round=%o",round;
	this_round_tame_dis:=0;
	while(this_round_tame_dis lt length_for_tame_or_wild ) do
		x1:=Random(tl,tb);
		x2:=Random(tl,tb);
		z:=Zp!((g1^x1)*(g2^x2));//z;
		//tame_initiations:=tame_initiations+1;	
		event_tame:=true;	
		dp:=false;
		seq_s_bar:=[99: i in [1..dis_pt_length_tt]];
		length_of_tame_walk_now:=0;
		while((not dp) and ((event_tame))) do
			u_bits_z:=u_ary_bits(z);//printf "\n step 16";
			//printf "\n z=%o, u_bits_z=%o",z,u_bits_z;
			//u_bits_z;
			pos:=total_pre_calculated;
			m:=Zp!mult_seq[pos];
			tame_points_walked_through:=tame_points_walked_through+1;
			present_m_m_hat_seq:=m_hat_seq[pos];
			value_s_bar:=s_bar(u_bits_z,present_m_m_hat_seq);//value_s_bar;
			seq_s_bar:=[seq_s_bar[i] : i in [2..#seq_s_bar]] cat [value_s_bar];
			dp:=seq_s_bar eq for_dp_seq_s_bar; //printf "\n step 15";
			j:=0;
			index_from_mult_seq:=[];
			while((j lt l1) and (not dp) and (value_s_bar ne -9999999999999) and (event_tame)) do
				random_index_from_table:=value_s_bar +1;//random_index_from_table;
				//printf "\n step 1";
				Append(~index_from_mult_seq,random_index_from_table);
				//printf "\n step 2";
				pos:=position(arrange_in_lex(index_from_mult_seq));//pos;
				//printf "\n step 3";
				//printf "\n index_from_mult_seq=%o, pos=%o",index_from_mult_seq,pos;
				present_m_m_hat_seq:=m_hat_seq[pos];//present_m_m_hat_seq;
				//printf "\n step 4";
				//printf "\n pos=%o, present_m_m_hat_seq=%o",pos,present_m_m_hat_seq;
				exp_now:=exponent_seq_for_first_r_terms[pos];
				//printf "\n step 5";
				tmp_x1:=x1+exp_now[1];//printf "\n step 6";
				tmp_x2:=x2+exp_now[2];//printf "\n step 7";
				event_1:=(tl le tmp_x1) and (tmp_x1 le tb);//printf "\n step 8"; 
				event_2:=(tl le tmp_x2) and (tmp_x2 le tb);  
				//printf "\n step 9";
				event_tame:=event_1 and event_2;
				if(not event_tame) then
					sum_tame_walk_length:=sum_tame_walk_length+length_of_tame_walk_now;
					break;
				end if;
				value_s_bar:=s_bar(u_bits_z,present_m_m_hat_seq);//value_s_bar;
				//printf "\n step 10";
				if((value_s_bar eq -9999999999999)) then
					break;
				end if;
				m:=Zp!mult_seq[pos];//printf "\n inside tame loop 2 m=%o",m;
				//printf "\n only tame m=%o",m;
				//printf "\n step 11";
				
				m_x1:=exp_now[1];
				m_x2:=exp_now[2];			
				seq_s_bar:=[seq_s_bar[i] : i in [2..#seq_s_bar]] cat [value_s_bar];
				//printf "\n step 12";
				dp:=seq_s_bar eq for_dp_seq_s_bar;
				length_of_tame_walk_now:=length_of_tame_walk_now+1;
				j:=j+1;
			end while;
			tame_points_walked_through:=tame_points_walked_through+j;
			if((value_s_bar eq -9999999999999)) then
				break;
			end if;
			z:=Zp!prod_mod_p(z,m);	
			//printf "\n in tame z=%o, m=%o",z,m;
			//printf "\n step 13";
			x1:=tmp_x1;
			x2:=tmp_x2;			
			if(dp) then
				this_round_tame_dis :=this_round_tame_dis +1;
				t_dis:=t_dis+1;
				Append(~tame_seq,z);
				Append(~exp_tame_seq,[x1,x2]);
				dp:=true;
				if(z in wild_seq) then
					printf "\n SOLVED \n ";
					index:=Index(wild_seq,z);
					exponents_wild:=exp_wild_seq[index];
					y1:=exponents_wild[1];
					y2:=exponents_wild[2];
					a1:=x1-y1;
					a2:=x2-y2;
					printf "\n a1=%o, a2=%o",a1,a2;
					discrete_log_solved:=1;
					break;
				end if;
				seq_s_bar:=[99: i in [1..dis_pt_length_tt]];//printf "\n step 14";
				//printf "\n dp=%o, event_tame=%o",dp,event_tame;
			end if;
			//printf "\n step 17";
		end while;
		//printf "\n step 18";	
	end while;
	if(discrete_log_solved eq 1) then
		break;
	end if;
	//printf "\n step 19";
	this_round_wild_dis:=0;
	while(this_round_wild_dis lt length_for_tame_or_wild ) do
		y1:=Random(wl,wb);
		y2:=Random(wl,wb);
		z:=Zp!h*((g1^y1)*(g2^y2));//printf "\n step 20";
		//wild_initiations:=wild_initiations+1;	
		event_wild:=true;	
		dp:=false;
		seq_s_bar:=[99: i in [1..dis_pt_length_tt]];
		length_of_wild_walk_now:=0;
		while((not dp) and ((event_wild))  ) do
			u_bits_z:=u_ary_bits(z);//printf "\n step 21";
			pos:=total_pre_calculated;
			m:=Zp!mult_seq[pos];// printf "\n inside wild loop 2 m=%o",m;
			wild_points_walked_through:=wild_points_walked_through+1;
			present_m_m_hat_seq:=m_hat_seq[pos];
			value_s_bar:=s_bar(u_bits_z,present_m_m_hat_seq);//printf "\n step 22";	
			seq_s_bar:=[seq_s_bar[i] : i in [2..#seq_s_bar]] cat [value_s_bar];
			dp:=seq_s_bar eq for_dp_seq_s_bar;
			j:=0;
			index_from_mult_seq:=[];//printf "\n (j lt l1) and (not dp) and (value_s_bar ne -9999999999999)=%o",(j lt l1) and (not dp) and (value_s_bar ne -9999999999999);
			while((j lt l1) and (not dp) and (value_s_bar ne -9999999999999) and (event_tame)) do
				random_index_from_table:=value_s_bar +1;//random_index_from_table;
				Append(~index_from_mult_seq,random_index_from_table);
				pos:=position(arrange_in_lex(index_from_mult_seq));//printf "\n pos=%o", pos;
				present_m_m_hat_seq:=m_hat_seq[pos];//printf "\n step 23";
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
				m:=Zp!mult_seq[pos];//printf "\n only wild m=%o",m;
				
				m_y1:=exp_now[1];
				m_y2:=exp_now[2];			
				seq_s_bar:=[seq_s_bar[i] : i in [2..#seq_s_bar]] cat [value_s_bar];
				dp:=seq_s_bar eq for_dp_seq_s_bar;//printf "\n step 24";
				length_of_wild_walk_now:=length_of_wild_walk_now+1;
				j:=j+1;
			end while;
			wild_points_walked_through:=wild_points_walked_through+j;
			if((value_s_bar eq -9999999999999)) then
				break;
			end if;
			z:=Zp!prod_mod_p(z,m);	//printf "\n in wild z=%o, m=%o",z,m;
			//printf "\n in wild z=%o",z;
			//printf "\n step 25 \n ";
			y1:=tmp_y1;
			y2:=tmp_y2;			
			if(dp) then
				this_round_wild_dis :=this_round_wild_dis +1;
				w_dis:=w_dis+1;
				Append(~wild_seq,z);
				Append(~exp_wild_seq,[y1,y2]);
				dp:=true;
				if(z in tame_seq) then
					printf "\n SOLVED \n ";
					index:=Index(tame_seq,z);
					exponents_tame:=exp_tame_seq[index];
					x1:=exponents_tame[1];
					x2:=exponents_tame[2];
					a1:=x1-y1;
					a2:=x2-y2;
					printf "\n a1=%o, a2=%o",a1,a2;
					discrete_log_solved:=1;
					break;
				end if;
				seq_s_bar:=[99: i in [1..dis_pt_length_tt]];
				//printf "\n dp=%o, event_wild=%o",dp,event_wild;
			end if;
		end while;	
	end while;
end while;
printf "\n ensuring validity \n";
printf "\n the multi dimensional dlogs are a1=%o, a2=%o",a1,a2;
printf "\n h eq (g1^a1)*(g2^a2)=%o", h eq (g1^a1)*(g2^a2);
					
points_walked_through:=tame_points_walked_through+wild_points_walked_through;
printf "\n points_walked_through=%o, points_walked_through/(tmp_itr_tt)=%o",points_walked_through,RealField(5)!(points_walked_through/(tmp_itr_tt));
printf "\n tame_points_walked_through=%o, #tame_seq=%o",tame_points_walked_through,#tame_seq;

printf "\n wild_points_walked_through=%o, #wild_seq=%o",wild_points_walked_through,#wild_seq;

