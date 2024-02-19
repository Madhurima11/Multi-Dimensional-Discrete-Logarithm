load "present_group.m";
printf "\n r1=%o, l1=%o",r1,l1;

printf "\n q=%o",q;
load "target_bases.m";
printf "\n g1=%o",g1;
printf "\n g2=%o",g2;
printf "\n h=%o",h;
printf "\n dis_pt_length=%o",dis_pt_length;


N:=q;
tl:=0;
tb:=N;
alpha:=0.1;
wl:=Ceiling(-(alpha*N)/2);
wb:=Floor((alpha*N)/2);

load "table.m";


wild_seq:=[];
exp_wild_seq:=[];

tame_seq:=[];
exp_tame_seq:=[];

pre_calculated_table_length:=#table;
printf "\n pre_calculated_table_length=%o",pre_calculated_table_length;

w_dis:=0; 
wild_initiations:=0;

t_dis:=0; 
tame_initiations:=0;

sum_tame_walk_length:=0;
sum_wild_walk_length:=0;

load "fun.m";
d1:=dis_pt_length;tmp_itr_gs:=2^d1;
printf "\n tmp_itr_gs=%o",tmp_itr_gs;
pow_2_pre_calculated_table_length:= Z!Log(2,pre_calculated_table_length);


tame_points_walked_through:=0;
wild_points_walked_through:=0;

discrete_log_solved:=0;

round:=0;
printf "\n time for full log finding in gs is";
time while(discrete_log_solved eq 0) do
	round:=round+1;
	printf "\n round=%o",round;
	this_round_tame_dis:=0;
	while(this_round_tame_dis lt length_for_tame_or_wild ) do
		//printf "\n hi 1";
		x1:=Random(tl,tb);
		x2:=Random(tl,tb);
		z:=Zp!((g1^x1)*(g2^x2));//printf "\n random choice of z=%o",z;
		z_in_integers:=Z!z;
		distinguished_point:=ModByPowerOf2(z_in_integers,dis_pt_length) eq 0;	
		tame_initiations:=tame_initiations+1;
		length_of_tame_walk_now:=0;
		event_tame:=true;
		tame_points_walked_through:=tame_points_walked_through+1;//printf "\n hi 2";
		while((not distinguished_point) and event_tame) do
			random_elm_index:=ModByPowerOf2(z_in_integers,pow_2_pre_calculated_table_length)+1;//random_elm_index;
			//printf "\n random_elm_index=%o",random_elm_index;
			//printf "\nz_in_integers=%o, pow_2_pre_calculated_table_length=%o",z_in_integers,pow_2_pre_calculated_table_length; 
			tame_points_walked_through:=tame_points_walked_through+1;
			random_elm_from_table:=table[random_elm_index];
			m:=random_elm_from_table[1]; //printf "\n m=%o",m;
			z:=Zp!prod_mod_p(z,m);	z_in_integers:=Z!z;//printf "\n z=%o",z;
			x1:=x1+random_elm_from_table[2];
			x2:=x2+random_elm_from_table[3];			
		
			event_1:=(tl le x1) and (x1 le tb); 
			event_2:=(tl le x2) and (x2 le tb); 
			event_tame:=event_1 and event_2;
			length_of_tame_walk_now:=length_of_tame_walk_now+1;
			if(not event_tame) then
				sum_tame_walk_length:=sum_tame_walk_length+length_of_tame_walk_now;
				break;	
			else
				distinguished_point:=ModByPowerOf2(z_in_integers,dis_pt_length) eq 0;	
			end if;	
			
		end while;
		if(distinguished_point) then
			this_round_tame_dis :=this_round_tame_dis +1;
			t_dis:=t_dis+1;
			Append(~tame_seq,z);
			Append(~exp_tame_seq,[x1,x2]);	
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
		end if;
	end while;
	if(discrete_log_solved eq 1) then
		break;
	end if;
	//printf "\n hi 3";
	this_round_wild_dis:=0;
	while(this_round_wild_dis lt length_for_tame_or_wild ) do
		y1:=Random(wl,wb);
		y2:=Random(wl,wb);
		z:=Zp!(h*((g1^y1)*(g2^y2)));
		z_in_integers:=Z!z;
		distinguished_point:=ModByPowerOf2(z_in_integers,dis_pt_length) eq 0;	
		wild_initiations:=wild_initiations+1;
		length_of_wild_walk_now:=0;
		event_wild:=true;
		wild_points_walked_through:=wild_points_walked_through+1;
		while((not distinguished_point) and event_wild) do
			random_elm_index:=ModByPowerOf2(z_in_integers,pow_2_pre_calculated_table_length)+1;
			wild_points_walked_through:=wild_points_walked_through+1;
			random_elm_from_table:=table[random_elm_index];
			m:=random_elm_from_table[1];
			z:=Zp!prod_mod_p(z,m);	z_in_integers:=Z!z;
			y1:=y1+random_elm_from_table[2];
			y2:=y2+random_elm_from_table[3];			
		
			event_1:=(wl le y1) and (y1 le wb); 
			event_2:=(wl le y2) and (y2 le wb); 
			event_wild:=event_1 and event_2;
			length_of_wild_walk_now:=length_of_wild_walk_now+1;
			if(not event_wild) then
				sum_wild_walk_length:=sum_wild_walk_length+length_of_wild_walk_now;
				break;	
			else
				distinguished_point:=ModByPowerOf2(z_in_integers,dis_pt_length) eq 0;	
			end if;	
			
		end while;
		if(distinguished_point) then
			this_round_wild_dis :=this_round_wild_dis +1;
			w_dis:=w_dis+1;
			Append(~wild_seq,z);
			Append(~exp_wild_seq,[y1,y2]);	
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
		end if;	
	end while;//printf "\n hi 4";
end while;

printf "\n ensuring validity \n";
printf "\n the multi dimensional dlogs are a1=%o, a2=%o",a1,a2;
					printf "\n h eq (g1^a1)*(g2^a2)=%o", h eq (g1^a1)*(g2^a2);
points_walked_through:=tame_points_walked_through+wild_points_walked_through;

printf "\n points_walked_through=%o, points_walked_through/(tmp_itr_gs)=%o",points_walked_through,RealField(5)!(points_walked_through/(tmp_itr_gs));
printf "\n tame_points_walked_through=%o, #tame_seq=%o",tame_points_walked_through,#tame_seq;

printf "\n wild_points_walked_through=%o, #wild_seq=%o",wild_points_walked_through,#wild_seq;

