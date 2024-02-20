
load "present_group.m";
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


load "table.m";


wild_seq:=[];
exp_wild_seq:=[];



pre_calculated_table_length:=#table;
printf "\n pre_calculated_table_length=%o",pre_calculated_table_length;

w_dis:=0; 


wild_initiations:=0;

sum_wild_walk_length:=0;



ii:=0;

load "fun.m";
d1:=dis_pt_length;tmp_itr_gs:=2^d1;
printf "\n tmp_itr_gs=%o",tmp_itr_gs;
pow_2_pre_calculated_table_length:= Z!Log(2,pre_calculated_table_length);

printf "\n time for %o wild points in GS is",length_for_tame_or_wild;
time while((#wild_seq lt length_for_tame_or_wild) ) do
	y1:=Random(wl,wb);
	y2:=Random(wl,wb);
	z:=Zp!(h*((g1^y1)*(g2^y2)));
	z_in_integers:=Z!z;
	distinguished_point:=ModByPowerOf2(z_in_integers,dis_pt_length) eq 0;	
	wild_initiations:=wild_initiations+1;
	length_of_wild_walk_now:=0;
	event_wild:=true;
	ii:=ii+1;
	while((not distinguished_point) and event_wild) do
		random_elm_index:=ModByPowerOf2(z_in_integers,pow_2_pre_calculated_table_length)+1;
		ii:=ii+1;
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
		w_dis:=w_dis+1;
		Append(~wild_seq,z);
		Append(~exp_wild_seq,[y1,y2]);	
		//break;
	end if;
end while;

printf "\n points_walked_through=%o, points/(2^dis_pt_length)=%o",ii,RealField(5)!(ii/(2^dis_pt_length));

printf "\n points/(wild_length*tmp_itr_gs)=%o",RealField(5)!(ii/(wild_length*tmp_itr_gs));

printf "\n w_dis=%o, wild_initiations=%o, sum_wild_walk_length=%o",w_dis,wild_initiations,sum_wild_walk_length;
printf "\n number of distinct elements in wild seq=%o",#SequenceToSet(wild_seq);
